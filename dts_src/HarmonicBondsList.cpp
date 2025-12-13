#include "HarmonicBondsList.h"
#include "State.h"
#include <algorithm>
#include <cctype>

// E = -k*(Field*Inc_direction)^2 --> we could make it k1*(Field*Inc_direction)-k2*(Field*Inc_direction)^2

HarmonicBondsList::HarmonicBondsList(State *pState, std::string filename) {
    m_FileName = filename;
    m_pState = pState;
}
HarmonicBondsList::~HarmonicBondsList()
{
    
}
void HarmonicBondsList::Initialize() {
    // Open the bond file
    std::ifstream input(m_FileName);
    if (!input) {
        std::cerr << "Error: Failed to open the input file: " << m_FileName << std::endl;
        exit(EXIT_FAILURE);
    }

    // Get active vertices
    const std::vector<vertex*>& activeVertices = m_pState->GetMesh()->GetActiveV();
    const size_t numVertices = activeVertices.size();

    // Default bond parameters (can be overridden per bond or from file header)
    double default_k = 10.0;
    double default_l0 = 0.472;

    std::string line;
    int bondId = 0;

    // Read the file line by line
    while (std::getline(input, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == ';') {
            // Try to parse default values from comments
            // Supports formats like: "; k=100 l0=0.5" or "; k_stretch = 1000" or "; l0 = 0.472"
            std::string lower_line = line;
            std::transform(lower_line.begin(), lower_line.end(), lower_line.begin(), ::tolower);
            
            // Look for k_stretch or k = value
            size_t k_pos = lower_line.find("k_stretch");
            if (k_pos == std::string::npos) {
                k_pos = lower_line.find(" k =");
                if (k_pos == std::string::npos) {
                    k_pos = lower_line.find(" k=");
                }
            }
            if (k_pos != std::string::npos) {
                size_t eq_pos = lower_line.find("=", k_pos);
                if (eq_pos != std::string::npos) {
                    std::string val_str = line.substr(eq_pos + 1);
                    // Remove any trailing comments or whitespace
                    size_t comment_pos = val_str.find(";");
                    if (comment_pos != std::string::npos) {
                        val_str = val_str.substr(0, comment_pos);
                    }
                    std::vector<std::string> val_tokens = Nfunction::Split(val_str);
                    if (!val_tokens.empty()) {
                        default_k = Nfunction::String_to_Double(val_tokens[0]);
                    }
                }
            }
            
            // Look for l0 = value
            size_t l0_pos = lower_line.find("l0");
            if (l0_pos != std::string::npos) {
                size_t eq_pos = lower_line.find("=", l0_pos);
                if (eq_pos != std::string::npos) {
                    std::string val_str = line.substr(eq_pos + 1);
                    // Remove any trailing comments or whitespace
                    size_t comment_pos = val_str.find(";");
                    if (comment_pos != std::string::npos) {
                        val_str = val_str.substr(0, comment_pos);
                    }
                    std::vector<std::string> val_tokens = Nfunction::Split(val_str);
                    if (!val_tokens.empty()) {
                        default_l0 = Nfunction::String_to_Double(val_tokens[0]);
                    }
                }
            }
            continue;
        }

        // Split the line into components
        std::vector<std::string> tokens = Nfunction::Split(line);
        if (tokens.size() < 2) {
            std::cerr << "Error: Malformed line in bond list (need at least 2 tokens: id1 id2): " << line << std::endl;
            exit(EXIT_FAILURE);
        }

        // Parse bond information
        int id1 = Nfunction::String_to_Int(tokens[0]);
        int id2 = Nfunction::String_to_Int(tokens[1]);
        // Support formats: "id1 id2", "id1 id2 k", or "id1 id2 k l0"
        double k = (tokens.size() >= 3) ? Nfunction::String_to_Double(tokens[2]) : default_k;
        double l0 = (tokens.size() >= 4) ? Nfunction::String_to_Double(tokens[3]) : default_l0;

        // Helper function to find vertex by index (handles both active vertex indices and original IDs)
        auto findVertex = [&](int id) -> vertex* {
            // First, try as active vertex index (0-based)
            if (id >= 0 && id < static_cast<int>(numVertices)) {
                return activeVertices[id];
            }
            // If out of range, try to map from original vertex ID to active vertex index
            const std::map<int, int>& originalToActiveMap = m_pState->GetMesh()->GetOriginalToActiveVertexMap();
            auto it = originalToActiveMap.find(id);
            if (it != originalToActiveMap.end()) {
                int active_index = it->second;
                if (active_index >= 0 && active_index < static_cast<int>(numVertices)) {
                    return activeVertices[active_index];
                } else {
                    std::cerr << "  Warning: Original ID " << id << " maps to active index " << active_index 
                              << " which is out of range (0-" << (numVertices-1) << ")" << std::endl;
                }
            } else {
                std::cerr << "  Warning: Original ID " << id << " not found in mapping. Map size: " 
                          << originalToActiveMap.size() << std::endl;
            }
            return nullptr;
        };

        vertex* v1 = findVertex(id1);
        vertex* v2 = findVertex(id2);

        if (v1 == nullptr || v2 == nullptr) {
            std::cerr << "Error: Vertex index out of range in bond list: " << line << std::endl;
            std::cerr << "  Active vertices: 0 to " << (numVertices - 1) << std::endl;
            std::cerr << "  Bond file uses indices: " << id1 << " " << id2 << std::endl;
            const std::map<int, int>& originalToActiveMap = m_pState->GetMesh()->GetOriginalToActiveVertexMap();
            std::cerr << "  Mapping contains " << originalToActiveMap.size() << " entries" << std::endl;
            if (!originalToActiveMap.empty()) {
                std::cerr << "  First few mappings: ";
                int count = 0;
                for (auto it = originalToActiveMap.begin(); it != originalToActiveMap.end() && count < 5; ++it, ++count) {
                    std::cerr << "[" << it->first << "->" << it->second << "] ";
                }
                std::cerr << std::endl;
                // Check if the requested IDs are in the map
                if (originalToActiveMap.find(id1) != originalToActiveMap.end()) {
                    std::cerr << "  ID " << id1 << " IS in mapping, maps to " << originalToActiveMap.at(id1) << std::endl;
                } else {
                    std::cerr << "  ID " << id1 << " is NOT in mapping" << std::endl;
                }
                if (originalToActiveMap.find(id2) != originalToActiveMap.end()) {
                    std::cerr << "  ID " << id2 << " IS in mapping, maps to " << originalToActiveMap.at(id2) << std::endl;
                } else {
                    std::cerr << "  ID " << id2 << " is NOT in mapping" << std::endl;
                }
            }
            std::cerr << "  Note: Bond file should use active vertex indices (0-based) or original vertex IDs." << std::endl;
            std::cerr << "  If using original IDs, they will be automatically mapped to active vertex indices." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Create and store the bond
        bond b(bondId++, v1, v2, k, l0);
        m_AllBonds.push_back(b);
    }

    // Associate bonds with their vertices
    for (bond& b : m_AllBonds) {
        b.GetV1()->AddBondToList(&b);
        b.GetV2()->AddBondToList(&b);
    }

}
std::string HarmonicBondsList::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+m_FileName;
    return state;
}
double HarmonicBondsList::GetTotalEnergy(){
    
    double en = 0;
    
    for (std::vector<bond>::iterator it = m_AllBonds.begin() ; it != m_AllBonds.end(); ++it){
        en += it->CalculateEnergy();
    }
    
    return en;
}
