#include "HarmonicBondsList.h"
#include "State.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <set>

// E = -k*(Field*Inc_direction)^2 --> we could make it k1*(Field*Inc_direction)-k2*(Field*Inc_direction)^2

HarmonicBondsList::HarmonicBondsList(State *pState, std::string filename) {
    m_FileName = filename;
    m_pState = pState;
}
HarmonicBondsList::~HarmonicBondsList()
{
    // Clean up angle pointers
    for (angle* a : m_AllAngles) {
        delete a;
    }
    m_AllAngles.clear();
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
    double default_k_angle = 0.0;  // Angle stiffness (0 = disabled)
    double default_theta0 = M_PI;  // Equilibrium angle (180 degrees = straight chain)
    bool default_k_angle_set = false;
    bool default_theta0_set = false;
    std::string angleFileName = "";  // Optional separate angle file

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
            
            // Look for k_angle = value
            size_t k_angle_pos = lower_line.find("k_angle");
            if (k_angle_pos != std::string::npos) {
                size_t eq_pos = lower_line.find("=", k_angle_pos);
                if (eq_pos != std::string::npos) {
                    std::string val_str = line.substr(eq_pos + 1);
                    size_t comment_pos = val_str.find(";");
                    if (comment_pos != std::string::npos) {
                        val_str = val_str.substr(0, comment_pos);
                    }
                    std::vector<std::string> val_tokens = Nfunction::Split(val_str);
                    if (!val_tokens.empty()) {
                        default_k_angle = Nfunction::String_to_Double(val_tokens[0]);
                        default_k_angle_set = true;
                    }
                }
            }
            
            // Look for theta0 = value
            size_t theta0_pos = lower_line.find("theta0");
            if (theta0_pos != std::string::npos) {
                size_t eq_pos = lower_line.find("=", theta0_pos);
                if (eq_pos != std::string::npos) {
                    std::string val_str = line.substr(eq_pos + 1);
                    size_t comment_pos = val_str.find(";");
                    if (comment_pos != std::string::npos) {
                        val_str = val_str.substr(0, comment_pos);
                    }
                    std::vector<std::string> val_tokens = Nfunction::Split(val_str);
                    if (!val_tokens.empty()) {
                        double theta0_val = Nfunction::String_to_Double(val_tokens[0]);
                        // Assume degrees if > 2*PI
                        if (theta0_val > 2.0 * M_PI) {
                            default_theta0 = theta0_val * M_PI / 180.0;
                        } else {
                            default_theta0 = theta0_val;
                        }
                        default_theta0_set = true;
                    }
                }
            }
            
            // Look for angle_file = value
            size_t angle_file_pos = lower_line.find("angle_file");
            if (angle_file_pos != std::string::npos) {
                size_t eq_pos = lower_line.find("=", angle_file_pos);
                if (eq_pos != std::string::npos) {
                    std::string val_str = line.substr(eq_pos + 1);
                    size_t comment_pos = val_str.find(";");
                    if (comment_pos != std::string::npos) {
                        val_str = val_str.substr(0, comment_pos);
                    }
                    std::vector<std::string> val_tokens = Nfunction::Split(val_str);
                    if (!val_tokens.empty()) {
                        angleFileName = val_tokens[0];
                    }
                }
            }
            
            continue;
        }

        // Check if this is an SMC bond (marked with ;SMC comment)
        bool is_smc_bond = false;
        std::string line_to_parse = line;
        size_t comment_pos = line_to_parse.find(";");
        if (comment_pos != std::string::npos) {
            std::string comment_part = line_to_parse.substr(comment_pos);
            if (comment_part.find("SMC") != std::string::npos) {
                is_smc_bond = true;
            }
            line_to_parse = line_to_parse.substr(0, comment_pos);
        }
        
        // Split the line into components
        std::vector<std::string> tokens = Nfunction::Split(line_to_parse);
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

        // Track SMC bonds (store in canonical order: smaller ID first)
        if (is_smc_bond) {
            if (id1 < id2) {
                m_SMCBonds.insert(std::make_pair(id1, id2));
            } else {
                m_SMCBonds.insert(std::make_pair(id2, id1));
            }
        }

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
    
    // Display bond parameters
    std::cout << "---> Using bond parameters: k_stretch=" << default_k 
              << ", l0=" << default_l0 << std::endl;
    std::cout << "---> Loaded " << m_AllBonds.size() << " bonds" << std::endl;
    
    // Load angles from separate file if provided, otherwise generate from bonds
    if (!angleFileName.empty()) {
        std::cout << "---> Loading angles from file: " << angleFileName << std::endl;
        LoadAnglesFromFile(angleFileName);
        std::cout << "---> Loaded " << m_AllAngles.size() << " angles from file" << std::endl;
    } else if (default_k_angle_set && default_k_angle > 0.0) {
        std::cout << "---> Generating angles from bonds (excluding SMC bonds)" << std::endl;
        std::cout << "---> Using angle parameters: k_angle=" << default_k_angle 
                  << ", theta0=" << default_theta0 << " rad (" << (default_theta0 * 180.0 / M_PI) << " deg)" << std::endl;
        CreateAnglesFromBonds(default_k_angle, default_theta0);
        std::cout << "---> Created " << m_AllAngles.size() << " angles from bonds" << std::endl;
    } else {
        std::cout << "---> No angles loaded (angle_file not specified and k_angle not set or zero)" << std::endl;
    }

}
std::string HarmonicBondsList::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+m_FileName;
    return state;
}
double HarmonicBondsList::GetTotalEnergy(){
    
    double en = 0;
    double bond_en = 0;
    double angle_en = 0;
    
    // Sum bond energies
    for (std::vector<bond>::iterator it = m_AllBonds.begin() ; it != m_AllBonds.end(); ++it){
        bond_en += it->CalculateEnergy();
    }
    
    // Sum angle energies
    for (std::vector<angle*>::iterator it = m_AllAngles.begin() ; it != m_AllAngles.end(); ++it){
        angle_en += (*it)->CalculateEnergy();
    }
    
    en = bond_en + angle_en;
    
    return en;
}

void HarmonicBondsList::LoadAnglesFromFile(const std::string& angleFileName) {
    std::ifstream input(angleFileName);
    if (!input) {
        std::cerr << "Error: Failed to open angle file: " << angleFileName << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Get active vertices
    const std::vector<vertex*>& activeVertices = m_pState->GetMesh()->GetActiveV();
    const size_t numVertices = activeVertices.size();
    
    // Default angle parameters
    double default_k_angle = 20.0;
    double default_theta0 = 0.0;
    
    std::string line;
    int angleId = 0;
    
    // Read the file line by line
    while (std::getline(input, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == ';') {
            // Try to parse default values from comments
            std::string lower_line = line;
            std::transform(lower_line.begin(), lower_line.end(), lower_line.begin(), ::tolower);
            
            // Look for k_angle = value
            size_t k_angle_pos = lower_line.find("k_angle");
            if (k_angle_pos != std::string::npos) {
                size_t eq_pos = lower_line.find("=", k_angle_pos);
                if (eq_pos != std::string::npos) {
                    std::string val_str = line.substr(eq_pos + 1);
                    size_t comment_pos = val_str.find(";");
                    if (comment_pos != std::string::npos) {
                        val_str = val_str.substr(0, comment_pos);
                    }
                    std::vector<std::string> val_tokens = Nfunction::Split(val_str);
                    if (!val_tokens.empty()) {
                        default_k_angle = Nfunction::String_to_Double(val_tokens[0]);
                    }
                }
            }
            
            // Look for theta0 = value
            size_t theta0_pos = lower_line.find("theta0");
            if (theta0_pos != std::string::npos) {
                size_t eq_pos = lower_line.find("=", theta0_pos);
                if (eq_pos != std::string::npos) {
                    std::string val_str = line.substr(eq_pos + 1);
                    size_t comment_pos = val_str.find(";");
                    if (comment_pos != std::string::npos) {
                        val_str = val_str.substr(0, comment_pos);
                    }
                    std::vector<std::string> val_tokens = Nfunction::Split(val_str);
                    if (!val_tokens.empty()) {
                        double theta0_val = Nfunction::String_to_Double(val_tokens[0]);
                        // Assume degrees if > 2*PI
                        if (theta0_val > 2.0 * M_PI) {
                            default_theta0 = theta0_val * M_PI / 180.0;
                        } else {
                            default_theta0 = theta0_val;
                        }
                    }
                }
            }
            continue;
        }
        
        // Strip inline comments
        std::string line_to_parse = line;
        size_t comment_pos = line_to_parse.find(";");
        if (comment_pos != std::string::npos) {
            line_to_parse = line_to_parse.substr(0, comment_pos);
        }
        
        // Split the line into components
        std::vector<std::string> tokens = Nfunction::Split(line_to_parse);
        if (tokens.size() < 3) {
            std::cerr << "Error: Malformed line in angle list (need at least 3 tokens: v1 v2 v3): " << line << std::endl;
            exit(EXIT_FAILURE);
        }
        
        // Parse angle information: v1 v2 v3 [k_angle] [theta0]
        int id1 = Nfunction::String_to_Int(tokens[0]);
        int id2 = Nfunction::String_to_Int(tokens[1]);
        int id3 = Nfunction::String_to_Int(tokens[2]);
        double k_angle = (tokens.size() >= 4) ? Nfunction::String_to_Double(tokens[3]) : default_k_angle;
        double theta0 = (tokens.size() >= 5) ? Nfunction::String_to_Double(tokens[4]) : default_theta0;
        
        // Convert theta0 from degrees if needed
        if (tokens.size() >= 5 && theta0 > 2.0 * M_PI) {
            theta0 = theta0 * M_PI / 180.0;
        }
        
        // Helper function to find vertex by index
        auto findVertex = [&](int id) -> vertex* {
            if (id >= 0 && id < static_cast<int>(numVertices)) {
                return activeVertices[id];
            }
            const std::map<int, int>& originalToActiveMap = m_pState->GetMesh()->GetOriginalToActiveVertexMap();
            auto it = originalToActiveMap.find(id);
            if (it != originalToActiveMap.end()) {
                int active_index = it->second;
                if (active_index >= 0 && active_index < static_cast<int>(numVertices)) {
                    return activeVertices[active_index];
                }
            }
            return nullptr;
        };
        
        vertex* v1 = findVertex(id1);
        vertex* v2 = findVertex(id2);
        vertex* v3 = findVertex(id3);
        
        if (v1 == nullptr || v2 == nullptr || v3 == nullptr) {
            std::cerr << "Error: Vertex index out of range in angle list: " << line << std::endl;
            exit(EXIT_FAILURE);
        }
        
        // Create and store the angle
        angle* a = new angle(angleId++, v1, v2, v3, k_angle, theta0);
        m_AllAngles.push_back(a);
        
        // Associate angle with all three vertices
        v1->AddAngleToList(a);
        v2->AddAngleToList(a);
        v3->AddAngleToList(a);
    }
}

void HarmonicBondsList::CreateAnglesFromBonds(double k_angle, double theta0) {
    // Build adjacency map from bonds (excluding SMC bonds)
    const std::vector<vertex*>& activeVertices = m_pState->GetMesh()->GetActiveV();
    std::map<int, std::set<int> > adjacency;  // Map from vertex index to neighbor indices
    
    // Build a map from vertex pointer to index for fast lookup
    std::map<vertex*, int> vertex_to_index;
    for (size_t i = 0; i < activeVertices.size(); i++) {
        vertex_to_index[activeVertices[i]] = i;
    }
    
    for (const bond& b : m_AllBonds) {
        // Skip SMC bonds when building adjacency map for angle creation
        int vid1 = b.GetV1()->GetVID();
        int vid2 = b.GetV2()->GetVID();
        std::pair<int, int> bond_pair = (vid1 < vid2) ? std::make_pair(vid1, vid2) : std::make_pair(vid2, vid1);
        
        if (m_SMCBonds.find(bond_pair) != m_SMCBonds.end()) {
            continue;  // Skip SMC bonds - don't include them in angle creation
        }
        
        // Find indices of vertices in activeVertices
        auto it1 = vertex_to_index.find(b.GetV1());
        auto it2 = vertex_to_index.find(b.GetV2());
        if (it1 != vertex_to_index.end() && it2 != vertex_to_index.end()) {
            int idx1 = it1->second;
            int idx2 = it2->second;
            adjacency[idx1].insert(idx2);
            adjacency[idx2].insert(idx1);
        }
    }
    
    // Find all chains (connected components) and create angles for triplets
    std::set<int> visited;
    int angleId = 0;
    
    for (const auto& pair : adjacency) {
        int start_idx = pair.first;
        
        // Skip if already visited
        if (visited.find(start_idx) != visited.end()) {
            continue;
        }
        
        // Find the chain starting from this vertex using DFS
        std::vector<int> chain;
        std::set<int> component_visited;
        std::vector<int> stack;
        stack.push_back(start_idx);
        
        while (!stack.empty()) {
            int current = stack.back();
            stack.pop_back();
            
            if (component_visited.find(current) != component_visited.end()) {
                continue;
            }
            
            component_visited.insert(current);
            chain.push_back(current);
            visited.insert(current);
            
            // Add neighbors to stack
            if (adjacency.find(current) != adjacency.end()) {
                for (int neighbor : adjacency[current]) {
                    if (component_visited.find(neighbor) == component_visited.end()) {
                        stack.push_back(neighbor);
                    }
                }
            }
        }
        
        // For each chain, create angles for consecutive triplets
        // Only create angles if chain has at least 3 vertices
        if (chain.size() >= 3) {
            // Check if chain is circular (first and last are connected)
            bool is_circular = (adjacency[chain[0]].find(chain.back()) != adjacency[chain[0]].end());
            
            if (is_circular && chain.size() >= 3) {
                // Circular chain: create angles for all triplets including wraparound
                for (size_t i = 0; i < chain.size(); i++) {
                    size_t prev_idx = (i == 0) ? chain.size() - 1 : i - 1;
                    size_t curr_idx = i;
                    size_t next_idx = (i == chain.size() - 1) ? 0 : i + 1;
                    
                    if (chain[prev_idx] < (int)activeVertices.size() && 
                        chain[curr_idx] < (int)activeVertices.size() && 
                        chain[next_idx] < (int)activeVertices.size()) {
                        vertex* v1 = activeVertices[chain[prev_idx]];
                        vertex* v2 = activeVertices[chain[curr_idx]];
                        vertex* v3 = activeVertices[chain[next_idx]];
                        
                        if (v1 && v2 && v3) {
                            angle* a = new angle(angleId++, v1, v2, v3, k_angle, theta0);
                            m_AllAngles.push_back(a);
                            
                            // Associate angle with ALL three vertices
                            v1->AddAngleToList(a);
                            v2->AddAngleToList(a);
                            v3->AddAngleToList(a);
                        }
                    }
                }
            } else {
                // Linear chain: create angles for triplets (skip first and last)
                for (size_t i = 1; i < chain.size() - 1; i++) {
                    vertex* v1 = activeVertices[chain[i-1]];
                    vertex* v2 = activeVertices[chain[i]];
                    vertex* v3 = activeVertices[chain[i+1]];
                    
                    if (v1 && v2 && v3) {
                        angle* a = new angle(angleId++, v1, v2, v3, k_angle, theta0);
                        m_AllAngles.push_back(a);
                        
                        // Associate angle with all three vertices
                        v1->AddAngleToList(a);
                        v2->AddAngleToList(a);
                        v3->AddAngleToList(a);
                    }
                }
            }
        }
    }
}
