#include "HarmonicBondsList.h"
#include "State.h"
#include <map>
#include <set>
#include <cmath>

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
    // Get active vertices
    const std::vector<vertex*>& activeVertices = m_pState->GetMesh()->GetActiveV();
    const size_t numVertices = activeVertices.size();

    // Create a map from vertex VID to vertex pointer
    // The bond file uses vertex VIDs, not array indices
    std::map<int, vertex*> vid_to_vertex;
    for (size_t i = 0; i < numVertices; i++) {
        vid_to_vertex[activeVertices[i]->GetVID()] = activeVertices[i];
    }

    std::string line;
    int bondId = 0;
    
    // Global bond parameters (defaults, can be overridden per bond)
    double global_k_stretch = 10.0;
    double global_l0 = 0.0;
    double global_k_angle = 0.0;  // Angle stiffness (0 = disabled)
    double global_theta0 = M_PI;  // Equilibrium angle (180 degrees = straight chain)
    bool global_k_stretch_set = false;
    bool global_l0_set = false;
    bool global_k_angle_set = false;
    bool global_theta0_set = false;

    // First pass: read global parameters from comments (both ; and #)
    std::ifstream input_first(m_FileName);
    if (!input_first) {
        std::cerr << "Error: Failed to open the input file: " << m_FileName << std::endl;
        exit(EXIT_FAILURE);
    }
    
    while (std::getline(input_first, line)) {
        // Look for global parameter declarations in comments (; or #)
        std::string comment;
        if (line.find(";") != std::string::npos) {
            comment = line.substr(line.find(";") + 1);
        } else if (line.find("#") != std::string::npos) {
            comment = line.substr(line.find("#") + 1);
        } else {
            continue;
        }
        
        // Look for "k_stretch=" or "k_stretch ="
        size_t k_stretch_pos = comment.find("k_stretch=");
        if (k_stretch_pos == std::string::npos) k_stretch_pos = comment.find("k_stretch =");
        if (k_stretch_pos != std::string::npos) {
            std::string k_stretch_str = comment.substr(k_stretch_pos + (comment[k_stretch_pos+9] == '=' ? 10 : 11));
            std::vector<std::string> k_stretch_tokens = Nfunction::Split(k_stretch_str);
            if (!k_stretch_tokens.empty()) {
                global_k_stretch = Nfunction::String_to_Double(k_stretch_tokens[0]);
                global_k_stretch_set = true;
            }
        }
        // Look for "l0=" or "l0 ="
        size_t l0_pos = comment.find("l0=");
        if (l0_pos == std::string::npos) l0_pos = comment.find("l0 =");
        if (l0_pos != std::string::npos) {
            std::string l0_str = comment.substr(l0_pos + (comment[l0_pos+2] == '=' ? 3 : 4));
            std::vector<std::string> l0_tokens = Nfunction::Split(l0_str);
            if (!l0_tokens.empty()) {
                global_l0 = Nfunction::String_to_Double(l0_tokens[0]);
                global_l0_set = true;
            }
        }
        // Look for "k_angle=" or "k_angle ="
        size_t k_angle_pos = comment.find("k_angle=");
        if (k_angle_pos == std::string::npos) k_angle_pos = comment.find("k_angle =");
        if (k_angle_pos != std::string::npos) {
            std::string k_angle_str = comment.substr(k_angle_pos + (comment[k_angle_pos+7] == '=' ? 8 : 9));
            std::vector<std::string> k_angle_tokens = Nfunction::Split(k_angle_str);
            if (!k_angle_tokens.empty()) {
                global_k_angle = Nfunction::String_to_Double(k_angle_tokens[0]);
                global_k_angle_set = true;
            }
        }
        // Look for "theta0=" or "theta0 ="
        size_t theta0_pos = comment.find("theta0=");
        if (theta0_pos == std::string::npos) theta0_pos = comment.find("theta0 =");
        if (theta0_pos != std::string::npos) {
            std::string theta0_str = comment.substr(theta0_pos + (comment[theta0_pos+7] == '=' ? 8 : 9));
            std::vector<std::string> theta0_tokens = Nfunction::Split(theta0_str);
            if (!theta0_tokens.empty()) {
                // theta0 can be in degrees or radians - assume degrees if > 2*PI
                double theta0_val = Nfunction::String_to_Double(theta0_tokens[0]);
                if (theta0_val > 2.0 * M_PI) {
                    global_theta0 = theta0_val * M_PI / 180.0;  // Convert degrees to radians
                } else {
                    global_theta0 = theta0_val;
                }
                global_theta0_set = true;
            }
        }
    }
    input_first.close();
    
    // Debug: print global parameters
    if (global_k_stretch_set || global_l0_set) {
        std::cout << "---> Using global bond parameters: k_stretch=" << global_k_stretch << ", l0=" << global_l0 << std::endl;
    }
    if (global_k_angle_set || global_theta0_set) {
        std::cout << "---> Using global angle parameters: k_angle=" << global_k_angle 
                  << ", theta0=" << global_theta0 << " rad (" << (global_theta0 * 180.0 / M_PI) << " deg)" << std::endl;
    }
    
    // Open file for reading bonds
    std::ifstream input(m_FileName);
    if (!input) {
        std::cerr << "Error: Failed to open the input file: " << m_FileName << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read the file line by line
    while (std::getline(input, line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        
        // Check if this is an SMC bond (marked with ;SMC or #SMC comment)
        bool is_smc_bond = false;
        std::string clean_line = line;
        if (line.find(";SMC") != std::string::npos || line.find("#SMC") != std::string::npos) {
            is_smc_bond = true;
            // Remove the comment part for parsing
            size_t comment_pos = line.find(";");
            if (comment_pos == std::string::npos) comment_pos = line.find("#");
            if (comment_pos != std::string::npos) {
                clean_line = line.substr(0, comment_pos);
            }
        }
        
        // Skip pure comment lines (lines that start with ; or # and don't have bond data)
        if (clean_line.empty() || (line[0] == ';' && !is_smc_bond) || (line[0] == '#' && !is_smc_bond)) {
            continue;
        }

        // Split the line into components
        std::vector<std::string> tokens = Nfunction::Split(clean_line);
        if (tokens.size() < 2) {
            std::cerr << "Error: Malformed line in bond list (need at least 2 tokens): " << line << " (tokens: " << tokens.size() << ")" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Parse bond information
        // The bond file uses vertex VIDs, not array indices
        int vid1 = Nfunction::String_to_Int(tokens[0]);
        int vid2 = Nfunction::String_to_Int(tokens[1]);
        
        // Track SMC bonds (store in canonical order: smaller VID first)
        if (is_smc_bond) {
            if (vid1 < vid2) {
                m_SMCBonds.insert(std::make_pair(vid1, vid2));
            } else {
                m_SMCBonds.insert(std::make_pair(vid2, vid1));
            }
            std::cout << "---> Marked SMC bond: " << vid1 << " - " << vid2 << std::endl;
        }
        
        // Use global parameters if not specified per bond
        // If global parameters weren't set and not provided per bond, use defaults
        double k = (tokens.size() > 2) ? Nfunction::String_to_Double(tokens[2]) : 
                   (global_k_stretch_set ? global_k_stretch : 10.0);
        double l0 = (tokens.size() > 3) ? Nfunction::String_to_Double(tokens[3]) : 
                    (global_l0_set ? global_l0 : 0.0);

        // Find vertices by VID (not array index)
        if (vid_to_vertex.find(vid1) == vid_to_vertex.end()) {
            std::cerr << "Error: Vertex VID " << vid1 << " not found in bond list: " << line << std::endl;
            exit(EXIT_FAILURE);
        }
        if (vid_to_vertex.find(vid2) == vid_to_vertex.end()) {
            std::cerr << "Error: Vertex VID " << vid2 << " not found in bond list: " << line << std::endl;
            exit(EXIT_FAILURE);
        }
        
        vertex* v1 = vid_to_vertex[vid1];
        vertex* v2 = vid_to_vertex[vid2];

        // Create and store the bond
        bond b(bondId++, v1, v2, k, l0);
        m_AllBonds.push_back(b);
    }

    // Associate bonds with their vertices
    for (bond& b : m_AllBonds) {
        b.GetV1()->AddBondToList(&b);
        b.GetV2()->AddBondToList(&b);
    }
    
    // Create angles for DNA chains if angle potential is enabled
    if (global_k_angle_set && global_k_angle > 0.0) {
        std::cout << "---> Creating angles for DNA chains with k_angle=" << global_k_angle 
                  << ", theta0=" << global_theta0 << " rad (" << (global_theta0 * 180.0 / M_PI) << " deg)" << std::endl;
        CreateAnglesForChains(global_k_angle, global_theta0);
        std::cout << "---> Created " << m_AllAngles.size() << " angles for DNA chains" << std::endl;
    } else {
        std::cout << "---> Angle potential disabled (k_angle=" << global_k_angle 
                  << ", set=" << global_k_angle_set << ")" << std::endl;
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
    
    // Debug output (only print occasionally to avoid spam)
    static int call_count = 0;
    if (++call_count % 1000 == 0 && m_AllAngles.size() > 0) {
        std::cout << "---> [HarmonicBondsList] Total energy: bonds=" << bond_en 
                  << ", angles=" << angle_en << ", total=" << en << std::endl;
    }
    
    return en;
}

void HarmonicBondsList::CreateAnglesForChains(double k_angle, double theta0) {
    // Build adjacency map from bonds (using vertex indices in activeVertices, not VID)
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
    
    // Find all chains (connected components)
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
                            
                            // Associate angle with ALL three vertices (v1, v2, v3)
                            // This ensures energy is recalculated when any vertex in the angle moves
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
                        
                        // Associate angle with middle vertex
                        v2->AddAngleToList(a);
                    }
                }
            }
        }
    }
}
