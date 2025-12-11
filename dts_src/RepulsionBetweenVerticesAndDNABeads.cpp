#include <unordered_set>
#include <cmath>
#include "RepulsionBetweenVerticesAndDNABeads.h"
#include "State.h"

// Constructor
RepulsionBetweenVerticesAndDNABeads::RepulsionBetweenVerticesAndDNABeads(State* pstate, std::string input_data) 
    : m_pState(pstate), m_Input_Data(input_data) {
}

RepulsionBetweenVerticesAndDNABeads::~RepulsionBetweenVerticesAndDNABeads() {
    // Destructor
}

void RepulsionBetweenVerticesAndDNABeads::Initialize() {
    std::vector<std::string> data = Nfunction::Split(m_Input_Data);
    
    // Prevent out-of-bounds access
    if (data.size() < 2) {
        std::cerr << "---> error: insufficient input data for "
                  << this->GetDefaultReadName() << ". Expected: EP R0\n";
        exit(-1);
    }

    // Parse parameters: EP (repulsion strength) and R0 (cutoff radius)
    m_EP = Nfunction::String_to_Double(data[0]);
    m_R0 = Nfunction::String_to_Double(data[1]);
    m_R0_2 = m_R0 * m_R0;

    m_pBox = m_pState->GetMesh()->GetBox();

    // Improved error messages and checks
    if (m_EP < 0) {
        std::cout << "---> error: assigned repulsion strength in "
                  << this->GetDefaultReadName() << " should be greater than or equal to zero.\n";
        exit(-1);
    }
    if (m_R0 < 0) {
        std::cout << "---> error: assigned cutoff radius in "
                  << this->GetDefaultReadName() << " should be greater than zero.\n";
        exit(-1);
    }

    // Separate vertices into membrane vertices and DNA beads
    const std::vector<vertex*>& allVertices = m_pState->GetMesh()->GetActiveV();
    m_pMembraneVertices.clear();
    m_pDNABeads.clear();
    
    for (auto* v : allVertices) {
        // CRITICAL FIX: Use GetBonds() to identify DNA beads first, don't call GetVLinkList() on DNA beads
        bool hasBonds = !v->GetBonds().empty();
        
        // Only check for links if it's not a DNA bead (to avoid hangs)
        bool hasLinks = false;
        if (!hasBonds) {
            // Not a DNA bead, safe to check for links
            hasLinks = !v->GetVLinkList().empty();
        }
        
        // Pure membrane vertex: has links but no bonds
        if (hasLinks && !hasBonds) {
            m_pMembraneVertices.push_back(v);
        }
        // Pure DNA bead: has bonds but no links (free polymer)
        // Note: We exclude vertices that are both (attached DNA) from repulsion
        if (hasBonds) {
            // DNA bead - check if it has links (attached DNA) or not (free DNA)
            // But we can't call GetVLinkList() safely, so assume pure DNA bead if it has bonds
            // In practice, DNA beads in our setup don't have links
            m_pDNABeads.push_back(v);
        }
    }
    
    std::cout << "---> Initialized " << this->GetDefaultReadName() 
              << " with " << m_pMembraneVertices.size() << " membrane vertices and "
              << m_pDNABeads.size() << " DNA beads.\n";
    // std::cerr << "---> [REPULSION DEBUG] Initialized with " << m_pMembraneVertices.size() 
    //           << " membrane vertices and " << m_pDNABeads.size() << " DNA beads." << std::endl;
    // if (m_pDNABeads.size() > 0) {
    //     std::cerr << "---> [REPULSION DEBUG] First DNA bead ID=" << m_pDNABeads[0]->GetVID() 
    //               << " at (" << m_pDNABeads[0]->GetXPos() << ", " 
    //               << m_pDNABeads[0]->GetYPos() << ", " << m_pDNABeads[0]->GetZPos() << ")" << std::endl;
    // }

    return;
}

bool RepulsionBetweenVerticesAndDNABeads::IsMembraneVertex(vertex* v) const {
    // CRITICAL FIX: Use the stored membrane vertex list instead of checking bonds/links at runtime
    // Check if this vertex is in our stored membrane vertex list
    for (auto* mem_vertex : m_pMembraneVertices) {
        if (mem_vertex == v) {
            return true;
        }
    }
    return false;
    
    // OLD CODE (unreliable - GetBonds() sometimes returns empty):
    // if (!v->GetBonds().empty()) {
    //     return false;  // Has bonds, so it's a DNA bead, not a membrane vertex
    // }
    // return true;  // Assume it's a membrane vertex if it has no bonds
}

bool RepulsionBetweenVerticesAndDNABeads::IsDNABead(vertex* v) const {
    // CRITICAL FIX: Use the stored DNA bead list instead of checking bonds at runtime
    // Bonds might not be accessible or might be cleared, but we know which vertices are DNA beads from initialization
    // Check if this vertex is in our stored DNA bead list
    for (auto* dna_bead : m_pDNABeads) {
        if (dna_bead == v) {
            return true;
        }
    }
    return false;
    
    // OLD CODE (unreliable - GetBonds() sometimes returns empty):
    // std::vector<bond*> bonds = v->GetBonds();
    // return !bonds.empty();
}

std::vector<vertex*> RepulsionBetweenVerticesAndDNABeads::FindNearbyDNABeads(vertex* pV) const {
    std::vector<vertex*> nearbyBeads;
    
    // If voxelization is available, use it for efficient neighbor search
    if (pV->GetVoxel() != nullptr) {
        // Get directly connected neighbor vertices (to exclude them)
        std::vector<vertex*> npvertex = pV->GetVNeighbourVertex();
        std::unordered_set<vertex*> npSet(npvertex.begin(), npvertex.end());

        // Loop through all 26 neighboring voxel cells around the vertex's current voxel
        for (int n = -1; n < 2; n++) {
            for (int m = -1; m < 2; m++) {
                for (int s = -1; s < 2; s++) {
                    // Safety check: ensure we can get the neighbor cell
                    Voxel<vertex>* neighbor_voxel = pV->GetVoxel()->GetANeighbourCell(n, m, s);
                    if (neighbor_voxel == nullptr) {
                        continue;  // Skip if neighbor voxel is null
                    }
                    // Get all vertices inside this neighboring cell
                    std::vector<vertex*> CV = neighbor_voxel->GetContentObjects();

                    // Add DNA beads that are NOT direct neighbors
                    // Check if vertex is a pure DNA bead (has bonds but no links)
                    for (vertex* v : CV) {
                        if (v == nullptr) {
                            continue;  // Skip null vertices
                        }
                        
                        // CRITICAL FIX: Validate DNA bead coordinates before adding
                        double x = v->GetXPos();
                        double y = v->GetYPos();
                        double z = v->GetZPos();
                        
                        // Check for NaN, infinity, or extremely large values (signs of corruption)
                        if (std::isnan(x) || std::isnan(y) || std::isnan(z) ||
                            std::isinf(x) || std::isinf(y) || std::isinf(z) ||
                            std::abs(x) > 1e10 || std::abs(y) > 1e10 || std::abs(z) > 1e10) {
                            continue;  // Skip corrupted DNA beads
                        }
                        
                        // CRITICAL FIX: Use GetBonds() to identify DNA beads, don't call GetVLinkList() on them
                        // DNA beads have bonds but no links (membrane vertices have links)
                        if (npSet.find(v) == npSet.end() && 
                            !v->GetBonds().empty()) {
                            // This is a DNA bead (has bonds)
                            // We don't check GetVLinkList() as it can cause hangs on DNA beads
                            nearbyBeads.push_back(v);
                        }
                    }
                }
            }
        }
    } else {
        // Fallback: check all DNA beads (less efficient but works)
        for (auto* v_dna : m_pDNABeads) {
            if (v_dna == nullptr || v_dna == pV) {
                continue;
            }
            
            // CRITICAL FIX: Validate DNA bead coordinates before calculating distance
            double x = v_dna->GetXPos();
            double y = v_dna->GetYPos();
            double z = v_dna->GetZPos();
            
            // Check for NaN, infinity, or extremely large values (signs of corruption)
            if (std::isnan(x) || std::isnan(y) || std::isnan(z) ||
                std::isinf(x) || std::isinf(y) || std::isinf(z) ||
                std::abs(x) > 1e10 || std::abs(y) > 1e10 || std::abs(z) > 1e10) {
                std::cerr << "---> WARNING: DNA bead ID=" << v_dna->GetVID() 
                          << " has corrupted coordinates (" << x << ", " << y << ", " << z 
                          << "). Skipping in FindNearbyDNABeads." << std::endl;
                continue;  // Skip corrupted DNA beads
            }
            
            double dist2 = pV->SquareDistanceFromAVertex(v_dna);
            
            // Check if distance calculation produced invalid result
            if (std::isnan(dist2) || std::isinf(dist2)) {
                continue;  // Skip if distance is invalid
            }
            
            if (dist2 <= m_R0_2) {
                nearbyBeads.push_back(v_dna);
            }
        }
    }

    return nearbyBeads;
}

double RepulsionBetweenVerticesAndDNABeads::CalculateRepulsionEnergy(vertex* v_membrane, vertex* v_dna) const {
    // CRITICAL FIX: Validate DNA bead coordinates before calculating distance
    // Corrupted coordinates can cause hangs or NaN values
    if (v_dna == nullptr || v_membrane == nullptr) {
        return 0;
    }
    
    // Check if DNA bead coordinates are valid (not corrupted/uninitialized)
    double x = v_dna->GetXPos();
    double y = v_dna->GetYPos();
    double z = v_dna->GetZPos();
    
    // Check for NaN, infinity, or extremely large values (signs of corruption)
    if (std::isnan(x) || std::isnan(y) || std::isnan(z) ||
        std::isinf(x) || std::isinf(y) || std::isinf(z) ||
        std::abs(x) > 1e10 || std::abs(y) > 1e10 || std::abs(z) > 1e10) {
        std::cerr << "---> WARNING: DNA bead ID=" << v_dna->GetVID() 
                  << " has corrupted coordinates (" << x << ", " << y << ", " << z 
                  << "). Skipping repulsion calculation." << std::endl;
        return 0;
    }
    
    double Dist2 = v_membrane->SquareDistanceFromAVertex(v_dna);
    
    // Check if distance calculation produced invalid result
    if (std::isnan(Dist2) || std::isinf(Dist2)) {
        std::cerr << "---> WARNING: Distance calculation produced NaN/Inf for DNA bead ID=" 
                  << v_dna->GetVID() << ". Skipping repulsion calculation." << std::endl;
        return 0;
    }
    
    if (Dist2 == 0) return 0;  // Prevent division by zero
    if (m_R0_2 < Dist2) return 0;  // Beyond cutoff radius
    
    // Simple repulsion: E = EP / Dist2
    // This creates a repulsive force that decreases with distance squared
    double E = m_EP / Dist2;
    
    // Check if energy calculation produced invalid result
    if (std::isnan(E) || std::isinf(E)) {
        std::cerr << "---> WARNING: Repulsion energy calculation produced NaN/Inf. "
                  << "Dist2=" << Dist2 << ", EP=" << m_EP << std::endl;
        return 0;
    }
    
    return E;
}

std::vector<vertex*> RepulsionBetweenVerticesAndDNABeads::FindNearbyMembraneVertices(vertex* pV) const {
    std::vector<vertex*> nearbyMembrane;
    
    // CRITICAL FIX: For DNA beads that may have escaped, always use fallback method
    // that checks ALL membrane vertices, not just those in nearby voxels
    // This ensures DNA beads far from membrane still feel repulsion
    // The voxel-based search is more efficient but can miss distant interactions
    
    // DEBUG: Check if DNA bead has voxel
    static int voxel_debug_count = 0;
    bool has_voxel = (pV->GetVoxel() != nullptr);
    if (voxel_debug_count < 5 && IsDNABead(pV)) {
        std::cerr << "---> DNA bead ID=" << pV->GetVID() 
                  << " has voxel: " << (has_voxel ? "YES" : "NO")
                  << ", total membrane vertices: " << m_pMembraneVertices.size() << std::endl;
        voxel_debug_count++;
    }
    
    // First try voxel-based search for efficiency (if DNA bead is close to membrane)
    if (pV->GetVoxel() != nullptr) {
        // CRITICAL FIX: DNA beads don't have membrane neighbors, so npSet will be empty
        // This is fine - we just want to avoid duplicates, and DNA beads won't have membrane neighbors anyway
        std::unordered_set<vertex*> npSet;
        // Note: We don't call GetVNeighbourVertex() on DNA beads as it might cause issues
        // DNA beads only have bonds, not membrane links

        // Loop through all 26 neighboring voxel cells around the vertex's current voxel
        for (int n = -1; n < 2; n++) {
            for (int m = -1; m < 2; m++) {
                for (int s = -1; s < 2; s++) {
                    // Safety check: ensure we can get the neighbor cell
                    Voxel<vertex>* neighbor_voxel = pV->GetVoxel()->GetANeighbourCell(n, m, s);
                    if (neighbor_voxel == nullptr) {
                        continue;  // Skip if neighbor voxel is null
                    }
                    // Get all vertices inside this neighboring cell
                    std::vector<vertex*> CV = neighbor_voxel->GetContentObjects();

                    // Add membrane vertices that are NOT direct neighbors
                    // Check if vertex is a pure membrane vertex (has links but no bonds)
                    for (vertex* v : CV) {
                        if (v == nullptr) {
                            continue;  // Skip null vertices
                        }
                        
                        // Validate membrane vertex coordinates before adding
                        double x = v->GetXPos();
                        double y = v->GetYPos();
                        double z = v->GetZPos();
                        
                        // Check for NaN, infinity, or extremely large values (signs of corruption)
                        if (std::isnan(x) || std::isnan(y) || std::isnan(z) ||
                            std::isinf(x) || std::isinf(y) || std::isinf(z) ||
                            std::abs(x) > 1e10 || std::abs(y) > 1e10 || std::abs(z) > 1e10) {
                            continue;  // Skip corrupted vertices
                        }
                        
                        // CRITICAL FIX: Check bonds first to avoid calling GetVLinkList() on DNA beads
                        // Membrane vertices have no bonds (DNA beads have bonds)
                        if (npSet.find(v) == npSet.end() && 
                            v->GetBonds().empty()) {
                            // This is likely a membrane vertex (no bonds)
                            // We don't check GetVLinkList() as it can cause hangs, but if it has no bonds, it's likely a membrane vertex
                            // Check distance to ensure it's within cutoff
                            double dist2 = pV->SquareDistanceFromAVertex(v);
                            if (!std::isnan(dist2) && !std::isinf(dist2) && dist2 <= m_R0_2) {
                                nearbyMembrane.push_back(v);
                            }
                        }
                    }
                }
            }
        }
    }
    
    // FIX: Always also check all membrane vertices to catch DNA beads that have escaped
    // This ensures DNA beads far from membrane still feel repulsion
    // Use a set to avoid duplicates from voxel search
    std::unordered_set<vertex*> foundSet(nearbyMembrane.begin(), nearbyMembrane.end());
    
    // DEBUG: Track how many membrane vertices are checked - ALWAYS print for DNA beads
    int checked_count = 0;
    int within_cutoff = 0;
    
    // if (IsDNABead(pV)) {
    //     std::cerr << "---> [FIND_NEARBY] Starting fallback search through " << m_pMembraneVertices.size() 
    //               << " membrane vertices for DNA bead ID=" << pV->GetVID() << std::endl;
    // }
    
    for (auto* v_membrane : m_pMembraneVertices) {
        checked_count++;
        if (v_membrane == nullptr || v_membrane == pV) {
            continue;
        }
        
        // Skip if already found via voxel search
        if (foundSet.find(v_membrane) != foundSet.end()) {
            continue;
        }
        
        // Validate membrane vertex coordinates before calculating distance
        double x = v_membrane->GetXPos();
        double y = v_membrane->GetYPos();
        double z = v_membrane->GetZPos();
        
        // Check for NaN, infinity, or extremely large values (signs of corruption)
        if (std::isnan(x) || std::isnan(y) || std::isnan(z) ||
            std::isinf(x) || std::isinf(y) || std::isinf(z) ||
            std::abs(x) > 1e10 || std::abs(y) > 1e10 || std::abs(z) > 1e10) {
            continue;  // Skip corrupted vertices
        }
        
        double dist2 = pV->SquareDistanceFromAVertex(v_membrane);
        
        // Check if distance calculation produced invalid result
        if (std::isnan(dist2) || std::isinf(dist2)) {
            continue;  // Skip if distance is invalid
        }
        
        if (dist2 <= m_R0_2) {
            nearbyMembrane.push_back(v_membrane);
            foundSet.insert(v_membrane);
            within_cutoff++;
        }
    }
    
    // DEBUG: ALWAYS print statistics for DNA beads
    // if (IsDNABead(pV)) {
    //     double x = pV->GetXPos();
    //     double y = pV->GetYPos();
    //     double z = pV->GetZPos();
    //     std::cerr << "---> [FIND_NEARBY] DNA bead ID=" << pV->GetVID() 
    //               << " at (" << x << ", " << y << ", " << z << "): "
    //               << "checked " << checked_count << " membrane vertices, "
    //               << "found " << within_cutoff << " within R0=" << sqrt(m_R0_2) << " nm" << std::endl;
    // }

    return nearbyMembrane;
}

double RepulsionBetweenVerticesAndDNABeads::GetVertexNonBondedEnergy(vertex* pvertex) const {
    // CRITICAL DEBUG: Print immediately with flush to ensure output appears
    // Use both cerr and cout to ensure we see output
    // std::cerr << "REPULSION_START" << std::endl << std::flush;
    // std::cout << "REPULSION_START_COUT" << std::endl << std::flush;
    
    double dE = 0;
    
    // DEBUG: ALWAYS print at the start to see what's happening
    int vid = -1;
    try {
        vid = pvertex->GetVID();
    } catch (...) {
        // std::cerr << "---> [REPULSION ERROR] Failed to get vertex ID!" << std::endl << std::flush;
        return 0;
    }
    
    bool is_dna = IsDNABead(pvertex);
    bool is_membrane = IsMembraneVertex(pvertex);
    int num_bonds = pvertex->GetBonds().size();
    
    // std::cerr << "---> [REPULSION] GetVertexNonBondedEnergy called for vertex ID=" << vid 
    //           << ", is_dna=" << is_dna << ", is_membrane=" << is_membrane 
    //           << ", num_bonds=" << num_bonds << std::endl << std::flush;
    
    // Two-way repulsion: both membrane vertices and DNA beads feel repulsion
    if (IsMembraneVertex(pvertex)) {
        // std::cerr << "---> [REPULSION] Processing as MEMBRANE vertex" << std::endl;
        // Membrane vertex: find nearby DNA beads and calculate repulsion
        std::vector<vertex*> nearbyBeads = FindNearbyDNABeads(pvertex);
        for (auto* v_dna : nearbyBeads) {
            dE += CalculateRepulsionEnergy(pvertex, v_dna);
        }
    } else if (IsDNABead(pvertex)) {
        // std::cerr << "---> [REPULSION] Processing as DNA BEAD" << std::endl;
        // DNA bead: find nearby membrane vertices and calculate repulsion
        // Note: The energy is the same (symmetric), but calculated from DNA bead's perspective
        std::vector<vertex*> nearbyMembrane = FindNearbyMembraneVertices(pvertex);
        
        // DEBUG: ALWAYS print for DNA beads to see what's happening - NO LIMIT
        // double x = pvertex->GetXPos();
        // double y = pvertex->GetYPos();
        // double z = pvertex->GetZPos();
        // std::cerr << "---> [REPULSION] DNA bead ID=" << pvertex->GetVID() 
        //           << " at (" << x << ", " << y << ", " << z 
        //           << ") found " << nearbyMembrane.size() 
        //           << " nearby membrane vertices (R0=" << sqrt(m_R0_2) << " nm, EP=" << m_EP << ")" << std::endl;
        
        // if (nearbyMembrane.empty()) {
        //     std::cerr << "---> [REPULSION] ERROR: DNA bead found NO nearby membrane vertices! Total membrane vertices: " 
        //               << m_pMembraneVertices.size() << std::endl;
        // }
        
        for (auto* v_membrane : nearbyMembrane) {
            // Use the same CalculateRepulsionEnergy function (it's symmetric)
            double repulsion = CalculateRepulsionEnergy(v_membrane, pvertex);
            dE += repulsion;
            
            // DEBUG: Print repulsion energy for ALL interactions
            // double dist = sqrt(v_membrane->SquareDistanceFromAVertex(pvertex));
            // std::cerr << "---> [REPULSION] DNA bead ID=" << pvertex->GetVID() 
            //           << " repulsion with membrane vertex ID=" << v_membrane->GetVID()
            //           << ": dist=" << dist << " nm, E=" << repulsion << std::endl;
        }
        
        // std::cerr << "---> [REPULSION] DNA bead ID=" << pvertex->GetVID() 
        //           << " total repulsion energy = " << dE << std::endl;
    } else {
        // std::cerr << "---> [REPULSION] Vertex is NEITHER membrane nor DNA (skipping)" << std::endl;
    }
    // If vertex is neither (e.g., attached DNA), return 0
    
    // std::cerr << "---> [REPULSION] Returning dE=" << dE << " for vertex ID=" << vid << std::endl;
    return dE;
}

// Implementation of CurrentState function
std::string RepulsionBetweenVerticesAndDNABeads::CurrentState() const {
    std::string state = GetBaseDefaultReadName() + " = " + this->GetDerivedDefaultReadName();
    state += " " + m_Input_Data;
    return state;
}

