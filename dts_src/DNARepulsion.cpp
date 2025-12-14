#include <unordered_set>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "DNARepulsion.h"
#include "State.h"
#include "bond.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Constructor
DNARepulsion::DNARepulsion(State* pstate, std::string input_data) 
    : m_pState(pstate), m_Input_Data(input_data) {
}

DNARepulsion::~DNARepulsion() {
    // Destructor
}

void DNARepulsion::Initialize() {
    std::vector<std::string> data = Nfunction::Split(m_Input_Data);
    
    // Prevent out-of-bounds access
    if (data.size() < 3) {
        std::cerr << "---> error: insufficient input data for "
                  << this->GetDefaultReadName() << ". Expected: EP R0_membrane R0_dna [power]\n";
        std::cerr << "  EP = repulsion strength\n";
        std::cerr << "  R0_membrane = cutoff radius for membrane-DNA interactions (nm)\n";
        std::cerr << "  R0_dna = cutoff radius for DNA-DNA interactions (nm)\n";
        std::cerr << "  Optional: power = exponent for power law (default 6, higher = steeper repulsion)\n";
        exit(-1);
    }

    // Parse parameters: EP (repulsion strength), R0_membrane (cutoff for membrane-DNA), R0_dna (cutoff for DNA-DNA)
    m_EP = Nfunction::String_to_Double(data[0]);
    m_R0 = Nfunction::String_to_Double(data[1]);
    m_R0_2 = m_R0 * m_R0;
    m_R0_DNA = Nfunction::String_to_Double(data[2]);
    m_R0_DNA_2 = m_R0_DNA * m_R0_DNA;
    
    // Optional parameter: power law exponent
    m_Power = (data.size() > 3) ? Nfunction::String_to_Double(data[3]) : 6.0;  // Default: 1/r^6 (LJ-like)

    m_pBox = m_pState->GetMesh()->GetBox();

    // Improved error messages and checks
    if (m_EP < 0) {
        std::cout << "---> error: assigned repulsion strength in "
                  << this->GetDefaultReadName() << " should be greater than or equal to zero.\n";
        exit(-1);
    }
    if (m_R0 < 0) {
        std::cout << "---> error: assigned membrane-DNA cutoff radius in "
                  << this->GetDefaultReadName() << " should be greater than zero.\n";
        exit(-1);
    }
    if (m_R0_DNA < 0) {
        std::cout << "---> error: assigned DNA-DNA cutoff radius in "
                  << this->GetDefaultReadName() << " should be greater than zero.\n";
        exit(-1);
    }

    // Separate vertices into membrane vertices and DNA beads
    const std::vector<vertex*>& allVertices = m_pState->GetMesh()->GetActiveV();
    m_pMembraneVertices.clear();
    m_pDNABeads.clear();
    m_MembraneVertexSet.clear();
    m_DNABeadSet.clear();
    m_NeighborExclusionCache.clear();
    
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
            m_MembraneVertexSet.insert(v);  // Add to hash set for O(1) lookup
        }
        // Pure DNA bead: has bonds but no links (free polymer)
        // Note: We exclude vertices that are both (attached DNA) from repulsion
        if (hasBonds) {
            // DNA bead - check if it has links (attached DNA) or not (free DNA)
            // But we can't call GetVLinkList() safely, so assume pure DNA bead if it has bonds
            // In practice, DNA beads in our setup don't have links
            m_pDNABeads.push_back(v);
            m_DNABeadSet.insert(v);  // Add to hash set for O(1) lookup
        }
    }
    
    std::cout << "---> Initialized " << this->GetDefaultReadName() 
              << " with " << m_pMembraneVertices.size() << " membrane vertices and "
              << m_pDNABeads.size() << " DNA beads.\n";
    std::cout << "--->   Membrane-DNA cutoff: " << m_R0 << " nm, DNA-DNA cutoff: " << m_R0_DNA << " nm\n";
    std::cout << "--->   Repulsion form: E(r) = EP * (R0/r)^power for r < R0, E(r) = 0 for r >= R0\n";
    std::cout << "--->   Power exponent: " << m_Power << " (higher = steeper repulsion at short distances)\n";
    std::cout << "--->   Performance optimizations: Hash sets enabled, neighbor caching enabled\n";

    return;
}

bool DNARepulsion::IsMembraneVertex(vertex* v) const {
    // Performance optimization: Use hash set for O(1) lookup instead of O(N) linear search
    return m_MembraneVertexSet.find(v) != m_MembraneVertexSet.end();
}

bool DNARepulsion::IsDNABead(vertex* v) const {
    // Performance optimization: Use hash set for O(1) lookup instead of O(N) linear search
    return m_DNABeadSet.find(v) != m_DNABeadSet.end();
}

bool DNARepulsion::ValidateCoordinates(vertex* v) const {
    if (v == nullptr) return false;
    
    double x = v->GetXPos();
    double y = v->GetYPos();
    double z = v->GetZPos();
    
    // Check for NaN, infinity, or extremely large values (signs of corruption)
    return !(std::isnan(x) || std::isnan(y) || std::isnan(z) ||
             std::isinf(x) || std::isinf(y) || std::isinf(z) ||
             std::abs(x) > 1e10 || std::abs(y) > 1e10 || std::abs(z) > 1e10);
}

const std::unordered_set<vertex*>& DNARepulsion::GetNeighborExclusionSet(vertex* pV) const {
    // Safety check: ensure vertex is valid
    if (pV == nullptr) {
        static std::unordered_set<vertex*> empty_set;
        return empty_set;
    }
    
    // Performance optimization: Cache neighbor exclusion sets to avoid recomputing
    // Thread-safety: Protect cache lookup with critical section
#ifdef _OPENMP
    const std::unordered_set<vertex*>* found_cached_set = nullptr;
    #pragma omp critical(dna_repulsion_cache)
    {
        auto it = m_NeighborExclusionCache.find(pV);
        if (it != m_NeighborExclusionCache.end()) {
            found_cached_set = &(it->second);  // Found in cache
        }
    }
    if (found_cached_set != nullptr) {
        return *found_cached_set;  // Return cached set
    }
#else
    auto it = m_NeighborExclusionCache.find(pV);
    if (it != m_NeighborExclusionCache.end()) {
        return it->second;  // Return cached set
    }
#endif
    
    // Compute and cache the exclusion set
    std::unordered_set<vertex*> neighboringBeads;
    
    // Safety check: only compute for DNA beads
    if (!IsDNABead(pV)) {
        // Not a DNA bead, return empty set and cache it
        m_NeighborExclusionCache[pV] = neighboringBeads;
        return m_NeighborExclusionCache[pV];
    }
    
    // Get directly bonded neighbors
    // Safety: Check if vertex is valid and has bonds before accessing
    if (pV == nullptr) {
        static std::unordered_set<vertex*> empty_set;
        return empty_set;
    }
    
    try {
        // Safety: Check if GetBonds() is safe to call
        std::vector<bond*> bonds;
        try {
            bonds = pV->GetBonds();
        } catch (...) {
            // If GetBonds() throws, return empty set
            m_NeighborExclusionCache[pV] = neighboringBeads;
            return m_NeighborExclusionCache[pV];
        }
        
        // Safety: Check if bonds vector is valid
        if (bonds.empty()) {
            // No bonds, return empty set and cache it
            m_NeighborExclusionCache[pV] = neighboringBeads;
            return m_NeighborExclusionCache[pV];
        }
        
        for (auto* b : bonds) {
            if (b == nullptr) continue;
            
            // Safety: Check bond vertices
            vertex* v1 = nullptr;
            vertex* v2 = nullptr;
            try {
                v1 = b->GetV1();
                v2 = b->GetV2();
            } catch (...) {
                continue;  // Skip if we can't get bond vertices
            }
            
            if (v1 == nullptr || v2 == nullptr) continue;
            
            vertex* other = (v1 == pV) ? v2 : v1;
            if (other == nullptr || other == pV) continue;
            
            neighboringBeads.insert(other);
            
            // Also get next-nearest neighbors (beads bonded to the direct neighbors)
            try {
                std::vector<bond*> neighborBonds;
                try {
                    neighborBonds = other->GetBonds();
                } catch (...) {
                    continue;  // Skip if we can't get neighbor bonds
                }
                
                for (auto* nb : neighborBonds) {
                    if (nb == nullptr) continue;
                    
                    // Safety: Check neighbor bond vertices
                    vertex* nb_v1 = nullptr;
                    vertex* nb_v2 = nullptr;
                    try {
                        nb_v1 = nb->GetV1();
                        nb_v2 = nb->GetV2();
                    } catch (...) {
                        continue;  // Skip if we can't get neighbor bond vertices
                    }
                    
                    if (nb_v1 == nullptr || nb_v2 == nullptr) continue;
                    
                    vertex* nextNeighbor = (nb_v1 == other) ? nb_v2 : nb_v1;
                    if (nextNeighbor != nullptr && nextNeighbor != pV) {
                        neighboringBeads.insert(nextNeighbor);
                    }
                }
            } catch (...) {
                // Skip if we can't process neighbor bonds
                continue;
            }
        }
    } catch (...) {
        // If anything fails, return empty set
        m_NeighborExclusionCache[pV] = neighboringBeads;
        return m_NeighborExclusionCache[pV];
    }
    
    // Cache the result (mutable allows modification of cache in const function)
    // Thread-safety: Use OpenMP critical section to protect both cache reads and writes
    // This ensures that even if multiple threads try to cache the same vertex simultaneously,
    // only one will write, and subsequent reads will get the cached value
#ifdef _OPENMP
    const std::unordered_set<vertex*>* result_set = nullptr;
    #pragma omp critical(dna_repulsion_cache)
    {
        // Double-check after acquiring lock (another thread might have cached it)
        auto it_check = m_NeighborExclusionCache.find(pV);
        if (it_check != m_NeighborExclusionCache.end()) {
            result_set = &(it_check->second);
        } else {
            m_NeighborExclusionCache[pV] = neighboringBeads;
            result_set = &(m_NeighborExclusionCache[pV]);
        }
    }
    // Return the cached set (pointer is valid because map doesn't rehash during our critical section)
    return *result_set;
#else
    m_NeighborExclusionCache[pV] = neighboringBeads;
    return m_NeighborExclusionCache[pV];
#endif
}

std::vector<vertex*> DNARepulsion::FindNearbyDNABeads(vertex* pV) const {
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

                    // Add DNA beads that are NOT direct neighbors and within cutoff
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
                            // Check distance to ensure it's within cutoff (voxels might be larger than cutoff)
                            double dist2 = pV->SquareDistanceFromAVertex(v);
                            if (std::isnan(dist2) || std::isinf(dist2)) {
                                continue;
                            }
                            if (dist2 <= m_R0_2) {  // Only add if within cutoff
                                nearbyBeads.push_back(v);
                            }
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

double DNARepulsion::CalculateRepulsionEnergy(vertex* v1, vertex* v2, double cutoff_radius_squared) const {
    // Performance optimization: Validate coordinates using helper function (coordinates already validated in FindNearby)
    // But keep check here for safety
    if (v1 == nullptr || v2 == nullptr) {
        return 0;
    }
    
    // Performance optimization: Coordinates should already be validated in FindNearby functions
    // Skip redundant validation here for speed
    
    double Dist2 = v1->SquareDistanceFromAVertex(v2);
    
    // Check if distance calculation produced invalid result
    if (std::isnan(Dist2) || std::isinf(Dist2) || Dist2 == 0) {
        return 0;
    }
    
    // Performance optimization: Early cutoff check before expensive sqrt
    // R0 is the cutoff radius - interactions beyond R0 are not calculated
    if (cutoff_radius_squared < Dist2) return 0;  // Beyond cutoff radius, no interaction
    
    double Dist = sqrt(Dist2);
    
    // Power law repulsion potential: E(r) = EP * (R0/r)^power for r < R0
    // This creates a steep repulsive potential that prevents strand crossings
    // Higher power values create steeper repulsion at short distances
    // At r = R0, E = EP (the repulsion strength parameter)
    // As r -> 0, E -> infinity (very strong repulsion prevents overlap)
    double R0 = sqrt(cutoff_radius_squared);
    double E = m_EP * std::pow(R0 / Dist, m_Power);
    
    // Check if energy calculation produced invalid result
    if (std::isnan(E) || std::isinf(E)) {
        std::cerr << "---> WARNING: Repulsion energy calculation produced NaN/Inf. "
                  << "Dist=" << Dist << ", EP=" << m_EP << ", Power=" << m_Power << std::endl;
        return 0;
    }
    
    return E;
}

std::vector<vertex*> DNARepulsion::FindNearbyMembraneVertices(vertex* pV) const {
    std::vector<vertex*> nearbyMembrane;
    
    // Performance optimization: Validate coordinates once at the start
    if (!ValidateCoordinates(pV)) {
        return nearbyMembrane;  // Early return if vertex coordinates are invalid
    }
    
    // Since DNA bead voxels are now updated when they move, we can trust voxel-based search
    // Only use fallback if there's no voxel (shouldn't happen, but safety check)
    
    if (pV->GetVoxel() != nullptr) {
        // Use fast voxel-based search - voxel should be up-to-date since we update it
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

                    // Add membrane vertices
                    for (vertex* v : CV) {
                        if (v == nullptr || v == pV) {
                            continue;  // Skip null vertices and self
                        }
                        
                        // Performance optimization: Use hash set for O(1) membrane vertex check
                        if (!IsMembraneVertex(v)) {
                            continue;  // Skip non-membrane vertices
                        }
                        
                        // Check distance to ensure it's within cutoff
                        double dist2 = pV->SquareDistanceFromAVertex(v);
                        if (std::isnan(dist2) || std::isinf(dist2)) {
                            continue;  // Skip if distance is invalid
                        }
                        
                        if (dist2 <= m_R0_2) {  // Use membrane-DNA cutoff
                            nearbyMembrane.push_back(v);
                        }
                    }
                }
            }
        }
    } else {
        // Fallback: Only if no voxel (shouldn't happen, but safety check)
        for (auto* v_membrane : m_MembraneVertexSet) {
            if (v_membrane == nullptr || v_membrane == pV) {
                continue;
            }
            
            // Performance optimization: Validate coordinates once
            if (!ValidateCoordinates(v_membrane)) {
                continue;  // Skip corrupted vertices
            }
            
            double dist2 = pV->SquareDistanceFromAVertex(v_membrane);
            
            // Check if distance calculation produced invalid result
            if (std::isnan(dist2) || std::isinf(dist2)) {
                continue;  // Skip if distance is invalid
            }
            
            if (dist2 <= m_R0_2) {  // Use membrane-DNA cutoff
                nearbyMembrane.push_back(v_membrane);
            }
        }
    }

    return nearbyMembrane;
}

std::vector<vertex*> DNARepulsion::FindNearbyDNABeadsForDNA(vertex* pV) const {
    std::vector<vertex*> nearbyBeads;
    
    // Performance optimization: Validate coordinates once at the start
    if (!ValidateCoordinates(pV)) {
        return nearbyBeads;  // Early return if vertex coordinates are invalid
    }
    
    // Performance optimization: Use cached neighbor exclusion set
    const std::unordered_set<vertex*>& neighboringBeads = GetNeighborExclusionSet(pV);
    
    // DNA bead voxels are updated when they move (for efficient nonbonded interaction searches)
    // Only use fallback if there's no voxel (shouldn't happen, but safety check)
    
    if (pV->GetVoxel() != nullptr) {
        // Use fast voxel-based search - voxel should be up-to-date since we update it
        // Loop through all 26 neighboring voxel cells around the vertex's current voxel
        for (int n = -1; n < 2; n++) {
            for (int m = -1; m < 2; m++) {
                for (int s = -1; s < 2; s++) {
                    Voxel<vertex>* neighbor_voxel = pV->GetVoxel()->GetANeighbourCell(n, m, s);
                    if (neighbor_voxel == nullptr) {
                        continue;
                    }
                    std::vector<vertex*> CV = neighbor_voxel->GetContentObjects();

                    for (vertex* v : CV) {
                        if (v == nullptr || v == pV) {
                            continue;
                        }
                        
                        // Performance optimization: Skip excluded neighbors early
                        if (neighboringBeads.find(v) != neighboringBeads.end()) {
                            continue;
                        }
                        
                        // Performance optimization: Use hash set for O(1) DNA bead check
                        if (!IsDNABead(v)) {
                            continue;  // Skip non-DNA beads
                        }
                        
                        // Check distance to ensure it's within DNA-DNA cutoff
                        double dist2 = pV->SquareDistanceFromAVertex(v);
                        if (std::isnan(dist2) || std::isinf(dist2)) {
                            continue;
                        }
                        
                        if (dist2 <= m_R0_DNA_2) {  // Use DNA-DNA cutoff
                            nearbyBeads.push_back(v);
                        }
                    }
                }
            }
        }
    } else {
        // Fallback: Only if no voxel (shouldn't happen, but safety check)
        for (auto* v_dna : m_DNABeadSet) {
            if (v_dna == nullptr || v_dna == pV) {
                continue;
            }
            
            // Performance optimization: Skip excluded neighbors early
            if (neighboringBeads.find(v_dna) != neighboringBeads.end()) {
                continue;
            }
            
            // Performance optimization: Validate coordinates once
            if (!ValidateCoordinates(v_dna)) {
                continue;
            }
            
            double dist2 = pV->SquareDistanceFromAVertex(v_dna);
            
            if (std::isnan(dist2) || std::isinf(dist2)) {
                continue;
            }
            
            if (dist2 <= m_R0_DNA_2) {  // Use DNA-DNA cutoff
                nearbyBeads.push_back(v_dna);
            }
        }
    }

    return nearbyBeads;
}

double DNARepulsion::GetVertexNonBondedEnergy(vertex* pvertex) const {
    // Safety check: ensure vertex is valid
    if (pvertex == nullptr) {
        return 0;
    }
    
    double dE = 0;
    
    int vid = -1;
    try {
        vid = pvertex->GetVID();
    } catch (...) {
        return 0;
    }
    
    // Repulsion: membrane-DNA and DNA-DNA
    if (IsMembraneVertex(pvertex)) {
        // Membrane vertex: find nearby DNA beads and calculate repulsion
        try {
            std::vector<vertex*> nearbyBeads = FindNearbyDNABeads(pvertex);
            for (auto* v_dna : nearbyBeads) {
                if (v_dna != nullptr) {
                    dE += CalculateRepulsionEnergy(pvertex, v_dna, m_R0_2);  // Use membrane-DNA cutoff
                }
            }
        } catch (...) {
            // If FindNearbyDNABeads fails, continue without repulsion
        }
    } else if (IsDNABead(pvertex)) {
        // DNA bead: find nearby membrane vertices and calculate repulsion
        try {
            std::vector<vertex*> nearbyMembrane = FindNearbyMembraneVertices(pvertex);
            for (auto* v_membrane : nearbyMembrane) {
                if (v_membrane != nullptr) {
                    // Use the same CalculateRepulsionEnergy function (it's symmetric)
                    double repulsion = CalculateRepulsionEnergy(v_membrane, pvertex, m_R0_2);  // Use membrane-DNA cutoff
                    dE += repulsion;
                }
            }
        } catch (...) {
            // If FindNearbyMembraneVertices fails, continue
        }
        
        // NEW: Also calculate repulsion with other DNA beads (excluding bonded neighbors)
        try {
            std::vector<vertex*> nearbyDNABeads = FindNearbyDNABeadsForDNA(pvertex);
            for (auto* v_dna : nearbyDNABeads) {
                if (v_dna != nullptr) {
                    // Use the same repulsion formula (symmetric) with DNA-DNA cutoff
                    double repulsion = CalculateRepulsionEnergy(pvertex, v_dna, m_R0_DNA_2);  // Use DNA-DNA cutoff
                    dE += repulsion;
                }
            }
        } catch (...) {
            // If FindNearbyDNABeadsForDNA fails, continue without DNA-DNA repulsion
        }
    }
    // If vertex is neither (e.g., attached DNA), return 0
    
    return dE;
}

// Implementation of CurrentState function
std::string DNARepulsion::CurrentState() const {
    std::string state = GetBaseDefaultReadName() + " = " + this->GetDerivedDefaultReadName();
    state += " " + m_Input_Data;
    return state;
}

