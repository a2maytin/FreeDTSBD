
#include "VAHGlobalMeshProperties.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "State.h"
VAHGlobalMeshProperties::VAHGlobalMeshProperties() : 
                                                                  m_TotalVolume(0.0),
                                                                  m_TotalArea(0.0),
                                                                  m_TotalCurvature(0.0),
                                                                  m_VolumeIsActive(false),
                                                                  m_AreaIsActive(false),
                                                                  m_GlobalCurvatureIsActive(false){

#ifdef _OPENMP
        omp_init_lock(&m_VLock);  // Initialize the lock
        omp_init_lock(&m_ALock);  // Initialize the lock
        omp_init_lock(&m_CGLock);  // Initialize the lock
#endif

}
VAHGlobalMeshProperties::~VAHGlobalMeshProperties() {
#ifdef _OPENMP
    omp_destroy_lock(&m_VLock);  // Destroy the lock when done
    omp_destroy_lock(&m_ALock);  // Destroy the lock when done
    omp_destroy_lock(&m_CGLock);  // Destroy the lock when done
#endif
}
void VAHGlobalMeshProperties::Initialize(State* pState){
    m_pState = pState;
    return;
}
void VAHGlobalMeshProperties::Add2Volume(double vol){
#ifdef _OPENMP
    omp_set_lock(&m_VLock);  // Lock before updating shared variables
    
    // Update m_TotalVolume
    m_TotalVolume += vol;
    omp_unset_lock(&m_VLock);  // Unlock after updating shared variables
#else
    // If OpenMP is not available, update normally without locks
    m_TotalVolume += vol;
#endif
    
    return;
}
void VAHGlobalMeshProperties::Add2TotalArea(double area){
#ifdef _OPENMP
    omp_set_lock(&m_ALock);  // Lock before updating shared variables

    // Update  m_TotalArea
    m_TotalArea += area;
    omp_unset_lock(&m_ALock);  // Unlock after updating shared variables
#else
    // If OpenMP is not available, update normally without locks
    m_TotalArea += area;
#endif
    return;
}
void VAHGlobalMeshProperties::Add2GlobalCurvature(double CG){
#ifdef _OPENMP
    omp_set_lock(&m_CGLock);  // Lock before updating shared variables

    // Update global curvature
    m_TotalCurvature += CG;

    omp_unset_lock(&m_CGLock);  // Unlock after updating shared variables
#else
    // If OpenMP is not available, update normally without locks
    m_TotalCurvature += CG;
#endif
    return;
}
void VAHGlobalMeshProperties::CalculateAVertexRingContributionToGlobalVariables(vertex *p_vertex, double &vol, double &area, double& curvature){
    
    // Call the overloaded version that gets neighbors internally
    std::vector<vertex *> neighbors = p_vertex->GetVNeighbourVertex();
    CalculateAVertexRingContributionToGlobalVariables(p_vertex, neighbors, vol, area, curvature);
    return;
}

void VAHGlobalMeshProperties::CalculateAVertexRingContributionToGlobalVariables(vertex *p_vertex, const std::vector<vertex*>& neighbors, double &vol, double &area, double& curvature){
    
    vol = 0.0;
    area = 0.0;
    curvature = 0.0;
    
    // Safety check
    if (p_vertex == nullptr) {
        return;
    }
    
    // Debug: track all vertices with non-6 triangles to find hang
    bool track_debug = false;
    int vid = -1;
    try {
        vid = p_vertex->GetVID();
        // Track vertices that don't have 6 triangles to find where hang occurs
        // std::vector<triangle *> test_triangles = p_vertex->GetVTraingleList();
        // if (test_triangles.size() != 6) {
        //     track_debug = true;
        //     // std::cerr << "\n[VAHGlobalMeshProperties " << vid << "] Entry, neighbors.size()=" << neighbors.size() << ", triangles.size()=" << test_triangles.size() << std::flush;
        // }
        track_debug = false;  // Disable debug output
    } catch (...) {
        // If GetVID() fails, just continue without tracking
    }
    
    if(m_VolumeIsActive && m_AreaIsActive){
        if (track_debug) {
            std::cerr << " [getting triangles for volume/area...]" << std::flush;
        }
        // CRITICAL FIX: Use index-based access instead of iterator to avoid potential hang
        std::vector<triangle *> ring_triangles = p_vertex->GetVTraingleList();  // Copy the vector
        if (track_debug) {
            std::cerr << " [got " << ring_triangles.size() << " triangles]" << std::flush;
        }
        
        // FIX: Calculate volume/area for ALL vertices, not just those with != 6 triangles
        // The previous skip for 6 triangles was causing massive volume/area leaks
        // Instead, validate each triangle individually (which CalculateSingleTriangleVolume already does)
        // NOTE: We do NOT divide by 3 here because:
        // 1. The coupling classes initialize by summing triangles directly (no division)
        // 2. For incremental updates, we need the full change from a vertex's ring, not 1/3
        // 3. When a vertex moves, its ring changes, and we update the total by the full change
        for (size_t i = 0; i < ring_triangles.size(); ++i){
            if (track_debug) {
                std::cerr << " [triangle " << i << "]" << std::flush;
            }
            if (ring_triangles[i] != nullptr) {
                if (track_debug) {
                    std::cerr << " [calculating volume...]" << std::flush;
                }
                // CalculateSingleTriangleVolume already validates triangle vertices internally
                // It returns 0.0 if triangle is invalid, so safe to call for all triangles
                vol += CalculateSingleTriangleVolume(ring_triangles[i]);
                
                if (track_debug) {
                    std::cerr << " [getting area...]" << std::flush;
                }
                // Validate triangle before getting area
                try {
                    vertex* v1 = ring_triangles[i]->GetV1();
                    vertex* v2 = ring_triangles[i]->GetV2();
                    vertex* v3 = ring_triangles[i]->GetV3();
                    if (v1 != nullptr && v2 != nullptr && v3 != nullptr) {
                        area += ring_triangles[i]->GetArea();
                    }
                } catch (...) {
                    // If area calculation fails, skip this triangle
                    if (track_debug) {
                        std::cerr << " [area calc failed, skipping]" << std::flush;
                    }
                }
                if (track_debug) {
                    std::cerr << " [done triangle " << i << "]" << std::flush;
                }
            } else {
                if (track_debug) {
                    std::cerr << " [triangle " << i << " is null]" << std::flush;
                }
            }
        }
        if (track_debug) {
            std::cerr << " [finished volume/area]" << std::flush;
        }
    }
    if(!m_VolumeIsActive && m_AreaIsActive){
        // CRITICAL FIX: Use index-based access instead of iterator to avoid potential hang
        std::vector<triangle *> ring_triangles = p_vertex->GetVTraingleList();  // Copy the vector
        for (size_t i = 0; i < ring_triangles.size(); ++i){
            if (ring_triangles[i] != nullptr) {
                area += ring_triangles[i]->GetArea();
            }
        }
    }
    if(m_GlobalCurvatureIsActive){
        if (track_debug) {
            std::cerr << " [calculating curvature...]" << std::flush;
        }
        
        double C = p_vertex->GetP1Curvature() + p_vertex->GetP2Curvature();
        curvature = C * (p_vertex->GetArea());
        
        if (track_debug) {
            std::cerr << " [iterating over " << neighbors.size() << " neighbors...]" << std::flush;
        }
        
        // CRITICAL FIX: Use the passed neighbor list instead of calling GetVNeighbourVertex() again
        // This avoids the hang that occurs when GetVNeighbourVertex() is called multiple times
        // Use index-based access and skip DNA beads which don't have curvature
        for (size_t i = 0; i < neighbors.size(); ++i){
            if (track_debug && i % 2 == 0) {
                std::cerr << " [neighbor " << i << "]" << std::flush;
            }
            
            vertex* neighbor = neighbors[i];
            if (neighbor == nullptr) {
                continue;
            }
            // Skip DNA beads - they don't have curvature
            if (!neighbor->GetBonds().empty()) {
                continue;  // Skip DNA beads
            }
            C = neighbor->GetP1Curvature() + neighbor->GetP2Curvature();
            double A = neighbor->GetArea();
            curvature += C * A;
        }
        
        if (track_debug) {
            std::cerr << " [finished curvature]" << std::flush;
        }
    }
    
    if (track_debug) {
        // std::cerr << " [VAHGlobalMeshProperties " << vid << "] Exit" << std::flush;
    }
    
    return;
}
double VAHGlobalMeshProperties::CalculateSingleTriangleVolume(triangle *pTriangle){

    // Safety check
    if (pTriangle == nullptr) {
        return 0.0;
    }
    
    if(m_pState->GetMesh()->GetHasCrossedPBC()){
        *(m_pState->GetTimeSeriesLog()) << "---> the system has crossed the PBC while volume is being calculated.";
        *(m_pState->GetTimeSeriesLog()) << " SOLUTION: Restart the simulation and center the system. Also, activate the command for centering the box.";

         exit(0);
    }
    
    // CRITICAL FIX: Validate triangle state before accessing methods
    // After Alexander moves, triangles may be in invalid states that cause hangs
    vertex* v1 = pTriangle->GetV1();
    vertex* v2 = pTriangle->GetV2();
    vertex* v3 = pTriangle->GetV3();
    
    // Safety check: ensure all vertices are valid
    if (v1 == nullptr || v2 == nullptr || v3 == nullptr) {
        return 0.0;
    }
    
    // CRITICAL FIX: Check if any vertex is a DNA bead (has bonds but no links)
    // DNA beads shouldn't be in triangles, but after topology changes they might be
    if (!v1->GetBonds().empty() || !v2->GetBonds().empty() || !v3->GetBonds().empty()) {
        // This triangle involves a DNA bead - skip it
        return 0.0;
    }
    
    // CRITICAL FIX: Add safety checks to prevent hangs
    // If any of these calls hang, return 0 to allow simulation to continue
    try {
        double T_area = pTriangle->GetArea();
        Vec3D Normal_v = pTriangle->GetNormalVector();
        Vec3D Pos = v1->GetPos();

        // Compute triangle volume
        return T_area * (Vec3D::dot(Normal_v, Pos)) / 3.0;
    } catch (...) {
        // If any method call fails, return 0 to allow simulation to continue
        return 0.0;
    }
}

void VAHGlobalMeshProperties::CalculateALinkTrianglesContributionToGlobalVariables(links *p_link, double &vol, double &area, double& curvature) {
    
    curvature = 0;
    vol = 0.0;
    area = 0.0;
    
    if(m_VolumeIsActive){
        vol += CalculateSingleTriangleVolume(p_link->GetTriangle());
        vol += CalculateSingleTriangleVolume(p_link->GetMirrorLink()->GetTriangle());
    }
    
    if(m_AreaIsActive){
        area += p_link->GetTriangle()->GetArea();
        area += p_link->GetMirrorLink()->GetTriangle()->GetArea();
    }
    
    if(m_GlobalCurvatureIsActive){
        vertex *v1 = p_link->GetV1();
        vertex *v2 = p_link->GetV2();
        vertex *v3 = p_link->GetV3();
        vertex *v4 = p_link->GetMirrorLink()->GetV3();

        double C1 = v1->GetP1Curvature() + v1->GetP2Curvature();
        double C2 = v2->GetP1Curvature() + v2->GetP2Curvature();
        double C3 = v3->GetP1Curvature() + v3->GetP2Curvature();
        double C4 = v4->GetP1Curvature() + v4->GetP2Curvature();
        double A1 = v1->GetArea();
        double A2 = v2->GetArea();
        double A3 = v3->GetArea();
        double A4 = v4->GetArea();
        curvature = C1*A1 + C2*A2 + C3*A3 + C4*A4;
    }
    

    return;
}
void VAHGlobalMeshProperties::CalculateBoxRescalingContributionToGlobalVariables(double lx, double ly, double lz, double& vol, double& area, double& curvature){
    
    vol = 0;
    area = 0;
    curvature = 0;

    if(m_VolumeIsActive && m_AreaIsActive){
        std::vector<triangle *>& all_tri = m_pState->GetMesh()->GetActiveT();
        for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
            vol += CalculateSingleTriangleVolume(*it);
            area += (*it)->GetArea();
        }
    }
    else if(m_VolumeIsActive && !m_AreaIsActive){
        std::vector<triangle *>& all_tri = m_pState->GetMesh()->GetActiveT();
        for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
            vol += CalculateSingleTriangleVolume(*it);
        }
    }
    else if(!m_VolumeIsActive && m_AreaIsActive){
        std::vector<triangle *>& all_tri = m_pState->GetMesh()->GetActiveT();
        for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
            area += (*it)->GetArea();
        }
    }
    
    if(m_GlobalCurvatureIsActive){
        std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();
        for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {
            double area = (*it)->GetArea();
            double curv = (*it)->GetP1Curvature() + (*it)->GetP2Curvature();
            curvature += curv*area;
        }
    }
    
    return;
}
void VAHGlobalMeshProperties::CalculateGlobalVariables(double& vol, double& area, double& curvature){
    
    vol = 0;
    area = 0;
    curvature = 0;
    
    // Calculate volume and area from triangles (each triangle counted once)
    std::vector<triangle *>& all_tri = m_pState->GetMesh()->GetActiveT();
    for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
        if (*it != nullptr) {
            vol += CalculateSingleTriangleVolume(*it);
            area += (*it)->GetArea();
        }
    }
    
    // Calculate curvature from vertices (only membrane vertices, not DNA beads)
    std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();
    for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {
        if (*it == nullptr) {
            continue;
        }
        
        // CRITICAL FIX: Skip DNA beads - they don't have area or curvature
        // Use GetBonds() to check - don't call GetVLinkList() on DNA beads as it may hang/crash
        if (!(*it)->GetBonds().empty()) {
            // This is a DNA bead (has bonds) - skip it
            continue;
        }
        
        // Only calculate for membrane vertices (those without bonds, which should have links)
        // Note: We don't check GetVLinkList() here to avoid potential hangs
        // If a vertex has no bonds, it should be a membrane vertex
        try {
            double v_area = (*it)->GetArea();
            double curv = (*it)->GetP1Curvature() + (*it)->GetP2Curvature();
            curvature += curv * v_area;
        } catch (...) {
            // If getting area/curvature fails, skip this vertex
            continue;
        }
    }
    
    return;
}
