

#include <stdio.h>
#include <thread>   // For std::this_thread::sleep_for
#include <chrono>   // For std::chrono::milliseconds
#ifdef _OPENMP
# include <omp.h>
#endif
#include "EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.h"
#include "State.h"

EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(State *pState)
    : m_pState(pState),
      m_pSurfV(pState->GetMesh()->GetSurfV()),
      m_pEdgeV(pState->GetMesh()->GetEdgeV()),
      m_Beta(pState->GetSimulation()->GetBeta()),
      m_DBeta(pState->GetSimulation()->GetDBeta()),
      m_MinLength2(pState->GetSimulation()->GetMinL2()),
      m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
      m_MinAngle(pState->GetSimulation()->GetMinAngle()){
          
          m_NumberOfMovePerStep_Surf = 1.0;
          m_NumberOfMovePerStep_Edge = 1.0;
          m_DR = 0.05;
          
      }
EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(State *pState, double rate_surf, double rate_edge, double dr)
    : m_pState(pState),
      m_pSurfV(pState->GetMesh()->GetSurfV()),
      m_pEdgeV(pState->GetMesh()->GetEdgeV()),
      m_Beta(pState->GetSimulation()->GetBeta()),
      m_DBeta(pState->GetSimulation()->GetDBeta()),
      m_MinLength2(pState->GetSimulation()->GetMinL2()),
      m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
      m_MinAngle(pState->GetSimulation()->GetMinAngle()) {
          
          m_NumberOfMovePerStep_Surf = rate_surf;
          m_NumberOfMovePerStep_Edge = rate_edge;
          m_DR = dr;
      }
EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::~EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(){

}
void EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::Initialize(){
    m_pBox = m_pState->GetMesh()->GetBox();
    m_Total_ThreadsNo = m_pState->GetThreads_Number();
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;

    return;
}
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveOneStep(int step){
    // Safety check: ensure Initialize() was called
    if (m_pBox == nullptr) {
        std::cerr << "Error: EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::Initialize() must be called before EvolveOneStep()\n";
        return false;
    }
 
#ifdef _OPENMP
    // OpenMP is available, continue with the program
    
    int no_surf_v = m_pSurfV.size();
    int no_edge_v = m_pEdgeV.size();
    int no_steps_edge = no_edge_v*m_NumberOfMovePerStep_Edge;
    int no_steps_surf = no_surf_v*m_NumberOfMovePerStep_Surf;
    
    // Debug: count DNA vs membrane vertices (removed to avoid potential race conditions)
    // static bool first_call = true;
    // if (first_call) {
    //     int dna_count = 0;
    //     int mem_count = 0;
    //     for (auto* v : m_pSurfV) {
    //         if (v->GetVLinkList().empty()) dna_count++;
    //         else mem_count++;
    //     }
    //     std::cout << "Debug: Surface vertices - DNA: " << dna_count << ", Membrane: " << mem_count << "\n";
    //     std::cout << "Debug: Will attempt " << no_steps_surf << " surface vertex moves per step\n";
    //     first_call = false;
    // }




    //omp_set_num_threads(m_Total_ThreadsNo);  // Set the number of threads
    double diff_energy = 0.0;
    int total_AcceptedMoves = 0;

    // Use OpenMP to parallelize the outer loop with dynamic scheduling and reduction
    #pragma omp parallel reduction(+:total_AcceptedMoves, diff_energy)
    {
        int local_AcceptedMoves = 0;  // Thread-local counters
        double local_diff_energy = 0.0;
        vertex* pvertex;

        // Parallel for with dynamic scheduling
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < no_steps_surf; i++) {

            // Implement a lock with backoff to avoid busy-waiting
            while (true) {
                int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_surf_v);  // Thread-local RNG
                // Safety check: ensure index is within bounds
                if (r_vid < 0 || r_vid >= no_surf_v) {
                    continue;  // Invalid index, try again
                }
                pvertex = m_pSurfV[r_vid];
                
                // Safety check: ensure vertex pointer is valid
                if (pvertex == nullptr) {
                    continue;  // Invalid vertex, try again
                }

                // Try to lock the vertex and its neighbors
                if (pvertex->CheckLockVertex()) {
                    if (!pvertex->CheckLockNeighbourVertex()) {
                        pvertex->UnlockVertex();  // Unlock if neighbors can't be locked
                    } else {
                        break;  // Successfully locked vertex and neighbors, break the loop
                    }
                }
            }
            /*int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_surf_v); // Thread-local RNG
            pvertex = m_pSurfV[r_vid];
        while (true) {
            if (pvertex->CheckLockVertex()) {
                if (!pvertex->CheckLockNeighbourVertex()) {
                    pvertex->UnlockVertex();  // Unlock the vertex if neighbors can't be locked
                } else {
                    break;  // Successfully locked vertex and neighbors, break the loop
                }
            }
        std::this_thread::sleep_for(std::chrono::nanoseconds(Backoff_Factor));
        }*/

            // Skip if the vertex is in the frozen group
            if (m_FreezGroupName == pvertex->GetGroupName()) {
                pvertex->UnlockVertex();
                pvertex->UnlockNeighbourVertex();
                continue;
            }

            // Generate random move vectors
            double dx = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            double dy = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            double dz = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

            // Check if the move is within the boundary
            if (!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(m_DR*dx, m_DR*dy, m_DR*dz, pvertex)) {
                pvertex->UnlockVertex();
                pvertex->UnlockNeighbourVertex();
                continue;
            }

            // Generate random thermal value
            double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            double tem_en = 0.0;

            // CRITICAL FIX: Use GetBonds() instead of GetVLinkList() to identify DNA beads
            // DNA beads have bonds, membrane vertices don't (in our setup)
            bool is_dna = !pvertex->GetBonds().empty();
            // Debug tracking removed for clean output
            // static int dna_attempts = 0;
            // static int dna_accepted = 0;

            // Evolve the vertex and update local counters if successful
            if (EvolveOneVertex(step, pvertex, m_DR * dx, m_DR * dy, m_DR * dz, thermal, tem_en)) {
                local_AcceptedMoves++;  // Increment local accepted moves
                local_diff_energy += tem_en;  // Accumulate local energy change
                // Debug output removed
            }

            // Unlock the vertex and its neighbors
            pvertex->UnlockVertex();
            pvertex->UnlockNeighbourVertex();
        }

        // Update the shared variables after the loop
        total_AcceptedMoves += local_AcceptedMoves;  // Reduction for accepted moves
        diff_energy += local_diff_energy;  // Reduction for energy difference
    }

    // Update the global variables once, after the parallel section
    m_AcceptedMoves += total_AcceptedMoves;  // Safely update accepted moves
    m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);  // Add energy
    m_NumberOfAttemptedMoves += no_steps_surf;  // Update the number of attempted moves
    
    // Debug: Print DNA bead statistics at end of step (removed for clean output)
    // static int step_count = 0;
    // step_count++;
    // if (step_count % 100 == 0) {
    //     static int total_dna_attempts = 0;
    //     static int total_dna_accepted = 0;
    //     #pragma omp critical
    //     {
    //         // These are thread-local, so we'd need to accumulate them differently
    //         // For now, just print a message
    //     }
    // }

omp_set_num_threads(m_Total_ThreadsNo); // Set the number of threads
 diff_energy = 0.0;

// Use OpenMP to parallelize the outer loop with dynamic scheduling
#pragma omp parallel for schedule(dynamic) num_threads(m_Total_ThreadsNo)
for (int i = 0; i < no_steps_edge; i++) {
    // Each thread has its own local counters
    int local_AcceptedMoves = 0;
    double local_diff_energy = 0.0;

    vertex *pvertex;

    // Try to lock a random vertex and its neighbors
    /*int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edge_v); // Thread-local RNG
    	pvertex = m_pEdgeV[r_vid];
    while (true) {
        if (pvertex->CheckLockVertex()) {
            if (!pvertex->CheckLockNeighbourVertex()) {
                pvertex->UnlockVertex();  // Unlock the vertex if neighbors can't be locked
            } else {
                break;  // Successfully locked vertex and neighbors, break the loop
            }
        }
    std::this_thread::sleep_for(std::chrono::nanoseconds(Backoff_Factor));
    }*/
    
    
    while (true) {
        int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edge_v);
        // Safety check: ensure index is within bounds
        if (r_vid < 0 || r_vid >= no_edge_v) {
            continue;  // Invalid index, try again
        }
        pvertex = m_pEdgeV[r_vid];
        
        // Safety check: ensure vertex pointer is valid
        if (pvertex == nullptr) {
            continue;  // Invalid vertex, try again
        }

        if (pvertex->CheckLockVertex()) {
            if (!pvertex->CheckLockNeighbourVertex()) {
                pvertex->UnlockVertex();  // Unlock the vertex if neighbors can't be locked
            } else {
                break;  // Successfully locked vertex and neighbors, break the loop
            }
        }
    }

    // Skip if vertex is in the frozen group
    if (m_FreezGroupName == pvertex->GetGroupName()) {
        pvertex->UnlockVertex();
        pvertex->UnlockNeighbourVertex();
        continue;
    }

    // Generate random move vectors
    double dx = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    double dy = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    double dz = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

    // Check if the move is within the boundary
    if (!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(m_DR*dx,m_DR*dy,m_DR*dz, pvertex)) {
        pvertex->UnlockVertex();
        pvertex->UnlockNeighbourVertex();
        continue;
    }

    double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    double tem_en = 0.0;

    // Evolve the vertex and update local counters if successful
    if (EvolveOneVertex(step, pvertex, m_DR * dx, m_DR * dy, m_DR * dz, thermal, tem_en)) {
        local_AcceptedMoves++;
        local_diff_energy += tem_en;
    }

    // Unlock the vertex and its neighbors
    pvertex->UnlockVertex();
    pvertex->UnlockNeighbourVertex();

    // Combine local results into shared variables
    #pragma omp atomic
    m_AcceptedMoves += local_AcceptedMoves;

    // Atomic update to shared diff_energy variable
    #pragma omp atomic
    diff_energy += local_diff_energy;
}

// Add total energy from the local diff_energy calculated by each thread
m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
m_NumberOfAttemptedMoves += no_steps_surf;



    return true;
#else
    std::cerr << "---> Error: OpenMP is not available, but the program requires it. Please recompile with OpenMP support.\n";
    exit(EXIT_FAILURE); // Exit with a non-zero value to indicate failure
    
    return false;
#endif
}
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp, double &changed_en){
    
    double old_energy = 0;
    double new_energy = 0;

//---> first checking if all the distances will be fine if we move the vertex
    // CRITICAL FIX: Use GetBonds() instead of GetVLinkList() to identify DNA beads
    // DNA beads have bonds, membrane vertices don't (in our setup)
    bool is_dna_vertex = !pvertex->GetBonds().empty();
    
    if(!VertexMoveIsFine(pvertex,dx,dy,dz,m_MinLength2,m_MaxLength2)) {  // this function could get a booling varaible to say, it crossed the voxel
        // Debug output removed - uncomment if needed
        // #pragma omp critical
        // {
        //     static int reject_count = 0;
        //     reject_count++;
        // }
        return 0;
    }

    //--- obtain vertices energy terms and make copies
    // Initialize variables that will be used later
    std::vector<triangle *> N_triangles;
    std::vector<links *> v_NLinks;
    std::vector<vertex *> vNeighbourV;
    
    // For DNA beads (no links), skip membrane-specific energy calculations
    if (!is_dna_vertex) {
    old_energy = pvertex->GetEnergy();
    old_energy += pvertex->GetBindingEnergy();
    pvertex->ConstantMesh_Copy();
    pvertex->Copy_VFsBindingEnergy();  // vector field
        vNeighbourV = pvertex->GetVNeighbourVertex();  // Copy the vector
    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
            // Skip DNA beads - use GetBonds() instead of GetVLinkList() to avoid hang
            if (!(*it)->GetBonds().empty()) {
                continue;  // Skip DNA beads
            }
        (*it)->ConstantMesh_Copy();
        old_energy += (*it)->GetEnergy();
        old_energy += (*it)->GetBindingEnergy();
        (*it)->Copy_VFsBindingEnergy();
    }
        N_triangles = pvertex->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
        (*it)->ConstantMesh_Copy();
    }
        v_NLinks = pvertex->GetVLinkList();  // Copy the vector
    for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
        (*it)->ConstantMesh_Copy();
        (*it)->GetNeighborLink1()->ConstantMesh_Copy();
    }
    //-- we need this to make sure all the links connected to this v is updated
    if(pvertex->GetVertexType() == 1){
        pvertex->GetPrecedingEdgeLink()->ConstantMesh_Copy();
    }
    } else {
        // DNA beads: only need bond energy, no membrane energy
        old_energy = 0;  // DNA beads don't have membrane energy
        // CRITICAL FIX: DNA beads MUST save their position before the move
        // Otherwise ReverseConstantMesh_Copy() won't work correctly if the move is rejected
        pvertex->ConstantMesh_Copy();  // Save old position for DNA beads
    }
    // find the links in which there interaction energy changes (only for membrane vertices)
    std::vector<links*> Affected_links;
    if (!is_dna_vertex) {
        Affected_links = GetEdgesWithInteractionChange(pvertex);
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy();
    }
    }
    //--> bond and nonbonded energy before move
    double bond_energy = -(pvertex->GetBondedEnergyOfVertex());  // Includes bonds + angles
    // Two-way repulsion: both membrane vertices and DNA beads feel repulsion
    double dE_nonbonded = 0.0;
    if (m_pState->GetNonbondedInteractionBetweenVertices() != nullptr) {
        // static int debug_count_old = 0;
        // bool is_dna = !pvertex->GetBonds().empty();
        // if (debug_count_old < 10 && is_dna) {
        //     std::cerr << "---> [METROPOLIS DEBUG] Calculating OLD nonbonded energy for DNA bead ID=" 
        //               << pvertex->GetVID() << std::endl;
        //     debug_count_old++;
        // }
        dE_nonbonded = -(m_pState->GetNonbondedInteractionBetweenVertices()->GetVertexNonBondedEnergy(pvertex));
        // if (debug_count_old <= 10 && is_dna) {
        //     std::cerr << "---> [METROPOLIS DEBUG] OLD nonbonded energy = " << dE_nonbonded << std::endl;
        // }
    } else {
        // static int null_count = 0;
        // if (null_count < 5) {
        //     std::cerr << "---> [METROPOLIS DEBUG] WARNING: GetNonbondedInteractionBetweenVertices() is NULL!" << std::endl;
        //     null_count++;
        // }
    }
    
    // Debug output removed - uncomment if needed
    
    // --- obtaining global variables that can change by the move. Note, this is not the total volume, only the one that can change.
     double old_Tvolume = 0;
     double old_Tarea = 0;
     double old_Tcurvature = 0;
     double new_Tvolume = 0;
     double new_Tarea = 0;
     double new_Tcurvature = 0;
//--->
    if (!is_dna_vertex && m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateAVertexRingContributionToGlobalVariables(pvertex, old_Tvolume, old_Tarea, old_Tcurvature);
    }
    //---> for now, only active nematic force: ForceonVerticesfromInclusions
    Vec3D Dx(dx,dy,dz);
    double dE_force_from_inc = 0;
    double dE_force_from_vector_fields = 0;
    double dE_force_on_vertex = 0;
    
    if (!is_dna_vertex) {
        dE_force_from_inc  = m_pState->GetForceonVerticesfromInclusions()->Energy_of_Force(pvertex, Dx);
        dE_force_from_vector_fields  = m_pState->GetForceonVerticesfromVectorFields()->Energy_of_Force(pvertex, Dx);
        dE_force_on_vertex  = m_pState->GetForceonVertices()->Energy_of_Force(pvertex, Dx);
    }

//----> Move the vertex;
        pvertex->PositionPlus(dx,dy,dz);
    //--- update triangles normal (only for membrane vertices)
    if (!is_dna_vertex) {
    for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
        (*it)->UpdateNormal_Area(m_pBox);
    }
    //  check new faces angles, if bad, reverse the trinagles
    if(!CheckFacesAfterAVertexMove(pvertex)){
        for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
        }
        pvertex->PositionPlus(-dx,-dy,-dz);
        return false;
        }
    }
//---->
    //--> calculate edge shape operator (only for membrane vertices)
    if (!is_dna_vertex) {
    for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
       // (*it)->UpdateNormal();
        //  (*it)->GetNeighborLink1()->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
        (*it)->GetNeighborLink1()->UpdateShapeOperator(m_pBox);
    }
    //-- we need this to make sure all the links connected to this v is updated
    if(pvertex->GetVertexType() == 1){
        pvertex->GetPrecedingEdgeLink()->UpdateEdgeVector(m_pBox);
    }
    // --> calculate vertex shape operator
    (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(pvertex);
    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
            // Skip DNA beads - use GetBonds() instead of GetVLinkList() to avoid hang
            if (!(*it)->GetBonds().empty()) {
                continue;  // Skip DNA beads
            }
        (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(*it);
    }
        //---> calculate new energies (only for membrane vertices)
    new_energy = (m_pState->GetEnergyCalculator())->SingleVertexEnergy(pvertex);
    new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(pvertex);

    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
            // Skip DNA beads - use GetBonds() instead of GetVLinkList() to avoid hang
            if (!(*it)->GetBonds().empty()) {
                continue;  // Skip DNA beads
            }
        new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(*it);
        new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(*it);
        }
    } else {
        // DNA beads: no membrane energy, only bond energy
        new_energy = 0;
    }
    //-- interaction energy should be calculated here (only for membrane vertices)
    if (!is_dna_vertex) {
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
        
        if(pvertex->GetNumberOfVF() != 0 ){
            for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++){
                
                new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
        }
            }
        }
    }
    //---> get energy for ApplyConstraintBetweenGroups
    double dE_Cgroup = 0;
    double dE_volume = 0;
    double dE_t_area = 0;
    double dE_g_curv = 0;

    if (!is_dna_vertex) {
        dE_Cgroup = m_pState->GetApplyConstraintBetweenGroups()->CalculateEnergyChange(pvertex, Dx);
//---> new global variables
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateAVertexRingContributionToGlobalVariables(pvertex, new_Tvolume, new_Tarea, new_Tcurvature);
    }
    //---> energy change of global variables
        dE_volume =  m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
        dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
        
        // Debug: check curvature values before calculation
        double dTarea = new_Tarea - old_Tarea;
        double dTcurvature = new_Tcurvature - old_Tcurvature;
        static int curv_debug_count = 0;
        if (std::isnan(dTcurvature) || std::isnan(new_Tcurvature) || std::isnan(old_Tcurvature)) {
            #pragma omp critical
            {
                curv_debug_count++;
                if (curv_debug_count <= 3) {
                    std::cout << "Debug: Curvature NaN! old_Tcurvature=" << old_Tcurvature 
                              << ", new_Tcurvature=" << new_Tcurvature 
                              << ", dTcurvature=" << dTcurvature << std::endl;
                }
            }
        }
        
        // Calculate curvature energy change
        dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(dTarea, dTcurvature);
        
        // If result is NaN, set to 0 to avoid NaN propagation
        // This can happen if curvature calculation produces invalid values
        if (std::isnan(dE_g_curv)) {
            static int nan_fix_count = 0;
            #pragma omp critical
            {
                nan_fix_count++;
                if (nan_fix_count <= 3) {
                    std::cout << "Debug: dE_g_curv was NaN, setting to 0. dTarea=" << dTarea 
                              << ", dTcurvature=" << dTcurvature << std::endl;
                }
            }
            dE_g_curv = 0.0;
        }
    }
    
    //--> bond and nonbonded energy after move (add to existing values)
    bond_energy += pvertex->GetBondedEnergyOfVertex(); // New bonded energy (bonds + angles)
    // Two-way repulsion: both membrane vertices and DNA beads feel repulsion
    if (m_pState->GetNonbondedInteractionBetweenVertices() != nullptr) {
        // static int debug_count_new = 0;
        // bool is_dna = !pvertex->GetBonds().empty();
        // if (debug_count_new < 10 && is_dna) {
        //     std::cerr << "---> [METROPOLIS DEBUG] Calculating NEW nonbonded energy for DNA bead ID=" 
        //               << pvertex->GetVID() << std::endl;
        //     debug_count_new++;
        // }
        double new_nonbonded = m_pState->GetNonbondedInteractionBetweenVertices()->GetVertexNonBondedEnergy(pvertex);
        dE_nonbonded += new_nonbonded; // New nonbonded energy
        // if (debug_count_new <= 10 && is_dna) {
        //     std::cerr << "---> [METROPOLIS DEBUG] NEW nonbonded energy = " << new_nonbonded 
        //               << ", total dE_nonbonded = " << dE_nonbonded << std::endl;
        // }
    }
    
    //--> only elatsic energy
    double diff_energy = new_energy - old_energy;
            changed_en = diff_energy;
    //std::cout<<diff_energy<<" dif en \n";
    //--> sum of all the energies
    // NOTE: bond_energy and dE_nonbonded were added for DNA support, but original only had:
    // double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_on_vertex + dE_force_from_inc + dE_force_from_vector_fields + dE_volume + dE_t_area + dE_g_curv;
    // Keeping bond_energy and dE_nonbonded for DNA functionality
    double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_on_vertex + dE_force_from_inc + dE_force_from_vector_fields + dE_volume + dE_t_area + dE_g_curv + bond_energy + dE_nonbonded;
    
    double U = m_Beta * tot_diff_energy - m_DBeta;
    
    //---> accept or reject the move
    if(U <= 0 || exp(-U) > temp ) {
        // move is accepted
        // NOTE: Energy is updated via reduction in the caller (EvolveOneStep) using changed_en
        // changed_en is set to diff_energy (elastic only) above, matching original behavior
        // AddToTotalEnergy is commented out here (as in original) - energy updated via reduction
        //---> if vertex is out of the voxel, update its voxel
        // Update voxels for both membrane vertices and DNA beads
        // DNA beads need voxel updates for efficient neighbor searches in repulsion calculations
        if (pvertex->GetVoxel() != nullptr && !pvertex->CheckVoxel()){
            // Try to update voxel - this may fail if vertex moved more than one voxel
            // In that case, the voxel will be stale but the fallback search will catch it
            pvertex->UpdateVoxelAfterAVertexMove();
        }
        //---> ApplyConstraintBetweenGroups
        if (!is_dna_vertex) {
        m_pState->GetApplyConstraintBetweenGroups()->AcceptMove();
        //---> global variables
        if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
            m_pState->GetVAHGlobalMeshProperties()->Add2Volume(new_Tvolume - old_Tvolume);
            m_pState->GetVAHGlobalMeshProperties()->Add2TotalArea(new_Tarea - old_Tarea);
            m_pState->GetVAHGlobalMeshProperties()->Add2GlobalCurvature(new_Tcurvature - old_Tcurvature);
            }
        }
        return true;
    }
    else {
        // move is rejected
        // Debug output removed - uncomment if needed
//---> reverse the changes that has been made to the system
        //---> reverse the changes that has been made to the system
        changed_en = 0;
        //---> reverse the triangles (only for membrane vertices)
        if (is_dna_vertex) {
            // DNA beads: just reverse position
            pvertex->ReverseConstantMesh_Copy();
            return false;
        }
        if (!is_dna_vertex) {
        for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
        }
        //---> reverse the links
        for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
            (*it)->GetNeighborLink1()->ReverseConstantMesh_Copy();
        }
        //--> the shape operator of these links has not been affected, therefore we only update the interaction energy
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
            (*it)->Reverse_VFInteractionEnergy();
        }
        //-- we need this to make sure all the links connected to this v is updated
        if(pvertex->GetVertexType() == 1){
            pvertex->GetPrecedingEdgeLink()->ReverseConstantMesh_Copy();
            }
        }
        //---> reverse the vertices
        pvertex->ReverseConstantMesh_Copy();
        if (!is_dna_vertex) {
        pvertex->Reverse_VFsBindingEnergy();
        for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
                // Skip DNA beads - use GetBonds() instead of GetVLinkList() to avoid hang
                if (!(*it)->GetBonds().empty()) {
                    continue;  // Skip DNA beads
                }
            (*it)->ReverseConstantMesh_Copy();
            (*it)->Reverse_VFsBindingEnergy();
            }
        }

        return false;
     }
    return true;
}
//---> this does not check the angle of the faces. Because this should be done after the move:
//waste of calculation if we do ith before the move. Unless, we store the values. That also not good because move could get rejected.
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2){
//--->  vertex new position if get accepted
    
        // Safety check: ensure pvertex and m_pBox are valid
        if (pvertex == nullptr || m_pBox == nullptr) {
            return false;
        }
    
        double new_x = pvertex->GetXPos() + dx;
        double new_y = pvertex->GetYPos() + dy;
        double new_z = pvertex->GetZPos() + dz;
        //-- if the adding crosses the box
        if (new_x >= (*m_pBox)(0)) {
            new_x = new_x - (*m_pBox)(0);
        } else if (new_x < 0) {
            new_x = new_x + (*m_pBox)(0);
        }
        if (new_y >= (*m_pBox)(1)) {
            new_y = new_y - (*m_pBox)(1);
        } else if (new_y < 0) {
            new_y = new_y + (*m_pBox)(1);
        }
        if (new_z >= (*m_pBox)(2)) {
            new_z = new_z - (*m_pBox)(2);
        } else if (new_z < 0) {
            new_z = new_z + (*m_pBox)(2);
        }
  
//--->  let check the distances with the nighbours
    // For vertices with links (membrane vertices), check link-based neighbors
    // For vertices without links (DNA beads), skip this check as they're not constrained by membrane edge lengths
    std::vector <vertex *> npvertex = pvertex->GetVNeighbourVertex();
    // CRITICAL FIX: Use GetBonds() instead of GetVLinkList() to identify DNA beads
    bool is_dna = !pvertex->GetBonds().empty();
    if (!is_dna) {
        // This is a membrane vertex - check distances to link-connected neighbors
    for (std::vector<vertex *>::iterator it = npvertex.begin() ; it != npvertex.end(); ++it){
            // Only check neighbors that are also membrane vertices (have links, no bonds)
            // DNA beads (with bonds) shouldn't constrain membrane vertex moves
            if ((*it)->GetBonds().empty()) {
        double dist2 = pvertex->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it);
                if(dist2 < mindist2 || dist2 > maxdist2) {
            return false;
                }
            }
        }
    }
    // For DNA beads (no links), we skip the Min_Max_Lenghts check as they're not part of the membrane mesh
//---> now check it within the voxel cells
//---> lets get the object voxel, note: the voxel could be different from the associated one as it is moving
        //-- Check if voxel is initialized (should always be, but safety check)
        // For DNA beads, skip voxel-based checks entirely as they're not constrained by membrane edge lengths
        if (is_dna) {
            return true;  // DNA beads are free to move without voxel-based constraints
        }
        
        Voxel<vertex>* pvoxel = pvertex->GetVoxel();
        if (pvoxel == nullptr || m_pBox == nullptr) {
            // Voxel or box not initialized - this shouldn't happen, but skip check if it does
            return true;
        }
        
        //-- obtain the object (vertex) new cell id, with respect to the current cell
        int NoX = pvoxel->GetXNoVoxel();
        int NoY = pvoxel->GetYNoVoxel();
        int NoZ = pvoxel->GetZNoVoxel();
    
        int new_nx = int(new_x/pvoxel->GetXSideVoxel((*m_pBox)(0)));
        int new_ny = int(new_y/pvoxel->GetYSideVoxel((*m_pBox)(1)));
        int new_nz = int(new_z/pvoxel->GetZSideVoxel((*m_pBox)(2)));
    
        int old_nx = pvoxel->GetXIndex();
        int old_ny = pvoxel->GetYIndex();
        int old_nz = pvoxel->GetZIndex();
    
        int i = Voxel<int>::Convert2LocalVoxelIndex(new_nx, old_nx, NoX);
        int j = Voxel<int>::Convert2LocalVoxelIndex(new_ny, old_ny, NoY);
        int k = Voxel<int>::Convert2LocalVoxelIndex(new_nz, old_nz, NoZ);

    
        //-- check if it has moved too far
        if(i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
            return false;
        }
  
        Voxel<vertex>* new_pvox = pvoxel->GetANeighbourCell(i, j, k);
        if (new_pvox == nullptr) {
            // Neighbor cell is null - shouldn't happen, but skip check if it does
            return true;
        }
  
        // For DNA beads (no links), completely skip the voxel-based Min_Max_Lenghts check
        // DNA beads should be free to move without constraint
        if (!is_dna) {
            // Only check constraints for membrane vertices
        for(int n=-1;n<2;n++)
        for(int m=-1;m<2;m++)
        for(int s=-1;s<2;s++){
                Voxel<vertex>* neighbor_voxel = new_pvox->GetANeighbourCell(n, m, s);
                if (neighbor_voxel == nullptr) {
                    continue;  // Skip if neighbor voxel is null
                }
                std::vector <vertex *> CV = neighbor_voxel->GetContentObjects();
            for (std::vector<vertex *>::iterator it = CV.begin() ; it != CV.end(); ++it){
                if(*it != pvertex){
                        // Only check distances to other membrane vertices (no bonds)
                        // DNA beads (with bonds) should not constrain membrane vertex moves
                        if ((*it)->GetBonds().empty()) {
                            double dist2 = pvertex->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it);
                            if(dist2 < mindist2 || dist2 > maxdist2) {
                        return false;
                            }
                        }
                    }
                }
            } ///   for(int s=-1;s<2;s++){
            }
        // DNA beads skip all Min_Max_Lenghts checks - they're free to move
    ///
    ///
    /*  */
    return true;
}
// finding the distance of the current vertex from the pv2; also considering the pbc conditions
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::CheckFacesAfterAVertexMove(vertex* p_vertex) {
    std::vector<links*> linkList = p_vertex->GetVLinkList();
    for (std::vector<links*>::iterator it = linkList.begin(); it != linkList.end(); ++it) {
        links* link = *it;
        if (!link->CheckFaceAngleWithMirrorFace(m_MinAngle) || !link->CheckFaceAngleWithNextEdgeFace(m_MinAngle)) {
            return false;
        }
    }
    return true;
}

//========
// this function can be deleted any time; it is for test cases only
double  EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::SystemEnergy()
{
    /*
    MESH* m_pMESH = m_pState->m_pMesh;
    std::vector<vertex *> ActiveV = m_pMESH->m_pActiveV;
    std::vector<triangle *> pActiveT = m_pMESH->m_pActiveT;
    std::vector<links *> mLink = m_pMESH->m_pHL;
    std::vector<links *>  pEdgeL = m_pMESH->m_pEdgeL;
    std::vector<vertex *> EdgeV  = m_pMESH->m_pEdgeV;
    double en = 0;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = mLink.begin() ; it != mLink.end(); ++it)
    {
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    
    for (std::vector<vertex *>::iterator it = ActiveV.begin() ; it != ActiveV.end(); ++it)
        (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);

    //====== edge links should be updated
    for (std::vector<links *>::iterator it = pEdgeL.begin() ; it != pEdgeL.end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = EdgeV.begin() ; it != EdgeV.end(); ++it)
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);
    
    en=m_pEnergyCalculator->TotalEnergy(ActiveV,mLink);
    //en=en+ m_pEnergyCalculator->TotalEnergy(EdgeV,pEdgeL);
   
    
    return en;
     */
    return 0;
}
std::string EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::CurrentState(){
    
    std::string state = AbstractVertexPositionIntegrator::GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(m_NumberOfMovePerStep_Surf) +" "+Nfunction::D2S(m_NumberOfMovePerStep_Edge);
    state = state +" "+ Nfunction::D2S(m_DR);
    return state;
}
std::vector<links*> EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::GetEdgesWithInteractionChange(vertex* p_vertex){
#if DEVELOPMENT_MODE == Enabled
    std::cout<<" DEVELOPMENT_MODE ID 665656: This function should be made much better\n ";
    // for example, precheck to select only links that both has includioon ...
#endif
    
    std::vector<links*> edge_with_interaction_change;
    std::vector<links *> all_temlinks;
    
    std::vector<vertex *> neighbor_vertices = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = neighbor_vertices.begin() ; it != neighbor_vertices.end(); ++it)
    {
            if((*it)->VertexOwnInclusion() || p_vertex->GetNumberOfVF() != 0 ) {  // due to vector fields
                std::vector<links *> ltem = (*it)->GetVLinkList();
                all_temlinks.insert(all_temlinks.end(), ltem.begin(), ltem.end());
                // note, this is even need it if the p_vertex vertex is not an edge vertex
                if((*it)->m_VertexType == 1){
                    all_temlinks.push_back((*it)->m_pPrecedingEdgeLink);
                }
            }
    }

    
    //-- now we remove the repeated links
    for (std::vector<links *>::iterator it = all_temlinks.begin() ; it != all_temlinks.end(); ++it)
    {
        bool Should_be_added = true;
        for (std::vector<links *>::iterator it2 = edge_with_interaction_change.begin() ; it2 != edge_with_interaction_change.end(); ++it2)
        {
            if( *it2 == *it){
                Should_be_added = false;
            }
            else if((*it2)->GetMirrorFlag() && (*it2)->GetMirrorLink() == *it){
                Should_be_added = false;
            }
        }
        if(Should_be_added == true)
            edge_with_interaction_change.push_back((*it));

    }
    
    return edge_with_interaction_change;
}



