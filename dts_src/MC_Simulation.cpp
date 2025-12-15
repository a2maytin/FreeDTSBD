

#ifdef _OPENMP
# include <omp.h>
#endif
#include <thread>
#include <ctime>
#include <iostream>
#include <string.h>
#include "MC_Simulation.h"
#include "State.h"
#include "SimDef.h"
/*
 List of skipped function due to lack of clarity based on current state of the code
 They need to be finished before calling this a new version.
 
 1) void TimeSeriesLogInformation::WriteStateInfo(){
 
 */






/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 MC simulation class, runs mc simulation if it is defined in the input file.
 */
MC_Simulation::MC_Simulation(State *pState) : m_pState(pState) {

}
MC_Simulation::~MC_Simulation(){
    
}
void MC_Simulation::Initialize(){
    
    return;
}
bool MC_Simulation::do_Simulation(){
#if DEBUG_MODE == Enabled
    std::cout<<" do_Simulation function is starting  \n";
#endif
    
//---> Voxelize the mesh for the first time; this should be done before any calculation
    m_pState->GetVoxelization()->Voxelize(m_pState->GetMesh()->GetActiveV());
    
#if DEBUG_MODE == Enabled
    std::cout<<" system has been voxelaized  \n";
#endif

//----> checking if the mesh is good, within the bond of the simulation type. For here, it should be within
        //CheckMesh();
 
//--- before simualtion lets have a frame of the initial system
        m_pState->GetVisualization()->WriteAFrame(0);
   // time_t startTime;
   // time(&startTime);
#if DEBUG_MODE == Enabled
    std::cout<<" We have reached simulation run loop!  \n";
#endif
    std::clock_t start = std::clock();
    
#ifdef _OPENMP
    double startwall_time = omp_get_wtime();
#endif
// Your OpenMP parallel code here

    if(!m_pState->GetMesh()->CheckMesh(m_MinLength2, m_MaxLength2, m_MinAngle, m_pState->GetVoxelization())){
        
        std::cout << "---> error: The  mesh quality is insufficient for running a simulation.\n";
        exit(0);
    }

    
    

    std::cout<<"------>   Simulation will be performed from "<<m_Initial_Step<<" to "<<m_Final_Step<<" steps\n";
for (int step = m_Initial_Step; step <= m_Final_Step; step++){
        
//---> System rotation to randomize multithreading bias direction (MUST happen before output)
    std::vector<vertex*>& all_vertices = m_pState->GetMesh()->GetActiveV();
    
    // Initialize on first step if rotation is enabled
    if (m_pState->GetSystemRotation()->IsEnabled() && step == m_Initial_Step) {
        // Calculate COM and center system at box center
        double xcm = 0.0, ycm = 0.0, zcm = 0.0;
        for (auto* v : all_vertices) {
            xcm += v->GetVXPos();
            ycm += v->GetVYPos();
            zcm += v->GetVZPos();
        }
        double n_vertices = static_cast<double>(all_vertices.size());
        Vec3D com(xcm/n_vertices, ycm/n_vertices, zcm/n_vertices);
        
        // Get box center
        Vec3D *pBox = m_pState->GetMesh()->GetBox();
        Vec3D box_center((*pBox)(0)/2.0, (*pBox)(1)/2.0, (*pBox)(2)/2.0);
        
        // Initialize rotation system with box center
        m_pState->GetSystemRotation()->Initialize(box_center);
        
        // Center the system at box center (like CenterMesh)
        for (auto* v : all_vertices) {
            v->UpdateVXPos(v->GetVXPos() - com(0) + box_center(0));
            v->UpdateVYPos(v->GetVYPos() - com(1) + box_center(1));
            v->UpdateVZPos(v->GetVZPos() - com(2) + box_center(2));
        }
        
        // Update geometry after centering
        std::vector<triangle*>& all_triangles = m_pState->GetMesh()->GetActiveT();
        std::vector<links*>& all_links = m_pState->GetMesh()->GetActiveL();
        for (auto* tri : all_triangles) {
            if (tri != nullptr) {
                tri->UpdateNormal_Area(pBox);
            }
        }
        for (auto* link : all_links) {
            if (link != nullptr) {
                link->UpdateShapeOperator(pBox);
            }
        }
        
        // Re-voxelize after centering
        m_pState->GetVoxelization()->ReassignMembersToVoxels(all_vertices);
        
        std::cout << "---> System rotation initialized: centered system at box center\n";
    }
    
    if(m_pState->GetSystemRotation()->UpdateRotation(step)){
        // Rotation was applied - rotate all vertices (system is already centered at origin)
        std::vector<triangle*>& all_triangles = m_pState->GetMesh()->GetActiveT();
        std::vector<links*>& all_links = m_pState->GetMesh()->GetActiveL();
        std::vector<inclusion*>& all_inclusions = m_pState->GetMesh()->GetInclusion();
        Vec3D *pBox = m_pState->GetMesh()->GetBox();
        
        // Rotate all vertices: move to origin, apply rotation, move back to box center
        Vec3D box_center = m_pState->GetSystemRotation()->GetBoxCenter();
        
        // Debug: Track bead with ID 1000 (if it exists)
        vertex* debug_vertex = nullptr;
        Vec3D debug_pos_before, debug_pos_after_rotation;
        for (auto* v : all_vertices) {
            if (v->GetVID() == 1000) {
                debug_vertex = v;
                debug_pos_before = Vec3D(v->GetVXPos(), v->GetVYPos(), v->GetVZPos());
                break;
            }
        }
        
        for (auto* v : all_vertices) {
            // Move to origin (relative to box center)
            Vec3D pos(v->GetVXPos() - box_center(0), v->GetVYPos() - box_center(1), v->GetVZPos() - box_center(2));
            // Apply NEW rotation (not cumulative) - this is the rotation that was just generated
            m_pState->GetSystemRotation()->ApplyNewRotationToVertex(pos);
            // Move back to box center
            v->UpdateVXPos(pos(0) + box_center(0));
            v->UpdateVYPos(pos(1) + box_center(1));
            v->UpdateVZPos(pos(2) + box_center(2));
            
            // Debug: Capture position after rotation for bead 1000
            if (v == debug_vertex) {
                debug_pos_after_rotation = Vec3D(v->GetVXPos(), v->GetVYPos(), v->GetVZPos());
            }
        }
        
        // Debug output
        if (debug_vertex != nullptr) {
            // Calculate what the position should be after rotation
            Vec3D pos_relative = debug_pos_before - box_center;
            m_pState->GetSystemRotation()->ApplyNewRotationToVertex(pos_relative);
            Vec3D expected_after = pos_relative + box_center;
            
            Vec3D diff = debug_pos_after_rotation - expected_after;
            double error = sqrt(diff(0)*diff(0) + diff(1)*diff(1) + diff(2)*diff(2));
            
            std::cout << "---> Rotation at step " << step << " - Bead ID 1000:\n";
            std::cout << "     Before rotation: (" << debug_pos_before(0) << ", " << debug_pos_before(1) << ", " << debug_pos_before(2) << ")\n";
            std::cout << "     After rotation:  (" << debug_pos_after_rotation(0) << ", " << debug_pos_after_rotation(1) << ", " << debug_pos_after_rotation(2) << ")\n";
            std::cout << "     Expected after:  (" << expected_after(0) << ", " << expected_after(1) << ", " << expected_after(2) << ")\n";
            std::cout << "     Rotation error:  " << error << "\n";
        }
        
        // Update all triangle normals and areas after rotation
        for (auto* tri : all_triangles) {
            if (tri != nullptr) {
                tri->UpdateNormal_Area(pBox);
            }
        }
        
        // Update all shape operators after rotation
        for (auto* link : all_links) {
            if (link != nullptr) {
                link->UpdateShapeOperator(pBox);
            }
        }
        
        // Reinitialize curvature calculations (updates vertex normals, areas, and transfer matrices)
        m_pState->GetCurvatureCalculator()->Initialize();
        
        // Now rotate inclusions and their directions (after transfer matrices are updated)
        for (auto* inc : all_inclusions) {
            if (inc == nullptr) continue;
            vertex* v = inc->Getvertex();
            if (v == nullptr) continue;
            
            // Rotate inclusion global direction
            Vec3D gdir = inc->GetGDirection();
            if (gdir(0) != 0.0 || gdir(1) != 0.0 || gdir(2) != 0.0) {
                // Rotate global direction using NEW rotation (not cumulative)
                m_pState->GetSystemRotation()->ApplyNewRotationToVertex(gdir);
                inc->UpdateGlobalDirection(gdir);
                // Update local direction from rotated global (transfer matrices are now updated)
                inc->UpdateLocalDirectionFromGlobal();
            }
        }
        
        // Rotate vector fields (after transfer matrices are updated)
        int nvf = m_pState->GetMesh()->GetNoVFPerVertex();
        if (nvf > 0) {
            for (auto* v : all_vertices) {
                for (int i = 0; i < nvf; i++) {
                    // Get global direction, rotate it, update
                    Vec3D gdir = v->GetVectorField(i)->GetGDirection();
                    if (gdir(0) != 0.0 || gdir(1) != 0.0 || gdir(2) != 0.0) {
                        // Rotate using NEW rotation (not cumulative)
                        m_pState->GetSystemRotation()->ApplyNewRotationToVertex(gdir);
                        v->GetVectorField(i)->UpdateGlobalDirection(gdir);
                        // Update local direction from rotated global (transfer matrices are now updated)
                        v->UpdateVFLocalDirectionFromGlobalDirection();
                    }
                }
            }
        }
        
        // Re-voxelize after rotation
        m_pState->GetVoxelization()->ReassignMembersToVoxels(all_vertices);
        std::cout << "---> Applied system rotation at step " << step << " to randomize bias direction\n";
    }

//----> write files (after rotation, so output can apply inverse rotation)
        //--- write visulaization frame
        m_pState->GetVisualization()->WriteAFrame(step);
        //--- write non-binary trejectory e.g., tsi, tsg
        m_pState->GetNonbinaryTrajectory()->WriteAFrame(step);
        //--- write binary trejectory e.g., bts
        m_pState->GetBinaryTrajectory()->WriteAFrame(step);
        //--- write into time seri file, e.g., energy, volume ...
        m_pState->GetTimeSeriesDataOutput()->WriteTimeSeriesDataOutput(step);
        //--- write check point for the state
        m_pState->GetRestart()->UpdateRestartState(step, m_pState->GetVertexPositionUpdate()->GetDR(), m_pState->GetDynamicBox()->GetDR());
    
//---> centering the simulation box
    if(m_CenteringFrequently != 0 && step%m_CenteringFrequently == 0){
        m_pState->GetMesh()->CenterMesh();
        m_pState->GetVoxelization()->ReassignMembersToVoxels(m_pState->GetMesh()->GetActiveV());
    } // [ if(GetBoxCentering()!=0 && step%GetBoxCentering()==0)]

//---> Run standard Integrators
        //--- run the vertex position update
        m_pState->GetVertexPositionUpdate()->EvolveOneStep(step); // we may need the final step as well to check if the update of move size should be done
        //--- run the link flip update
        m_pState->GetAlexanderMove()->EvolveOneStep(step);
        //--- run the inclusion update
        m_pState->GetInclusionPoseUpdate()->EvolveOneStep(step);
        //--- run vector fields
        m_pState->GetVectorFieldsRotationUpdate()->EvolveOneStep(step);

//----> Run the supplementary integrators
       //--- update the box side
         m_pState->GetDynamicBox()->ChangeBoxSize(step); // we may need the final step as well to check if the update of move size should be done
        //--- update edge of mesh open edge
         m_pState->GetOpenEdgeEvolution()->Move(step);
        //--- update the mesh topology
        m_pState->GetDynamicTopology()->MCMove(step);
        //---- convert inclusions
        m_pState->GetInclusionConversion()->Exchange(step);
    
        //---- NonequilibriumCommands
        m_pState->GetNonequilibriumCommands()->Run(step);

//----> print info about the simulation, e.g., rate,
   // time_t currentTime;
   // time(&currentTime);
    if(!CheckMesh(step)){
        std::cout<<"---> error, the mesh does not meet the requirment for MC sim \n";
    }
    if (step%100 == 0) {
        PrintRate(step, true, true);
    }

} // for(int step=GetInitialStep(); step<GetFinalStep(); step++)
    std::clock_t end = std::clock();
    

    
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    std::cout<<"---- Simulation has ended ----\n";
    std::cout<<" The run took: "<<Nfunction::ConvertSecond2Time(elapsed_secs)<<"thread time and  ";
#ifdef _OPENMP
    double endwall_time = omp_get_wtime();
    endwall_time = endwall_time - startwall_time;
    std::cout<<Nfunction::ConvertSecond2Time(endwall_time)<<" wall time \n";
#endif

    m_pState->GetCurvatureCalculator()->Initialize();
    double Final_energy = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
    double energy_leak = Final_energy - m_pState->GetEnergyCalculator()->GetEnergy();
    std::cout << std::fixed << std::setprecision(4);
    if(fabs(energy_leak) > 0.0001){
        
        std::cout<<"---> possible source of code error: energy leak... "<<energy_leak<<" with real energy of "<<Final_energy<<"  and stored energy of "<<m_pState->GetEnergyCalculator()->GetEnergy()<<"\n";
    }

    
    
    double vol = 0;
    double g_c = 0;
    double t_a = 0;
    m_pState->GetVAHGlobalMeshProperties()->CalculateGlobalVariables(vol,t_a,g_c);
    vol -= m_pState->GetVAHGlobalMeshProperties()->GetTotalVolume();
    g_c -= m_pState->GetVAHGlobalMeshProperties()->GetTotalMeanCurvature();
    t_a -= m_pState->GetVAHGlobalMeshProperties()->GetTotalArea();

    if (m_pState->GetVAHGlobalMeshProperties()->VolumeIsActive())
    if(fabs(vol) > 0.0001 ){
        std::cout<<fabs(vol)<<" volume leak\n";
    }
    if (m_pState->GetVAHGlobalMeshProperties()->GlobalCurvatureIsActive())
    if(fabs(g_c) > 0.0001)
    {
        std::cout<<fabs(g_c)<<" global curvature leak\n";
    }
    if (m_pState->GetVAHGlobalMeshProperties()->AreaIsActive())
    if(fabs(t_a) > 0.0001)
    {
        std::cout<<fabs(t_a)<<" total area leak\n";
    }
        
    return true;
}
void  MC_Simulation::PrintRate(int step, bool clean, bool clear){
    
    double  vmove_rate =  100 * (m_pState->GetVertexPositionUpdate()->GetAcceptanceRate(clean));
    double  emove_rate =  100 * (m_pState->GetAlexanderMove()->GetAcceptanceRate(clean));
    double  imove_rate =  100 * (m_pState->GetInclusionPoseUpdate()->GetAcceptanceRate(clean));
    double  bmove_rate =  100 * (m_pState->GetDynamicBox()->GetAcceptanceRate(clean));
    double  vfmove_rate = 100 * (m_pState->GetVectorFieldsRotationUpdate()->GetAcceptanceRate(clean));

    std::cout<<"Step = "<<step<<"/"<<m_Final_Step<<std::flush;
    std::cout << std::fixed << std::setprecision(1);
    std::cout<<" Rates: "<<std::flush;
    std::cout<<" vertex move = "<<vmove_rate<<"%"<<std::flush;
    std::cout<<"; alexander move = "<<emove_rate<<"%"<<std::flush;
    if(m_pState->GetMesh()->GetInclusion().size() != 0){
        std::cout<<"; inclusion move = "<<imove_rate<<"%"<<std::flush;
    }
    if(m_pState->GetMesh()->GetNoVFPerVertex() != 0){
        std::cout<<"; vector fields move = "<<vfmove_rate<<"%"<<std::flush;
    }
    if(m_pState->GetDynamicBox()->GetDerivedDefaultReadName() != "No"){
        std::cout<<"; Box Move = "<<bmove_rate<<"%"<<std::flush;
    }
    if(clear){
        std::cout << '\r';
        std::cout << "\033[K";
    }
}
std::string MC_Simulation::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    state = state + "\n Min_Max_Lenghts = "+Nfunction::D2S(m_MinLength2)+" "+Nfunction::D2S(m_MaxLength2);
    state = state + "\n MinfaceAngle = "+Nfunction::D2S(m_MinAngle);
    state = state + "\n Temprature = "+Nfunction::D2S(m_Beta)+" "+Nfunction::D2S(m_DBeta);
    state = state + "\n Box_Centering_F = "+Nfunction::D2S(m_CenteringFrequently);
    state = state + "\n Set_Steps = "+Nfunction::D2S(m_Initial_Step)+" "+Nfunction::D2S(m_Final_Step);
    
    return state;
}
bool MC_Simulation::CheckMesh(int step){
    
    if(m_CheckMeshFrequently == 0 || step%m_CheckMeshFrequently == 0){
        return true;
    }
    
    // this is only the  edges
    const std::vector<links *>& all_links = m_pState->GetMesh()->GetActiveL();
    for (std::vector<links *>::const_iterator it = all_links.begin() ; it != all_links.end(); ++it){

        vertex *p_v1 = (*it)->GetV1();
        vertex *p_v2 = (*it)->GetV2();
        
        // Skip if either vertex is a bonded (linearly connected untriangulated) vertex
        // Bonded vertices shouldn't have links, but we check bonds first to avoid calling GetVLinkList which can hang
        bool is_bonded_v1 = !p_v1->GetBonds().empty();
        bool is_bonded_v2 = !p_v2->GetBonds().empty();
        if (is_bonded_v1 || is_bonded_v2) continue;
        
        double dist2 = p_v1->SquareDistanceFromAVertex(p_v2);
        if(dist2 < m_MinLength2 || dist2 > m_MaxLength2){
            
            return false;
        }
    }
   
    // all vertices should also have a distance larger then sqrt(m_MinLength2)


    
    return true;
}
