


#include "CouplingGlobalCurvatureToHarmonicPotential.h"
#include "Nfunction.h"
#include "State.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This object affects the simulation results in every step. It is a global coupling, penalize the deviations of total mean curvature  from a target value
 It can also model the energy associated with the differnce in the inner and outer monolayer area
 The energy coupling is as below
 E=k/(2A)*(sum(2h-gc0)av)^2
 E=k/(2Ah^2)*(DA-DA0)^2
 */
CouplingGlobalCurvatureToHarmonicPotential::CouplingGlobalCurvatureToHarmonicPotential(VAHGlobalMeshProperties *VHA,  double G_kappa, double GlobalC0) :
                AbstractGlobalCurvature(VHA),
                m_K(G_kappa/2.0),
                m_gC0(GlobalC0) {

}

CouplingGlobalCurvatureToHarmonicPotential::~CouplingGlobalCurvatureToHarmonicPotential()
{
    
}
// A function to initialize total area, total curvature and energy asscoiated with the coupling
void CouplingGlobalCurvatureToHarmonicPotential::Initialize(State* pState) {
    m_pState = pState;
    m_TotalCurvature = 0;
    m_TotalArea = 0;
    m_CalculatedGlobalVariable = true;
    std::vector<vertex *> all_vertex = m_pState->GetMesh()->GetActiveV();
    for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {
        double area = (*it)->GetArea();
        double curv = (*it)->GetP1Curvature() + (*it)->GetP2Curvature();

        m_TotalArea += area;
        m_TotalCurvature += curv*area;
    }


    return;
}
double CouplingGlobalCurvatureToHarmonicPotential::GetCouplingEnergy(){
    
    double dh=(m_TotalCurvature-m_gC0*m_TotalArea);
    
    return m_K/(m_TotalArea)*dh*dh;
}
// when a vertex moves or a link flips or any other changes, we can see the changes in the total area and total mean curavture:
// this function can tell us how much this changes cost energy but it does not update the change since it can be rejected.
double CouplingGlobalCurvatureToHarmonicPotential::CalculateEnergyChange(double D_area, double D_curvature)
{
    double dh2 = m_TotalCurvature+D_curvature-m_gC0*(m_TotalArea+D_area);
    double dh1 = m_TotalCurvature-m_gC0*m_TotalArea;
    
    // Avoid division by zero or very small numbers that could cause NaN
    // If the new area is too small, use the old area to avoid numerical instability
    double new_area = m_TotalArea + D_area;
    double old_area = m_TotalArea;
    
    // Use a small threshold to avoid numerical issues
    const double MIN_AREA = 1e-10;
    if (new_area < MIN_AREA) {
        new_area = old_area;  // Fall back to old area if new area is too small
    }
    if (old_area < MIN_AREA) {
        // If old area is also too small, return 0 (shouldn't happen in normal simulations)
        return 0.0;
    }
    
    dh2 = dh2 * dh2 / new_area;
    dh1 = dh1 * dh1 / old_area;

    double de = dh2 - dh1;
    
    // Check for NaN or infinity before returning
    if (std::isnan(de) || std::isinf(de)) {
        return 0.0;  // Return 0 if calculation produces invalid result
    }
    
    return m_K * de;
}
std::string CouplingGlobalCurvatureToHarmonicPotential::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " " + Nfunction::D2S(2*m_K) + " " + Nfunction::D2S(m_gC0);
    return state;
}

