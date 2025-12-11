#include <stdio.h>
#include "VertexHarmonicBounds.h"
#include "MESH.h"

VertexHarmonicBounds::VertexHarmonicBounds() {
    m_pVertexBond.clear();
    m_pVertexAngle.clear();
}

VertexHarmonicBounds::~VertexHarmonicBounds() {

}
void VertexHarmonicBounds::AddBondToList(bond *b){
    
    m_pVertexBond.push_back(b);
    return;
}
void VertexHarmonicBounds::AddAngleToList(angle *a){
    
    m_pVertexAngle.push_back(a);
    return;
}
double VertexHarmonicBounds::GetBondEnergyOfVertex(){
    
    double en = 0;
    
    for (std::vector<bond*>::iterator it = m_pVertexBond.begin() ; it != m_pVertexBond.end(); ++it){
        en += (*it)->CalculateEnergy();
    }
    
    return en;
}
double VertexHarmonicBounds::GetAngleEnergyOfVertex(){
    
    double en = 0;
    
    for (std::vector<angle*>::iterator it = m_pVertexAngle.begin() ; it != m_pVertexAngle.end(); ++it){
        en += (*it)->CalculateEnergy();
    }
    
    return en;
}
double VertexHarmonicBounds::GetBondedEnergyOfVertex(){
    
    return GetBondEnergyOfVertex() + GetAngleEnergyOfVertex();
}
