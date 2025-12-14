#if !defined(AFX_HarmonicBondsList_H_334B21B8_INCLUDED_)
#define AFX_HarmonicBondsList_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include <set>
#include "vertex.h"
#include "bond.h"
#include "angle.h"
#include "AbstractBondedPotentialBetweenVertices.h"
class State;
class HarmonicBondsList : public AbstractBondedPotentialBetweenVertices {

public:
    HarmonicBondsList(State *pState, std::string filename);
    ~HarmonicBondsList();
    void Initialize();
  //  double GetBondEnergyOfVertex(vertex *pvertex);
    double GetTotalEnergy();
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "HarmonicBondsList";}
    inline static std::string GetDefaultReadName() {return "HarmonicBondsList";}
    
private:
    std::string m_FileName;
    std::vector<bond> m_AllBonds;
    std::vector<angle*> m_AllAngles;  // Use pointers to avoid invalidation
    std::set<std::pair<int, int> > m_SMCBonds;  // Set of SMC bond pairs (vid1, vid2) to exclude from angle creation

    State *m_pState;
    
    // Helper functions
    void LoadAnglesFromFile(const std::string& angleFileName);
    void CreateAnglesFromBonds(double k_angle, double theta0);  // Generate angles from bonds (excluding SMC)
};


#endif
