#ifndef VertexHarmonicBounds_H
#define VertexHarmonicBounds_H

#include "SimDef.h"
#include "Vec3D.h"
#include "bond.h"
#include "angle.h"


class MESH;
class VertexHarmonicBounds {

public:

    VertexHarmonicBounds();
    ~VertexHarmonicBounds();
    inline std::vector<bond *> GetBonds()   {return m_pVertexBond;}
    inline std::vector<angle *> GetAngles()   {return m_pVertexAngle;}
    double GetBondEnergyOfVertex();
    double GetAngleEnergyOfVertex();
    double GetBondedEnergyOfVertex();  // Returns bond + angle energy
    void AddBondToList(bond* b);
    void AddAngleToList(angle* a);

protected:
    std::vector<bond *> m_pVertexBond;
    std::vector<angle *> m_pVertexAngle;


};

#endif // VertexHarmonicBounds_H
