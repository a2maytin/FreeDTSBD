#ifndef ANGLE_H
#define ANGLE_H

#include "SimDef.h"
#include "Vec3D.h"

/*
 * -----------------------------------------------------------------------------
 * Angle Class - Simulation Component
 * -----------------------------------------------------------------------------
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Date: December 2024
 *
 * Description:
 * The `angle` class represents an angle potential in a simulation. It models
 * the bending energy between three consecutive vertices (i-1, i, i+1) with
 * properties such as stiffness (k_angle) and equilibrium angle (theta0).
 * The class provides methods to calculate angle energy.
 *  E = k_angle/2 * (theta - theta0)^2
 *
 * -----------------------------------------------------------------------------
 * License:
 * Copyright (c) Weria Pezeshkian, 2024. All rights reserved.
 * -----------------------------------------------------------------------------
 */

class vertex;

class angle {
public:
    // Constructors and Destructor
    angle(int id, vertex* v1, vertex* v2, vertex* v3, double k_angle, double theta0);
    explicit angle(int id);
    ~angle();

    // Accessor Methods
    int GetID() const { return m_ID; }
    vertex* GetV1() const { return m_V1; }  // First vertex (i-1)
    vertex* GetV2() const { return m_V2; }  // Middle vertex (i)
    vertex* GetV3() const { return m_V3; }  // Third vertex (i+1)

    // Public Member Functions
    void UpdateV(vertex* v1, vertex* v2, vertex* v3);
    void UpdateAngleK(double k_angle);
    void UpdateAngleTheta0(double theta0);
    double CalculateEnergy() const;

private:
    // Member Variables
    int m_ID;
    vertex* m_V1;
    vertex* m_V2;
    vertex* m_V3;
    double m_K_Angle;  // Angle stiffness (k_angle/2 stored for efficiency)
    double m_Theta0;   // Equilibrium angle in radians
    Vec3D *m_pBox; 

};

#endif // ANGLE_H

