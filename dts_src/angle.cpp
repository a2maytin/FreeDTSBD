#include "angle.h"
#include "vertex.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <iostream>
#include <limits>

// Constructor with initialization of angle parameters
angle::angle(int id, vertex* v1, vertex* v2, vertex* v3, double k_angle, double theta0)
    : m_ID(id), m_V1(v1), m_V2(v2), m_V3(v3), m_K_Angle(k_angle / 2), m_Theta0(theta0) {
    if (m_V1 != nullptr) {
        m_pBox = m_V1->GetBox();
    } else {
        m_pBox = nullptr;
    }
}

// Constructor for uninitialized angle
angle::angle(int id)
    : m_ID(id), m_V1(nullptr), m_V2(nullptr), m_V3(nullptr), m_pBox(nullptr), m_K_Angle(0.0), m_Theta0(0.0) {}

// Destructor
angle::~angle() = default;

// Update vertices associated with the angle
void angle::UpdateV(vertex* v1, vertex* v2, vertex* v3) {
    m_V1 = v1;
    m_V2 = v2;
    m_V3 = v3;
}

// Update the angle stiffness
void angle::UpdateAngleK(double k_angle) {
    m_K_Angle = k_angle / 2;
}

// Update the angle's equilibrium angle
void angle::UpdateAngleTheta0(double theta0) {
    m_Theta0 = theta0;
}

// Calculate the angle's energy
double angle::CalculateEnergy() const {
    // Safety checks
    if (m_V1 == nullptr || m_V2 == nullptr || m_V3 == nullptr || m_pBox == nullptr) {
        return 0.0;
    }
    
    // Get positions of the three vertices
    Vec3D pos1 = m_V1->GetPos();
    Vec3D pos2 = m_V2->GetPos();
    Vec3D pos3 = m_V3->GetPos();
    
    // Calculate bond vectors (for polymer bending angle)
    // For three consecutive beads i-1, i, i+1:
    // bond1 = vector from i-1 to i (pos2 - pos1)
    // bond2 = vector from i to i+1 (pos3 - pos2)
    // The bending angle is the angle between these bond vectors
    // For a straight chain, bond1 and bond2 point in the same direction, so angle = 0
    Vec3D bond1 = pos2 - pos1;  // Vector from v1 to v2
    Vec3D bond2 = pos3 - pos2;  // Vector from v2 to v3
    
    // Apply periodic boundary conditions to bond vectors
    for (int i = 0; i < 3; i++) {
        if (fabs(bond1(i)) > (*m_pBox)(i) / 2.0) {
            if (bond1(i) < 0)
                bond1(i) = (*m_pBox)(i) + bond1(i);
            else if (bond1(i) > 0)
                bond1(i) = bond1(i) - (*m_pBox)(i);
        }
        if (fabs(bond2(i)) > (*m_pBox)(i) / 2.0) {
            if (bond2(i) < 0)
                bond2(i) = (*m_pBox)(i) + bond2(i);
            else if (bond2(i) > 0)
                bond2(i) = bond2(i) - (*m_pBox)(i);
        }
    }
    
    // Calculate the bending angle between the two bond vectors
    double norm1 = bond1.norm();
    double norm2 = bond2.norm();
    
    if (norm1 < 1e-10 || norm2 < 1e-10) {
        // Degenerate case: vectors are too short
        return 0.0;
    }
    
    // Calculate cosine of the bending angle
    // For a straight chain, bond1 and bond2 point in same direction, so cos_theta = 1, theta = 0
    double cos_theta = Vec3D::dot(bond1, bond2) / (norm1 * norm2);
    
    // Clamp to [-1, 1] to avoid numerical issues
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    
    // Calculate the angle in radians
    double theta = acos(cos_theta);
    
    // Calculate energy: E = k_angle/2 * (theta - theta0)^2
    double delta_theta = theta - m_Theta0;
    double energy = m_K_Angle * delta_theta * delta_theta;
    
    // Safety check: cap energy to prevent numerical instability
    // If energy is extremely large, something is wrong
    if (energy > 1e10) {
        std::cerr << "Warning: Very large angle energy detected: " << energy 
                  << " (theta=" << theta << ", theta0=" << m_Theta0 
                  << ", delta=" << delta_theta << ")" << std::endl;
        // Cap at reasonable value to prevent explosion
        energy = 1e10;
    }
    
    // Check for NaN or Inf
    if (!std::isfinite(energy)) {
        std::cerr << "Warning: Non-finite angle energy detected: " << energy << std::endl;
        return 0.0;
    }
    
    return energy;
}

