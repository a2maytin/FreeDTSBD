#ifndef SYSTEM_ROTATION_H
#define SYSTEM_ROTATION_H

#include "Tensor2.h"
#include "Vec3D.h"
#include <random>

/**
 * @brief SystemRotation: Manages periodic random rotations to randomize multithreading bias direction
 * 
 * This class addresses the issue where multithreading bias causes conformational changes
 * to prefer the 111 direction. By periodically rotating the entire system by random angles,
 * the bias direction is randomized relative to the system, while output coordinates remain
 * consistent by applying the inverse rotation.
 * 
 * Usage:
 * - Initialize with rotation period (how often to rotate, in steps)
 * - Call UpdateRotation(step) periodically to check if rotation is needed
 * - Call ApplyRotationToVertex() when reading vertex positions for simulation
 * - Call ApplyInverseRotationToVertex() when writing vertex positions to output
 */
class SystemRotation {
public:
    SystemRotation(int rotation_period = 0);
    ~SystemRotation();
    
    /**
     * @brief Initialize: store the box center
     * This should be called once at the beginning of the simulation
     * @param box_center The box center (where system is centered)
     */
    void Initialize(const Vec3D& box_center);
    
    /**
     * @brief Check if rotation should be applied at this step, and apply if needed
     * @param step Current simulation step
     * @return true if rotation was applied, false otherwise
     */
    bool UpdateRotation(int step);
    
    /**
     * @brief Get the box center
     */
    const Vec3D& GetBoxCenter() const { return m_BoxCenter; }
    
    /**
     * @brief Apply rotation to a vertex position (for simulation use)
     * @param pos Input position (will be modified)
     */
    void ApplyRotationToVertex(Vec3D& pos) const;
    
    /**
     * @brief Apply the new rotation matrix to a vertex position (use this when rotating the system)
     * This applies only the new rotation, not the cumulative rotation
     * @param pos Input position (will be modified)
     */
    void ApplyNewRotationToVertex(Vec3D& pos) const;
    
    /**
     * @brief Apply inverse rotation to a vertex position (for output)
     * @param pos Input position (will be modified)
     */
    void ApplyInverseRotationToVertex(Vec3D& pos) const;
    
    /**
     * @brief Get rotated position without modifying input
     */
    Vec3D GetRotatedPosition(const Vec3D& pos) const;
    
    /**
     * @brief Get inverse-rotated position without modifying input (for output)
     */
    Vec3D GetInverseRotatedPosition(const Vec3D& pos) const;
    
    /**
     * @brief Check if rotation is enabled
     */
    bool IsEnabled() const { return m_RotationPeriod > 0; }
    
    /**
     * @brief Check if any rotation has been applied
     */
    bool HasRotation() const { return m_HasRotation; }
    
    /**
     * @brief Get current rotation matrix (for debugging)
     */
    const Tensor2& GetRotationMatrix() const { return m_RotationMatrix; }
    
    /**
     * @brief Get inverse rotation matrix (for debugging)
     */
    const Tensor2& GetInverseRotationMatrix() const { return m_InverseRotationMatrix; }
    
    /**
     * @brief Get the new rotation matrix that was just generated (before cumulative update)
     * This should be called right after UpdateRotation returns true, to apply the new rotation
     * @return The new rotation matrix to apply to vertices
     */
    const Tensor2& GetNewRotationMatrix() const { return m_NewRotationMatrix; }
    
    /**
     * @brief Restore rotation state from restart file
     * Used when restarting a simulation to restore the rotation matrices and state
     * @param rotation_matrix The cumulative rotation matrix to restore
     * @param inverse_rotation_matrix The inverse rotation matrix to restore
     * @param last_rotation_step The last step when rotation was applied
     * @param box_center The box center
     */
    void RestoreRotationState(const Tensor2& rotation_matrix, 
                              const Tensor2& inverse_rotation_matrix,
                              int last_rotation_step,
                              const Vec3D& box_center);
    
    /**
     * @brief Get the last rotation step (for saving to restart file)
     */
    int GetLastRotationStep() const { return m_LastRotationStep; }
    
    /**
     * @brief Check if rotation system is initialized
     */
    bool IsInitialized() const { return m_Initialized; }
    
private:
    /**
     * @brief Generate a random rotation matrix (uniformly distributed on SO(3))
     * @return The newly generated rotation matrix
     */
    Tensor2 GenerateRandomRotation();
    
    int m_RotationPeriod;              // How often to rotate (0 = disabled)
    int m_LastRotationStep;            // Last step when rotation was applied
    Tensor2 m_RotationMatrix;          // Cumulative rotation matrix (total rotation applied)
    Tensor2 m_InverseRotationMatrix;   // Inverse of rotation matrix (for output)
    Tensor2 m_NewRotationMatrix;       // New rotation matrix (just generated, before cumulative update)
    Vec3D m_BoxCenter;                 // Box center (where system is centered)
    bool m_Initialized;                // Whether Initialize() has been called
    bool m_HasRotation;                // Whether any rotation has been applied
    std::mt19937 m_RandomGenerator;    // Random number generator
    std::uniform_real_distribution<double> m_UniformDist;  // Uniform distribution [0,1)
};

#endif // SYSTEM_ROTATION_H
