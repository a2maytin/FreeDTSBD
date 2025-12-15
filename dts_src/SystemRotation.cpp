#include "SystemRotation.h"
#include <cmath>
#include <ctime>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SystemRotation::SystemRotation(int rotation_period) 
    : m_RotationPeriod(rotation_period),
      m_LastRotationStep(-1),
      m_RotationMatrix('I'),  // Start with identity
      m_InverseRotationMatrix('I'),  // Start with identity
      m_BoxCenter(0.0, 0.0, 0.0),
      m_Initialized(false),
      m_HasRotation(false),
      m_RandomGenerator(std::time(nullptr)),
      m_UniformDist(0.0, 1.0)
{
}

SystemRotation::~SystemRotation() {
}

void SystemRotation::Initialize(const Vec3D& box_center) {
    m_BoxCenter = box_center;
    m_Initialized = true;
}

bool SystemRotation::UpdateRotation(int step) {
    if (!IsEnabled() || !m_Initialized) {
        return false;
    }
    
    // Check if it's time to rotate
    if (m_RotationPeriod > 0 && step > 0 && step % m_RotationPeriod == 0 && step != m_LastRotationStep) {
        // Generate new rotation matrix
        m_NewRotationMatrix = GenerateRandomRotation();
        
        // Store the new rotation matrix (will be applied to vertices, then cumulative will be updated)
        // The cumulative update happens after vertices are rotated, via UpdateCumulativeRotation()
        
        // For now, we'll update cumulative here, but the actual rotation applied to vertices should be m_NewRotationMatrix
        // TODO: Move cumulative update to after vertices are rotated
        Tensor2 R = m_NewRotationMatrix;
        
        // Multiply with existing rotation (cumulative)
        // IMPORTANT: We want to apply the new rotation AFTER the existing rotations
        // So if we have R_old already applied, and we want to apply R_new, we do:
        // R_total = R_new * R_old (this means: first apply R_old, then R_new)
        // When we apply R_total to a vector: R_total * v = R_new * (R_old * v)
        // The inverse is: R_total^T = (R_new * R_old)^T = R_old^T * R_new^T
        // When we apply R_total^T: R_total^T * v = R_old^T * (R_new^T * v) (undoes R_new first, then R_old)
        m_RotationMatrix = R * m_RotationMatrix;
        
        // Re-orthonormalize to fix numerical drift
        // Rotation matrices must have orthonormal columns (or rows) and determinant = 1
        // We'll use Gram-Schmidt on the columns, then ensure right-handedness
        Vec3D col0(m_RotationMatrix(0, 0), m_RotationMatrix(1, 0), m_RotationMatrix(2, 0));
        Vec3D col1(m_RotationMatrix(0, 1), m_RotationMatrix(1, 1), m_RotationMatrix(2, 1));
        Vec3D col2(m_RotationMatrix(0, 2), m_RotationMatrix(1, 2), m_RotationMatrix(2, 2));
        
        // Normalize first column
        double norm0 = sqrt(col0(0)*col0(0) + col0(1)*col0(1) + col0(2)*col0(2));
        if (norm0 > 1e-10) {
            col0 = Vec3D(col0(0)/norm0, col0(1)/norm0, col0(2)/norm0);
        }
        
        // Make second column orthogonal to first, then normalize
        double dot01 = col0(0)*col1(0) + col0(1)*col1(1) + col0(2)*col1(2);
        col1 = Vec3D(col1(0) - dot01*col0(0), col1(1) - dot01*col0(1), col1(2) - dot01*col0(2));
        double norm1 = sqrt(col1(0)*col1(0) + col1(1)*col1(1) + col1(2)*col1(2));
        if (norm1 > 1e-10) {
            col1 = Vec3D(col1(0)/norm1, col1(1)/norm1, col1(2)/norm1);
        }
        
        // Third column = cross product of first two (ensures right-handedness and det = 1)
        col2 = Vec3D(
            col0(1)*col1(2) - col0(2)*col1(1),
            col0(2)*col1(0) - col0(0)*col1(2),
            col0(0)*col1(1) - col0(1)*col1(0)
        );
        
        // Reconstruct matrix from orthonormalized columns
        m_RotationMatrix.put(0, 0, col0(0)); m_RotationMatrix.put(0, 1, col1(0)); m_RotationMatrix.put(0, 2, col2(0));
        m_RotationMatrix.put(1, 0, col0(1)); m_RotationMatrix.put(1, 1, col1(1)); m_RotationMatrix.put(1, 2, col2(1));
        m_RotationMatrix.put(2, 0, col0(2)); m_RotationMatrix.put(2, 1, col1(2)); m_RotationMatrix.put(2, 2, col2(2));
        
        // Calculate inverse (transpose for rotation matrices)
        // For rotation matrices, R^T = R^-1, so this correctly inverts the cumulative rotation
        // (R_new * R_old)^T = R_old^T * R_new^T
        m_InverseRotationMatrix = m_RotationMatrix.Transpose();
        
        // Verify determinant is 1 (for rotation matrices) - will be printed in debug output below
        
        // Verify rotation matrix is orthogonal (R * R^T should be identity)
        // This is a sanity check - rotation matrices should be orthogonal
        Tensor2 should_be_identity = m_RotationMatrix * m_InverseRotationMatrix;
        double max_deviation = 0.0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double expected = (i == j) ? 1.0 : 0.0;
                double deviation = fabs(should_be_identity(i, j) - expected);
                if (deviation > max_deviation) max_deviation = deviation;
            }
        }
        if (max_deviation > 1e-10) {
            std::cerr << "ERROR: Rotation matrix is not orthogonal! Max deviation: " << max_deviation << std::endl;
            std::cerr << "This indicates a serious problem with the rotation matrix calculation.\n";
        }
        
        // Also verify that R^T * R is identity (both directions)
        Tensor2 should_be_identity2 = m_InverseRotationMatrix * m_RotationMatrix;
        double max_deviation2 = 0.0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double expected = (i == j) ? 1.0 : 0.0;
                double deviation = fabs(should_be_identity2(i, j) - expected);
                if (deviation > max_deviation2) max_deviation2 = deviation;
            }
        }
        if (max_deviation2 > 1e-10) {
            std::cerr << "ERROR: Inverse rotation matrix is incorrect! Max deviation: " << max_deviation2 << std::endl;
        }
        
        m_HasRotation = true;
        m_LastRotationStep = step;
        
        // Debug: Print rotation matrix and determinant
        double det = m_RotationMatrix(0, 0) * (m_RotationMatrix(1, 1)*m_RotationMatrix(2, 2) - m_RotationMatrix(1, 2)*m_RotationMatrix(2, 1))
                   - m_RotationMatrix(0, 1) * (m_RotationMatrix(1, 0)*m_RotationMatrix(2, 2) - m_RotationMatrix(1, 2)*m_RotationMatrix(2, 0))
                   + m_RotationMatrix(0, 2) * (m_RotationMatrix(1, 0)*m_RotationMatrix(2, 1) - m_RotationMatrix(1, 1)*m_RotationMatrix(2, 0));
        
        // Calculate angle of rotation around z-axis for debugging
        // For a rotation around z-axis, the angle is atan2(R(1,0), R(0,0))
        double angle_deg = atan2(m_RotationMatrix(1, 0), m_RotationMatrix(0, 0)) * 180.0 / M_PI;
        if (angle_deg < 0) angle_deg += 360.0;
        
        std::cout << "---> Rotation matrix at step " << step << " (determinant = " << det << ", z-axis angle = " << angle_deg << "°):\n";
        for (int i = 0; i < 3; i++) {
            std::cout << "     [";
            for (int j = 0; j < 3; j++) {
                std::cout << " " << m_RotationMatrix(i, j);
            }
            std::cout << " ]\n";
        }
        std::cout << "---> Inverse rotation matrix:\n";
        for (int i = 0; i < 3; i++) {
            std::cout << "     [";
            for (int j = 0; j < 3; j++) {
                std::cout << " " << m_InverseRotationMatrix(i, j);
            }
            std::cout << " ]\n";
        }
        
        return true;
    }
    
    return false;
}

Tensor2 SystemRotation::GenerateRandomRotation() {
    // Generate a random rotation matrix using the method of Arvo (1992)
    // This generates uniformly distributed rotations on SO(3)
    
    // Generate three random numbers
    double u1 = m_UniformDist(m_RandomGenerator);
    double u2 = m_UniformDist(m_RandomGenerator);
    double u3 = m_UniformDist(m_RandomGenerator);
    
    // Convert to rotation angles
    double theta = 2.0 * M_PI * u1;  // Rotation angle around axis
    double phi = acos(2.0 * u2 - 1.0);  // Polar angle for axis
    double psi = 2.0 * M_PI * u3;  // Rotation around z-axis
    
    // Build rotation axis from spherical coordinates
    Vec3D axis(sin(phi) * cos(psi), sin(phi) * sin(psi), cos(phi));
    
    // Rodrigues' rotation formula to build rotation matrix
    // R = I + sin(θ)[k]× + (1-cos(θ))[k]×²
    // where [k]× is the cross product matrix of the axis
    
    Tensor2 I('I');  // Identity matrix
    Tensor2 K_cross;  // Cross product matrix [k]×
    K_cross.put(0, 0, 0.0);
    K_cross.put(0, 1, -axis(2));
    K_cross.put(0, 2, axis(1));
    K_cross.put(1, 0, axis(2));
    K_cross.put(1, 1, 0.0);
    K_cross.put(1, 2, -axis(0));
    K_cross.put(2, 0, -axis(1));
    K_cross.put(2, 1, axis(0));
    K_cross.put(2, 2, 0.0);
    
    Tensor2 K_cross_sq = K_cross * K_cross;  // [k]×²
    
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    
    // Build rotation matrix: R = I + sin(θ)[k]× + (1-cos(θ))[k]×²
    Tensor2 R = I + (K_cross * sin_theta) + (K_cross_sq * (1.0 - cos_theta));
    
    return R;
}

void SystemRotation::ApplyRotationToVertex(Vec3D& pos) const {
    if (!IsEnabled()) {
        return;
    }
    // Apply the cumulative rotation matrix (for output inverse rotation)
    pos = m_RotationMatrix * pos;
}

void SystemRotation::ApplyNewRotationToVertex(Vec3D& pos) const {
    if (!IsEnabled()) {
        return;
    }
    // Apply only the new rotation matrix (not cumulative)
    pos = m_NewRotationMatrix * pos;
}

void SystemRotation::ApplyInverseRotationToVertex(Vec3D& pos) const {
    if (!IsEnabled()) {
        return;
    }
    pos = m_InverseRotationMatrix * pos;
}

Vec3D SystemRotation::GetRotatedPosition(const Vec3D& pos) const {
    if (!IsEnabled()) {
        return pos;
    }
    return m_RotationMatrix * pos;
}

Vec3D SystemRotation::GetInverseRotatedPosition(const Vec3D& pos) const {
    if (!IsEnabled()) {
        return pos;
    }
    return m_InverseRotationMatrix * pos;
}

