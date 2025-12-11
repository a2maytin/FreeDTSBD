#ifndef AFX_RepulsionBetweenVerticesAndDNABeads_H_334B21B8_INCLUDED_
#define AFX_RepulsionBetweenVerticesAndDNABeads_H_334B21B8_INCLUDED_

#include "vertex.h"
#include "AbstractNonbondedInteractionBetweenVertices.h"

/**
 * @brief Class implementing repulsion between membrane vertices and DNA beads.
 * 
 * This class enforces repulsion between membrane vertices (vertices with triangles/links)
 * and DNA beads (vertices with bonds). The repulsion follows a simple inverse distance
 * squared law within a cutoff radius.
 * 
 * Author: Adapted from existing repulsion implementations
 * Copyright: Based on FreeDTS codebase
 */
class RepulsionBetweenVerticesAndDNABeads : public AbstractNonbondedInteractionBetweenVertices {
public:
    RepulsionBetweenVerticesAndDNABeads(State* pstate, std::string input_data);
    ~RepulsionBetweenVerticesAndDNABeads();

    double GetVertexNonBondedEnergy(vertex* pvertex) const override;
    std::string CurrentState() const override;

    void Initialize() override;
    inline std::string GetDerivedDefaultReadName() const override { return "RepulsionBetweenVerticesAndDNABeads"; }
    inline static std::string GetDefaultReadName() { return "RepulsionBetweenVerticesAndDNABeads"; }

private:
    std::string m_Input_Data;
    State *m_pState;
    
    std::vector<vertex*> m_pMembraneVertices;  // Vertices that are part of the membrane
    std::vector<vertex*> m_pDNABeads;           // Vertices that are DNA beads (have bonds)
    double m_EP;                                // Repulsion strength parameter
    double m_R0;                                // Cutoff radius
    double m_R0_2;                              // Cutoff radius squared (for efficiency)
    Vec3D *m_pBox;
    
    /**
     * @brief Check if a vertex is a membrane vertex (has triangles/links).
     * @param v Pointer to vertex to check.
     * @return True if vertex is part of membrane mesh.
     */
    bool IsMembraneVertex(vertex* v) const;
    
    /**
     * @brief Check if a vertex is a DNA bead (has bonds).
     * @param v Pointer to vertex to check.
     * @return True if vertex is a DNA bead.
     */
    bool IsDNABead(vertex* v) const;
    
    /**
     * @brief Find nearby DNA beads for a given membrane vertex using voxelization.
     * @param pV Pointer to membrane vertex.
     * @return Vector of nearby DNA bead vertices.
     */
    std::vector<vertex*> FindNearbyDNABeads(vertex* pV) const;
    
    /**
     * @brief Find nearby membrane vertices for a given DNA bead using voxelization.
     * @param pV Pointer to DNA bead vertex.
     * @return Vector of nearby membrane vertices.
     */
    std::vector<vertex*> FindNearbyMembraneVertices(vertex* pV) const;
    
    /**
     * @brief Calculate repulsion energy between a membrane vertex and a DNA bead.
     * @param v_membrane Pointer to membrane vertex.
     * @param v_dna Pointer to DNA bead vertex.
     * @return Repulsion energy.
     */
    double CalculateRepulsionEnergy(vertex* v_membrane, vertex* v_dna) const;

};

#endif

