#ifndef AFX_DNARepulsion_H_334B21B8_INCLUDED_
#define AFX_DNARepulsion_H_334B21B8_INCLUDED_

#include <unordered_set>
#include <unordered_map>
#include "vertex.h"
#include "AbstractNonbondedInteractionBetweenVertices.h"

// Forward declaration
class State;

/**
 * @brief Class implementing repulsion involving DNA beads.
 * 
 * This class enforces repulsion between:
 * - Membrane vertices (vertices with triangles/links) and DNA beads (vertices with bonds)
 * - DNA beads and other DNA beads (excluding bonded neighbors)
 * 
 * The repulsion follows a soft cosine potential (LAMMPS-style) that smoothly goes to zero at the cutoff radius.
 * 
 * Author: Adapted from existing repulsion implementations
 * Copyright: Based on FreeDTS codebase
 */
class DNARepulsion : public AbstractNonbondedInteractionBetweenVertices {
public:
    DNARepulsion(State* pstate, std::string input_data);
    ~DNARepulsion();

    double GetVertexNonBondedEnergy(vertex* pvertex) const override;
    std::string CurrentState() const override;

    void Initialize() override;
    inline std::string GetDerivedDefaultReadName() const override { return "DNARepulsion"; }
    inline static std::string GetDefaultReadName() { return "DNARepulsion"; }

private:
    std::string m_Input_Data;
    State *m_pState;
    
    std::vector<vertex*> m_pMembraneVertices;  // Vertices that are part of the membrane
    std::vector<vertex*> m_pDNABeads;           // Vertices that are DNA beads (have bonds)
    
    // Performance optimization: Hash sets for O(1) lookup instead of O(N) linear search
    std::unordered_set<vertex*> m_MembraneVertexSet;  // Fast lookup for membrane vertices
    std::unordered_set<vertex*> m_DNABeadSet;          // Fast lookup for DNA beads
    
    // Performance optimization: Cache neighbor exclusion sets per DNA bead
    // Maps vertex pointer to set of excluded neighbors (bonded + next-nearest)
    mutable std::unordered_map<vertex*, std::unordered_set<vertex*> > m_NeighborExclusionCache;
    
    double m_EP;                                // Repulsion strength parameter
    double m_R0;                                // Cutoff radius for membrane-DNA interactions
    double m_R0_2;                              // Cutoff radius squared for membrane-DNA (for efficiency)
    double m_R0_DNA;                            // Cutoff radius for DNA-DNA interactions
    double m_R0_DNA_2;                          // Cutoff radius squared for DNA-DNA (for efficiency)
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
     * @brief Find nearby DNA beads for a given DNA bead, excluding bonded neighbors.
     * @param pV Pointer to DNA bead vertex.
     * @return Vector of nearby DNA bead vertices (excluding bonded neighbors).
     */
    std::vector<vertex*> FindNearbyDNABeadsForDNA(vertex* pV) const;
    
    /**
     * @brief Calculate repulsion energy between two vertices (membrane-DNA or DNA-DNA).
     * @param v1 Pointer to first vertex.
     * @param v2 Pointer to second vertex.
     * @param cutoff_radius_squared Cutoff radius squared to use for this interaction.
     * @return Repulsion energy.
     */
    double CalculateRepulsionEnergy(vertex* v1, vertex* v2, double cutoff_radius_squared) const;
    
    /**
     * @brief Get cached neighbor exclusion set for a DNA bead (bonded + next-nearest neighbors).
     * @param pV Pointer to DNA bead vertex.
     * @return Set of excluded neighbor vertices.
     */
    const std::unordered_set<vertex*>& GetNeighborExclusionSet(vertex* pV) const;
    
    /**
     * @brief Validate vertex coordinates (helper function to reduce redundant checks).
     * @param v Pointer to vertex.
     * @return True if coordinates are valid.
     */
    bool ValidateCoordinates(vertex* v) const;

};

#endif

