# Implementation Plan: Linearly Connected Untriangulated Vertices

## Overview
Re-implement support for linearly connected untriangulated vertices (previously called "DNA beads"). These are vertices that:
- Are connected via harmonic bonds in linear chains
- Have no triangulation (no links/triangles)
- Should be excluded from membrane-specific calculations (curvature, area, volume)
- Should only have bond energy and optional nonbonded interactions

## Naming Convention
- **Terminology**: Use "bonded vertex" or "untriangulated vertex" instead of "DNA"
- **Variable names**: `is_bonded_vertex`, `is_untriangulated_vertex`, or check `GetBonds().empty()`
- **Comments**: Use "linearly connected untriangulated vertex" or "bonded vertex"

## Files to Modify

### 1. Core Vertex Evolution (`EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.cpp`)
**Changes:**
- Add check: `bool is_bonded_vertex = !pvertex->GetBonds().empty();`
- **Skip membrane energy calculations** for bonded vertices:
  - Skip `old_energy` accumulation from neighbors/triangles/links
  - Skip `new_energy` calculation via `SingleVertexEnergy()`
  - Skip curvature updates
  - Skip triangle/link shape operator updates
- **Add bond energy** for bonded vertices:
  - Before move: `bond_energy = -(pvertex->GetBondEnergyOfVertex());`
  - After move: `bond_energy += pvertex->GetBondEnergyOfVertex();`
  - Add to `tot_diff_energy` only for bonded vertices
- **Skip global mesh properties** (volume/area/curvature) for bonded vertices
- **Skip force calculations** (from inclusions, vector fields, etc.) for bonded vertices
- **Skip face angle checks** for bonded vertices (they have no faces)

**Key Pattern:**
```cpp
bool is_bonded_vertex = !pvertex->GetBonds().empty();

if (!is_bonded_vertex) {
    // Standard membrane vertex processing (identical to v2.1)
    // ... all existing membrane energy calculations ...
} else {
    // Bonded vertex: only bond energy, no membrane energy
    old_energy = 0;
    new_energy = 0;
    // ... bond energy calculations ...
}
```

### 2. Global Mesh Properties (`VAHGlobalMeshProperties.cpp`)
**Changes:**
- In `CalculateAVertexRingContributionToGlobalVariables()`:
  - Skip volume/area/curvature calculations if vertex has bonds
  - Check: `if (!p_vertex->GetBonds().empty()) return;` or skip the calculations

**Rationale:** Bonded vertices don't contribute to membrane volume/area/curvature.

### 3. Curvature Calculation (`CurvatureByShapeOperatorType1.cpp`)
**Changes:**
- In `UpdateVertexCurvature()`:
  - Early return if vertex has no links: `if (pvertex->GetVLinkList().empty()) return true;`
  - This naturally excludes bonded vertices (they have no links)

**Rationale:** Bonded vertices have no triangulation, so no curvature.

### 4. Mesh Quality Checks (`MC_Simulation.cpp`)
**Changes:**
- In `CheckMesh()`:
  - Skip distance checks for bonded vertices
  - Check: `if (!p_v1->GetBonds().empty() || !p_v2->GetBonds().empty()) continue;`

**Rationale:** Bonded vertices don't need to satisfy membrane edge length constraints.

### 5. Mesh Distance Checks (`MESH.cpp`)
**Changes:**
- In distance checking functions (e.g., `CheckMesh()`):
  - Skip checks between bonded vertices
  - Skip checks if either vertex is bonded

**Rationale:** Bonded vertices can be closer together than membrane vertices.

### 6. Non-OpenMP Vertex Evolution (`EvolveVerticesByMetropolisAlgorithm.cpp`)
**Changes:**
- Same pattern as `EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.cpp`
- Add bonded vertex checks and skip membrane calculations
- Add bond energy calculations

### 7. Alexander Moves (`AlexanderMoveByMetropolisAlgorithmWithOpenMP.cpp`)
**Changes:**
- **No changes needed** - Alexander moves only affect membrane edges
- Bonded vertices have no links, so they're automatically excluded

### 8. State Initialization (`State.cpp`)
**Changes:**
- **Already implemented** - `HarmonicBondsList` initialization is present
- Verify that `Initialize()` is called on `HarmonicBondsList` after mesh is loaded

## Implementation Order

### Phase 1: Core Infrastructure (Already Done)
- ✅ `bond` class exists
- ✅ `VertexHarmonicBounds` class exists  
- ✅ `HarmonicBondsList` class exists
- ✅ `State.cpp` reads `HarmonicBondsList` from input

### Phase 2: Vertex Evolution (Priority 1)
1. Modify `EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.cpp`
   - Add bonded vertex detection
   - Skip membrane energy for bonded vertices
   - Add bond energy calculations
   - Test with simple case

2. Modify `EvolveVerticesByMetropolisAlgorithm.cpp` (non-OpenMP version)
   - Same changes as above

### Phase 3: Global Properties (Priority 2)
3. Modify `VAHGlobalMeshProperties.cpp`
   - Skip volume/area/curvature for bonded vertices

4. Modify `CurvatureByShapeOperatorType1.cpp`
   - Skip curvature for vertices with no links

### Phase 4: Mesh Checks (Priority 3)
5. Modify `MC_Simulation.cpp`
   - Skip mesh quality checks for bonded vertices

6. Modify `MESH.cpp`
   - Skip distance checks for bonded vertices

### Phase 5: Testing & Validation
7. Test with simple linear chain
8. Verify energy conservation
9. Verify acceptance rates match expected values
10. Test with mixed system (membrane + bonded vertices)

## Key Principles

1. **Preserve v2.1 behavior**: When `GetBonds().empty() == true`, code should be identical to v2.1
2. **Minimal changes**: Only add conditional checks, don't refactor existing code
3. **Clear separation**: Bonded vertices have completely different energy landscape
4. **No membrane properties**: Bonded vertices don't contribute to membrane volume/area/curvature
5. **Bond energy only**: For bonded vertices, only calculate bond energy (and optional nonbonded)

## Testing Strategy

1. **Unit test**: Single bonded vertex pair, verify bond energy calculation
2. **Integration test**: Linear chain of 10 bonded vertices
3. **Mixed test**: Membrane with embedded linear chain
4. **Energy test**: Verify total energy is conserved
5. **Acceptance rate test**: Verify acceptance rates are reasonable

## Potential Issues to Watch For

1. **Voxelization**: Bonded vertices should still be voxelized for neighbor finding
2. **Nonbonded interactions**: Optional feature for repulsion between bonded vertices and membrane
3. **Initialization order**: Bonds must be initialized after mesh is loaded
4. **File I/O**: Bond file format must be correct (vertex_id1 vertex_id2 k l0)

## Notes

- The `GetBondEnergyOfVertex()` method already exists in `VertexHarmonicBounds`
- The infrastructure is mostly in place, we just need to add the conditional logic
- This is a clean re-implementation, avoiding the problematic "safety checks" that caused issues before
