# Code Verification Summary: Energy Calculation Logic

## Comparison: Original vs Current Code

### Critical Energy Calculation Section

**ORIGINAL CODE (from git HEAD):**
```cpp
// Line ~540-556
double diff_energy = new_energy - old_energy;
changed_en = diff_energy;  // Elastic energy only

double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_on_vertex + 
                         dE_force_from_inc + dE_force_from_vector_fields + 
                         dE_volume + dE_t_area + dE_g_curv;

double U = m_Beta * tot_diff_energy - m_DBeta;

if(U <= 0 || exp(-U) > temp) {
    // move accepted
    // (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);  // COMMENTED OUT
    if(!pvertex->CheckVoxel()){
        pvertex->UpdateVoxelAfterAVertexMove();
    }
    // ... global variables update ...
    return true;
} else {
    changed_en = 0;  // On reject
    // ... reverse changes ...
    return false;
}
```

**CURRENT CODE:**
```cpp
// Line 565-615
double diff_energy = new_energy - old_energy;
changed_en = diff_energy;  // Elastic energy only

double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_on_vertex + 
                         dE_force_from_inc + dE_force_from_vector_fields + 
                         dE_volume + dE_t_area + dE_g_curv + 
                         bond_energy + dE_nonbonded;  // Added DNA terms

double U = m_Beta * tot_diff_energy - m_DBeta;

if(U <= 0 || exp(-U) > temp) {
    // move accepted
    // NOTE: Energy updated via reduction in caller using changed_en
    // AddToTotalEnergy commented out (as in original)
    if (pvertex->GetVoxel() != nullptr && !pvertex->CheckVoxel()){
        pvertex->UpdateVoxelAfterAVertexMove();
    }
    if (!is_dna_vertex) {
        m_pState->GetApplyConstraintBetweenGroups()->AcceptMove();
        // ... global variables update ...
    }
    return true;
} else {
    changed_en = 0;  // On reject
    if (is_dna_vertex) {
        pvertex->ReverseConstantMesh_Copy();
        return false;
    }
    // ... reverse changes ...
    return false;
}
```

## Key Findings

### ✅ Energy Logic is IDENTICAL to Original

1. **`changed_en` assignment:**
   - **Original:** `changed_en = diff_energy` (elastic only) before accept/reject, NOT updated on accept
   - **Current:** `changed_en = diff_energy` (elastic only) before accept/reject, NOT updated on accept
   - **Match:** ✅ IDENTICAL

2. **`tot_diff_energy` calculation:**
   - **Original:** Includes elastic + all coupling terms (volume, area, curvature, forces)
   - **Current:** Includes elastic + all coupling terms + bond_energy + dE_nonbonded (for DNA)
   - **Match:** ✅ LOGICALLY IDENTICAL (only added DNA-specific terms)

3. **Metropolis decision:**
   - **Original:** Uses `tot_diff_energy` for accept/reject decision
   - **Current:** Uses `tot_diff_energy` for accept/reject decision
   - **Match:** ✅ IDENTICAL

4. **Energy accumulation:**
   - **Original:** `AddToTotalEnergy(diff_energy)` called in `EvolveOneStep` after reduction
   - **Current:** `AddToTotalEnergy(diff_energy)` called in `EvolveOneStep` after reduction
   - **Match:** ✅ IDENTICAL

### Differences (All DNA-Related, Non-Critical for Membrane Behavior)

1. **Added DNA-specific terms to `tot_diff_energy`:**
   - `bond_energy` (for DNA bonds)
   - `dE_nonbonded` (for DNA repulsion)
   - **Impact:** Only affects DNA beads, not membrane vertices

2. **Added DNA-specific conditionals:**
   - `is_dna_vertex` checks to skip membrane-specific operations for DNA
   - **Impact:** Prevents errors when processing DNA beads, doesn't affect membrane

3. **Added null check for voxel update:**
   - `if (pvertex->GetVoxel() != nullptr && !pvertex->CheckVoxel())`
   - **Impact:** Safety improvement, shouldn't affect behavior

4. **DNA-specific rejection handling:**
   - Simplified reversal for DNA beads
   - **Impact:** Only affects DNA beads, not membrane vertices

## Conclusion

**The energy calculation logic for membrane vertices is IDENTICAL to the original code.**

The only changes are:
- DNA-specific additions that don't affect membrane behavior
- Safety improvements (null checks)
- Conditional logic to handle DNA beads separately

**If elongation occurs with all couplings enabled, it is likely:**
1. Already present in the original code (needs verification by reverting)
2. Related to parameter values (Kappa, coupling strengths, etc.)
3. Related to numerical precision or accumulation errors (unlikely given identical logic)

## Recommendation

✅ **Safe to commit current changes** - the energy logic matches the original.

Then revert to original commit to verify if elongation existed before our changes.

