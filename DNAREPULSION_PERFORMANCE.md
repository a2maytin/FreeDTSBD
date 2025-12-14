# DNARepulsion Computational Cost Analysis

## How FreeDTS Handles Nonbonded Interactions

### 1. **Voxelization-Based Neighbor Search** (Primary Method)

FreeDTS uses **spatial voxelization** to efficiently find nearby vertices:

- **Voxel Grid**: The simulation box is divided into a 3D grid of voxels
- **Typical Voxel Size**: ~1.05-1.2 nm (configurable via `Voxel_Size` in input file)
- **Neighbor Search**: For each vertex, checks only **26 neighboring voxel cells** (3×3×3 grid minus center cell)
- **Complexity**: 
  - Finding voxel: O(1)
  - Checking 26 cells: O(1) 
  - Processing vertices in those cells: O(N_voxel) where N_voxel = average vertices per 26-cell region

### 2. **Fallback Method** (If Voxel Missing)

If a vertex doesn't have a voxel assigned:
- **Brute Force**: Checks ALL DNA beads or membrane vertices
- **Complexity**: O(N) where N = total number of DNA beads or membrane vertices
- **Note**: This is much slower and should be avoided

### 3. **Current Implementation Details**

#### For DNA-DNA Interactions:
1. **Neighbor Exclusion**: Uses cached sets to exclude bonded + next-nearest neighbors (O(1) lookup)
2. **Voxel Search**: Checks 26 neighboring voxel cells
3. **Distance Filtering**: Calculates distance² for each candidate (cheap: ~10-20 FLOPs)
4. **Energy Calculation**: Only for vertices within cutoff (power law: ~20-30 FLOPs)

#### For DNA-Membrane Interactions:
1. **Voxel Search**: Checks 26 neighboring voxel cells
2. **Distance Filtering**: Calculates distance² for each candidate
3. **Energy Calculation**: Only for vertices within cutoff

### 4. **Performance Characteristics**

#### Cost Per Vertex Move:
- **Neighbor Finding**: ~26 voxel cells × (average vertices per cell)
  - Typical: 26 × 5-20 vertices = 130-520 distance calculations
- **Distance Calculations**: ~10-20 FLOPs each (squared distance, no sqrt needed for cutoff check)
- **Energy Calculations**: Only for vertices within cutoff
  - Typical: 5-50 neighbors within cutoff
  - Cost: ~20-30 FLOPs per energy calculation (power law with sqrt)

#### Total Cost Per Move:
- **Best Case** (sparse system): ~100-200 distance checks, ~10-20 energy calculations = **~2000-4000 FLOPs**
- **Worst Case** (dense system): ~500-1000 distance checks, ~50-100 energy calculations = **~15000-30000 FLOPs**

### 5. **Effect of Increasing Cutoff Radius**

**YES, increasing the cutoff radius WILL slow things down**, but the scaling is **linear**, not quadratic:

#### Why Linear Scaling?
- **Voxel search is fixed**: Still only checks 26 neighboring voxel cells (doesn't change)
- **More neighbors pass distance check**: As cutoff increases, more vertices in those 26 cells are within cutoff
- **More energy calculations**: Each additional neighbor within cutoff requires one energy calculation

#### Scaling Formula:
```
Cost ∝ (cutoff_radius / voxel_size)³ × density
```

For example:
- **Cutoff = 1.0 nm, Voxel = 1.2 nm**: ~10-20 neighbors within cutoff
- **Cutoff = 2.0 nm, Voxel = 1.2 nm**: ~40-80 neighbors within cutoff (4× more)
- **Cutoff = 3.0 nm, Voxel = 1.2 nm**: ~90-180 neighbors within cutoff (9× more)

**Important Limitation**: The current implementation only checks 26 immediate neighboring voxels. If your cutoff radius is **larger than ~1.5× voxel_size**, you may miss some interactions! The code does have a distance check, so it correctly filters, but it won't find neighbors beyond those 26 cells.

### 6. **Optimizations in Place**

1. **Hash Sets**: O(1) lookup for vertex type (membrane vs DNA)
2. **Cached Neighbor Exclusion**: Bonded/next-nearest neighbors cached per vertex
3. **Early Distance Filtering**: Uses distance² (no sqrt) for cutoff check
4. **Voxel-Based Search**: Reduces search space from O(N) to O(N_voxel)
5. **OpenMP Thread Safety**: Critical sections protect neighbor cache

### 7. **Recommendations**

1. **Keep cutoff radius ≤ 1.5× voxel_size** for optimal performance
   - If you need larger cutoff, increase voxel size in input file: `Voxel_Size lx ly lz`
   
2. **For DNA-DNA interactions**: 
   - Typical cutoff: 1.0-1.5 nm works well
   - Larger cutoffs (2.0+ nm) will significantly slow down simulation
   
3. **For DNA-membrane interactions**:
   - Can use larger cutoff (2.0-3.0 nm) since membrane vertices are typically sparser
   
4. **Power Law Exponent**:
   - Higher power (12) = steeper repulsion = fewer close contacts = potentially faster
   - Lower power (3-6) = gentler repulsion = more close contacts = potentially slower

### 8. **Typical Performance**

For a system with:
- 2000 DNA beads
- 1000 membrane vertices
- Cutoff = 1.5 nm
- Voxel size = 1.2 nm

**Expected cost per vertex move**: ~5000-10000 FLOPs for nonbonded interactions
**Relative to other energies**: Nonbonded interactions are typically 10-30% of total energy calculation time

### 9. **Potential Improvements** (Future)

1. **Multi-level voxelization**: Check more voxel layers for larger cutoffs
2. **Neighbor lists**: Update less frequently than every move
3. **Verlet lists**: Only rebuild when particles move significant distance
4. **Cell lists**: Alternative to voxelization for very large systems
