# Spatial Drift Analysis - Potential Causes

## Summary
The drift occurs in Phase 4.12 (random selection) but NOT in Phase 4.10 (sequential selection), indicating the issue is related to how random selection interacts with parallel execution.

## Confirmed Issues

### 1. **Non-Thread-Safe Shared RNG** ⚠️ CRITICAL
- **Location**: `State.cpp:52` - `m_RandomNumberGenerator(new RNG(1234))`
- **Problem**: Single shared instance of `RNG` used by ALL threads
- **RNG Implementation**: Uses `std::default_random_engine` which is NOT thread-safe
- **Impact**: 
  - Multiple threads calling `IntRNG()` concurrently cause race conditions
  - RNG state corruption leads to correlated/identical random numbers
  - Multiple threads may select the SAME vertex simultaneously
  - Creates lock contention hotspots
- **Evidence**: 
  - 30 calls to `GetRandomNumberGenerator()` in `EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.cpp`
  - Both surface and edge vertex loops use the same shared RNG
  - Phase 4.12 shows only ~47% of selected vertices are processed (many skipped due to locking)

### 2. **Lock Contention Hotspots**
- **Location**: `vertex.cpp:532-575` - `CheckLockVertex()` and `CheckLockNeighbourVertex()`
- **Problem**: 
  - With random selection, multiple threads compete for the same vertices
  - Lock contention creates processing bias (some vertices processed more/less often)
  - Position-dependent bias if vertices are spatially ordered in array
- **Evidence**:
  - Phase 4.10 (sequential): 100% of selected vertices processed (std_dev=0.0)
  - Phase 4.12 (random): Only ~47% processed, std_dev=308.2 (high variance)

### 3. **Vertex Array Ordering** ✅ CONFIRMED SPATIAL CORRELATION
- **Location**: `MESH.cpp:130-141` - Vertices added to `m_pActiveV` in order from `m_Vertex`
- **Problem**: 
  - DNA vertices loaded sequentially from `.q` file in chain order (ID 0, 1, 2, ... 1999)
  - Vertices stored in array in same order as loaded
  - **Array index = Vertex ID ≈ Spatial position along chain**
  - Random selection uses array index, creating spatial correlation
- **Evidence**:
  - `.q` file shows sequential IDs (0-1999) with spatially connected coordinates
  - First vertices: ~(106, 109, 112) - middle of box
  - Last vertices: ~(104-106, 108-109, 112) - similar region
  - Vertices form a chain, so array index correlates with chain position
  - **This creates spatial correlation between array index and physical position**

## Potential Issues (Need Verification)

### 4. **OpenMP Schedule Dynamic**
- **Location**: `#pragma omp for schedule(dynamic)`
- **Potential Problem**: 
  - Dynamic scheduling may cluster threads around certain array regions
  - If certain vertices take longer (due to locking), threads may cluster
  - Could create spatial correlations if vertices are spatially ordered
- **Status**: Unlikely to be primary cause, but could amplify other issues

### 5. **Voxelization (Membrane Only)**
- **Location**: `EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.cpp:1590-1593`
- **Status**: 
  - Only affects membrane vertices (skipped for DNA)
  - Doesn't lock vertices, just spatial partitioning
  - **NOT a cause for DNA-only simulations**

### 6. **Other Parallel Operations**
- **Location**: `MC_Simulation.cpp:102-108`
- **Status**: 
  - Alexander moves, inclusion updates run sequentially (not simultaneously)
  - Don't interfere with vertex locking
  - **NOT a cause**

## Root Cause Hypothesis

**Primary Cause**: Non-thread-safe shared RNG + Random selection = Lock contention hotspots

**Mechanism**:
1. Multiple threads call `IntRNG()` concurrently on shared RNG
2. Race conditions cause RNG state corruption
3. Correlated/identical random numbers generated
4. Multiple threads select same vertex simultaneously
5. Lock contention creates processing bias
6. **Vertices are spatially ordered in array (array index ≈ chain position)**
7. **Lock contention hotspots correlate with array indices**
8. **Array index correlation → Spatial position correlation → Position-dependent drift**

**Why Sequential Selection Works**:
- Sequential selection (`vid = i % no_surf_v`) doesn't use RNG for selection
- Threads process different vertices (no contention)
- 100% of selected vertices are processed
- No position-dependent bias

## Recommended Fixes

### Fix 1: Thread-Local RNG (RECOMMENDED)
- Create thread-local RNG instances
- Each thread has its own RNG with unique seed
- Eliminates race conditions and correlations

### Fix 2: Thread-Safe RNG Wrapper
- Add mutex protection around RNG calls
- Slower but safer
- Maintains single RNG instance

### Fix 3: Use Sequential Selection
- Keep Phase 4.10 approach (sequential selection)
- Eliminates lock contention
- But loses "true randomness" of selection

## Testing Recommendations

1. **✅ Vertex Array Ordering VERIFIED**:
   - DNA vertices ARE spatially ordered (chain order in array)
   - Array index correlates with chain position
   - This confirms spatial correlation hypothesis

2. **Test Thread-Local RNG**:
   - Implement thread-local RNG
   - Compare Phase 4.12 results with thread-local vs shared RNG
   - Should eliminate drift if RNG is the cause

3. **Track RNG State**:
   - Add logging to track RNG calls per thread
   - Check for identical random numbers across threads
   - Verify race conditions

4. **Lock Contention Analysis**:
   - Track which vertices are locked most often
   - Check for spatial clustering of lock contention
   - Verify position-dependent processing bias
