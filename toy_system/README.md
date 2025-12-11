# Toy System for Testing DNA Bead Repulsion

This is a minimal test system to verify the DNA bead repulsion functionality.

## Files

1. **membrane.q** - Topology file containing:
   - 9 membrane vertices forming a 3x3 grid
   - 3 DNA beads in a chain
   - 4 triangles connecting membrane vertices

2. **dna_bonds.txt** - Bond file connecting DNA beads:
   - Bead 9 connected to bead 10 (equilibrium distance: 2.0)
   - Bead 10 connected to bead 11 (equilibrium distance: 2.0)

3. **topology.top** - Topology reference file

4. **input.dts** - Simulation input file with:
   - DNA bonds enabled
   - Repulsion between vertices and DNA beads enabled
   - Basic simulation parameters

## Compilation (First Time Only)

If you haven't compiled DTS yet, you need to compile it first (this includes the new repulsion code):

```bash
# From the project root directory
cd /home/andrew/Desktop/Projects/FreeDTSBD
bash compile.sh
```

This will create the DTS executable in the project root.

## Running the Simulation

```bash
# From the toy_system directory
cd toy_system
./run_test.sh

# Or manually:
../DTS -in input.dts -top topology.top
```

**Note**: The command uses `-in` not `-f` for the input file flag.

## What to Expect

- DNA beads should maintain their bond distances (~2.0)
- DNA beads should be repelled from membrane vertices
- Both DNA and membrane should move during simulation
- Check output files in VTU_F/ and TrajTSI/ directories

## Adjusting Parameters

- **Repulsion strength**: Change `1.0` in `RepulsionBetweenVerticesAndDNABeads` line
- **Cutoff radius**: Change `5.0` in `RepulsionBetweenVerticesAndDNABeads` line
- **Bond stiffness**: Change `10.0` in `dna_bonds.txt`
- **Equilibrium distance**: Change `2.0` in `dna_bonds.txt`

## Visualization

### For VMD (GROMACS .gro format):
```bash
# Convert a single frame
../CNV -in TrajTSI/dts1.tsi -o output.gro

# Or use the conversion script to convert all frames
./convert_to_vmd.sh

# Then view in VMD:
vmd VMD_format/dts1.gro
# Or load multiple frames:
vmd VMD_format/*.gro
```

### For ParaView (VTU format):
```bash
# Convert to VTU format
../CNV -in TrajTSI/dts1.tsi -o output.vtu

# Then open in ParaView
paraview output.vtu
```

**Note**: 
- `.gro` files are GROMACS format and work with VMD
- `.vtu` files are ParaView format
- `.tsi` files are the native FreeDTS trajectory format

