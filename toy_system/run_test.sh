#!/bin/bash
# Simple script to run the toy system test

echo "Running toy system test..."
echo ""

# Try to find DTS executable
# First check if DTS is in PATH
if command -v DTS &> /dev/null; then
    DTS_CMD="DTS"
elif [ -f "../DTS" ]; then
    # DTS is in parent directory (project root)
    DTS_CMD="../DTS"
elif [ -f "./DTS" ]; then
    # DTS is in current directory
    DTS_CMD="./DTS"
else
    echo "Error: DTS executable not found!"
    echo "Please either:"
    echo "  1. Add DTS to your PATH, or"
    echo "  2. Place DTS in the parent directory (../DTS), or"
    echo "  3. Place DTS in the current directory (./DTS)"
    exit 1
fi

echo "Using DTS: $DTS_CMD"
echo ""

# Run the simulation
$DTS_CMD -in input.dts -top topology.top

echo ""
echo "Simulation complete! Check output files:"
echo "  - VTU_F/ directory for visualization files"
echo "  - TrajTSI/ directory for trajectory files"
echo "  - dts.log for simulation log"

