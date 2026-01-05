#!/bin/bash
# Script to convert trajectory files to VMD-compatible format

echo "Converting trajectory files to VMD format (.gro)..."
echo ""

# Find CNV executable
if command -v CNV &> /dev/null; then
    CNV_CMD="CNV"
elif [ -f "../CNV" ]; then
    CNV_CMD="../CNV"
elif [ -f "./CNV" ]; then
    CNV_CMD="./CNV"
else
    echo "Error: CNV executable not found!"
    echo "Please compile the code first: bash ../compile.sh"
    exit 1
fi

# Create output directory
mkdir -p VMD_format

# Convert each trajectory file
if [ -d "TrajTSI" ]; then
    # Sort files numerically (increasing order)
    for tsi_file in $(ls TrajTSI/*.tsi 2>/dev/null | sort -V); do
        if [ -f "$tsi_file" ]; then
            base_name=$(basename "$tsi_file" .tsi)
            output_file="VMD_format/${base_name}.gro"
            echo "Converting $tsi_file -> $output_file"
            $CNV_CMD -in "$tsi_file" -o "$output_file" 2>&1 | grep -v "^$" || true
        fi
    done
    echo ""
    echo "Conversion complete! VMD-compatible files are in VMD_format/ directory"
    echo ""
    echo "To view in VMD:"
    echo "  vmd VMD_format/dts1.gro"
    echo ""
    echo "Or load multiple frames:"
    echo "  vmd VMD_format/*.gro"
else
    echo "Error: TrajTSI directory not found!"
    echo "Please run the simulation first."
    exit 1
fi

