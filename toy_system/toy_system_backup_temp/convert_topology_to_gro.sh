#!/bin/bash
# Script to convert topology.top system to .gro format
# Reads topology.top, loads all .q files, and writes a combined .gro file

if [ $# -lt 1 ]; then
    echo "Usage: $0 <output.gro> [topology.top]"
    echo "  output.gro: Output GROMACS .gro file"
    echo "  topology.top: Topology file (default: topology.top)"
    exit 1
fi

OUTPUT_GRO="$1"
TOPOLOGY_FILE="${2:-topology.top}"

if [ ! -f "$TOPOLOGY_FILE" ]; then
    echo "Error: Topology file not found: $TOPOLOGY_FILE"
    exit 1
fi

echo "Converting topology system to .gro format..."
echo "Topology file: $TOPOLOGY_FILE"
echo "Output file: $OUTPUT_GRO"
echo ""

# Temporary file for processing
TMP_COORDS=$(mktemp)
TMP_BOX=$(mktemp)

# Read topology file and process each .q file
total_vertices=0
max_box_x=0
max_box_y=0
max_box_z=0

while IFS= read -r line || [ -n "$line" ]; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    
    # Parse: filename group_id
    qfile=$(echo "$line" | awk '{print $1}')
    group_id=$(echo "$line" | awk '{print $2}')
    
    if [ -z "$qfile" ]; then
        continue
    fi
    
    if [ ! -f "$qfile" ]; then
        echo "Warning: .q file not found: $qfile (skipping)"
        continue
    fi
    
    echo "Processing: $qfile (group $group_id)"
    
    # Read .q file
    # Format: 
    # Line 1: box_x box_y box_z
    # Line 2: num_vertices
    # Lines 3+: id x y z domain
    
    {
        read -r box_line
        read -r num_vertices
        
        # Extract box dimensions
        box_x=$(echo "$box_line" | awk '{print $1}')
        box_y=$(echo "$box_line" | awk '{print $2}')
        box_z=$(echo "$box_line" | awk '{print $3}')
        
        # Update max box dimensions
        if (( $(echo "$box_x > $max_box_x" | bc -l) )); then max_box_x=$box_x; fi
        if (( $(echo "$box_y > $max_box_y" | bc -l) )); then max_box_y=$box_y; fi
        if (( $(echo "$box_z > $max_box_z" | bc -l) )); then max_box_z=$box_z; fi
        
        # Read vertices
        vertex_count=0
        while IFS= read -r vertex_line && [ $vertex_count -lt $num_vertices ]; do
            # Skip empty lines
            [[ -z "$vertex_line" ]] && continue
            
            # Parse: id x y z domain
            id=$(echo "$vertex_line" | awk '{print $1}')
            x=$(echo "$vertex_line" | awk '{print $2}')
            y=$(echo "$vertex_line" | awk '{print $3}')
            z=$(echo "$vertex_line" | awk '{print $4}')
            
            # Determine atom type based on group
            # Format matches FreeDTS: resid resname atomname atomid x y z vx vy vz
            # resname = "Ver", atomname = "C" (membrane) or "D" (DNA)
            if [ "$group_id" = "0" ]; then
                resname="Ver"
                atomname="C"  # Membrane vertex
            else
                resname="Ver"
                atomname="D"  # DNA bead
            fi
            
            # Write to temp file (format: resid resname atomname atomid x y z vx vy vz)
            # resid = group_id, atomid = total_vertices, velocities = 0.0
            printf "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" \
                $group_id "$resname" "$atomname" $total_vertices $x $y $z 0.0 0.0 0.0 >> "$TMP_COORDS"
            
            total_vertices=$((total_vertices + 1))
            vertex_count=$((vertex_count + 1))
        done
    } < "$qfile"
    
done < "$TOPOLOGY_FILE"

# Write .gro file
{
    # Header (title line)
    echo "Generated from topology.top by FreeDTS"
    
    # Number of atoms
    printf "%5d\n" $total_vertices
    
    # Coordinates (already formatted in TMP_COORDS)
    cat "$TMP_COORDS"
    
    # Box dimensions (in nm)
    # Format: v1x v2y v3z v1y v1z v2x v2z v3x v3y (triclinic, but we use rectangular)
    printf "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n" \
        $max_box_x $max_box_y $max_box_z \
        0.0 0.0 0.0 \
        0.0 0.0 0.0
    
} > "$OUTPUT_GRO"

# Cleanup
rm -f "$TMP_COORDS" "$TMP_BOX"

echo ""
echo "Conversion complete!"
echo "  Total vertices: $total_vertices"
echo "  Box dimensions: ${max_box_x} x ${max_box_y} x ${max_box_z} nm"
echo "  Output file: $OUTPUT_GRO"
echo ""
echo "To view in VMD:"
echo "  vmd $OUTPUT_GRO"

