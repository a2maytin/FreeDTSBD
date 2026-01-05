#!/usr/bin/env python3
"""
Script to shift x, y, z coordinates in dna.q file by specified amounts.
"""

import sys
import os

def shift_coordinates(input_file, output_file, shift_x=0.0, shift_y=0.0, shift_z=0.0):
    """
    Shift x, y, z coordinates in dna.q file by specified values.
    
    Format:
    Line 1: Box dimensions (3 floats)
    Line 2: Number of atoms
    Lines 3+: Atom data (id x y z 0)
    
    Args:
        input_file: Input file path
        output_file: Output file path
        shift_x: Shift amount for x coordinates (default 0.0)
        shift_y: Shift amount for y coordinates (default 0.0)
        shift_z: Shift amount for z coordinates (default 0.0)
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 2:
        print(f"Error: File {input_file} has too few lines")
        return False
    
    # Parse header
    box_dims = lines[0].strip().split()
    num_atoms = int(lines[1].strip())
    
    print(f"Box dimensions: {box_dims}")
    print(f"Number of atoms: {num_atoms}")
    print(f"Shifting coordinates:")
    print(f"  x: {shift_x:+f}")
    print(f"  y: {shift_y:+f}")
    print(f"  z: {shift_z:+f}")
    
    # Process atom lines
    output_lines = []
    output_lines.append(lines[0])  # Box dimensions (unchanged)
    output_lines.append(lines[1])  # Number of atoms (unchanged)
    
    atoms_shifted = 0
    for i in range(2, len(lines)):
        line = lines[i].strip()
        if not line:
            output_lines.append('\n')
            continue
        
        parts = line.split()
        if len(parts) >= 4:
            atom_id = parts[0]
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            rest = ' '.join(parts[4:]) if len(parts) > 4 else ''
            
            # Shift coordinates
            x_new = x + shift_x
            y_new = y + shift_y
            z_new = z + shift_z
            
            # Format output line (preserve original formatting as much as possible)
            if rest:
                output_lines.append(f"{atom_id}  {x_new:15.10f}  {y_new:15.10f}  {z_new:15.10f}  {rest}\n")
            else:
                output_lines.append(f"{atom_id}  {x_new:15.10f}  {y_new:15.10f}  {z_new:15.10f}  0\n")
            atoms_shifted += 1
        else:
            # Keep line as-is if it doesn't match expected format
            output_lines.append(line + '\n')
    
    # Write output
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Shifted {atoms_shifted} atoms")
    print(f"Output written to {output_file}")
    return True

if __name__ == "__main__":
    input_file = "dna.q"
    output_file = "dna.q"
    shift_x = 0.0
    shift_y = 0.0
    shift_z = 0.0
    
    # Allow command line arguments
    # Usage: python shift_dna_x.py [input_file] [output_file] [shift_x] [shift_y] [shift_z]
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    if len(sys.argv) > 3:
        shift_x = float(sys.argv[3])
    if len(sys.argv) > 4:
        shift_y = float(sys.argv[4])
    if len(sys.argv) > 5:
        shift_z = float(sys.argv[5])
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    
    # If output is same as input, create backup
    if input_file == output_file:
        backup_file = input_file + ".backup"
        print(f"Creating backup: {backup_file}")
        import shutil
        shutil.copy2(input_file, backup_file)
    
    success = shift_coordinates(input_file, output_file, shift_x, shift_y, shift_z)
    if not success:
        sys.exit(1)
