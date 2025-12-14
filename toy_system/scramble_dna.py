#!/usr/bin/env python3
"""
Scramble DNA vertex indices to test if array ordering causes spatial drift.
Keeps the same coordinates but randomizes the order in the array.
"""

import random
import sys

def scramble_dna_q(input_file, output_file, seed=42):
    """Scramble vertex indices in .q file"""
    random.seed(seed)  # For reproducibility
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    box_line = lines[0].strip()
    num_vertices = int(lines[1].strip())
    
    # Parse vertices
    vertices = []
    for i in range(2, 2 + num_vertices):
        parts = lines[i].strip().split()
        vid = int(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        group = int(parts[4]) if len(parts) > 4 else 0
        vertices.append((vid, x, y, z, group))
    
    # Create random permutation: perm[i] = old vertex ID at new position i
    # Example: perm[0] = 795 means old vertex 795 is at new position 0
    perm = list(range(num_vertices))
    random.shuffle(perm)
    
    # Create mapping: old_vertex_id -> new_vertex_id (dictionary)
    old_to_new_dict = {old_id: new_pos for new_pos, old_id in enumerate(perm)}
    
    # Write scrambled file
    with open(output_file, 'w') as f:
        f.write(box_line + '\n')
        f.write(f'{num_vertices}\n')
        
        # Write vertices in new order (new index, but keep original coordinates)
        for new_idx in range(num_vertices):
            old_idx = perm[new_idx]  # Get old vertex ID at this new position
            vid, x, y, z, group = vertices[old_idx]
            # Write with new index, but keep original coordinates
            f.write(f'{new_idx}  {x}  {y}  {z}  {group}\n')
        
        # Write triangle count (should be 0 for DNA)
        if len(lines) > 2 + num_vertices:
            f.write(lines[2 + num_vertices])
        else:
            f.write('0\n')
    
    return old_to_new_dict, perm

def scramble_bonds(input_file, output_file, old_to_new_dict):
    """Update bonds file with new vertex indices"""
    # old_to_new_dict is a dictionary: old_vertex_id -> new_vertex_id
    # Use it directly to map bond indices
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    with open(output_file, 'w') as f:
        # Write header (first 8 lines are comments)
        header_end = 0
        for i, line in enumerate(lines):
            if line.strip().startswith(';') or not line.strip():
                f.write(line)
                header_end = i + 1
            else:
                break
        
        # Process bonds
        for line in lines[header_end:]:
            line = line.strip()
            if not line:
                f.write('\n')
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    v1_old = int(parts[0])
                    v2_old = int(parts[1])
                    
                    # Map to new indices
                    v1_new = old_to_new_dict.get(v1_old, v1_old)
                    v2_new = old_to_new_dict.get(v2_old, v2_old)
                    
                    # Write bond with new indices
                    if len(parts) > 2:
                        # Include additional parameters if present
                        f.write(f'{v1_new} {v2_new} {" ".join(parts[2:])}\n')
                    else:
                        f.write(f'{v1_new} {v2_new}\n')
                except ValueError:
                    # Not a bond line, just copy it
                    f.write(line + '\n')
            else:
                f.write(line + '\n')

if __name__ == '__main__':
    input_q = 'dna_single_chromosome.q'
    output_q = 'dna_single_chromosome_scrambled.q'
    input_bonds = 'dna_bonds_single_chromosome.txt'
    output_bonds = 'dna_bonds_single_chromosome_scrambled.txt'
    
    print(f"Reading {input_q}...")
    old_to_new_dict, perm = scramble_dna_q(input_q, output_q, seed=42)
    print(f"Created {output_q} with scrambled vertex indices")
    
    print(f"Reading {input_bonds}...")
    scramble_bonds(input_bonds, output_bonds, old_to_new_dict)
    print(f"Created {output_bonds} with updated bond indices")
    
    print("\nScrambling complete!")
    print(f"  Original: vertices 0-1999 in chain order")
    print(f"  Scrambled: same coordinates, random array order")
    print(f"  Example mapping: old vertex 0 -> new position {old_to_new_dict[0]}")
    print(f"                   old vertex 1 -> new position {old_to_new_dict[1]}")
