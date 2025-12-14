#!/usr/bin/env python3
"""
Generate angles file from bonds file.
For each triplet of consecutive vertices in the bonds (excluding SMC bonds),
creates an angle entry.

Usage:
    python generate_angles_from_bonds.py <bonds_file> [output_angles_file] [k_angle] [theta0]
    
If output_angles_file is not specified, it will be bonds_file with .txt replaced by _angles.txt
"""

import sys
import os

def parse_bonds_file(bonds_file):
    """Parse bonds file and return bonds list and SMC bonds set."""
    bonds = []
    smc_bonds = set()
    default_k_angle = 20.0
    default_theta0 = 0.0
    
    with open(bonds_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                # Check for default parameters in comments
                if 'k_angle' in line.lower():
                    parts = line.split('=')
                    if len(parts) > 1:
                        try:
                            default_k_angle = float(parts[1].split()[0])
                        except:
                            pass
                if 'theta0' in line.lower():
                    parts = line.split('=')
                    if len(parts) > 1:
                        try:
                            default_theta0 = float(parts[1].split()[0])
                        except:
                            pass
                continue
            
            # Strip inline comments
            if ';' in line:
                comment_pos = line.find(';')
                is_smc = 'SMC' in line[comment_pos:]
                line = line[:comment_pos].strip()
                if is_smc:
                    # This is an SMC bond
                    parts = line.split()
                    if len(parts) >= 2:
                        v1, v2 = int(parts[0]), int(parts[1])
                        smc_bonds.add((min(v1, v2), max(v1, v2)))
                    continue
            
            parts = line.split()
            if len(parts) >= 2:
                v1, v2 = int(parts[0]), int(parts[1])
                bonds.append((v1, v2))
    
    return bonds, smc_bonds, default_k_angle, default_theta0

def build_adjacency(bonds, smc_bonds):
    """Build adjacency map from bonds, excluding SMC bonds."""
    adjacency = {}
    for v1, v2 in bonds:
        bond_pair = (min(v1, v2), max(v1, v2))
        if bond_pair not in smc_bonds:
            if v1 not in adjacency:
                adjacency[v1] = set()
            if v2 not in adjacency:
                adjacency[v2] = set()
            adjacency[v1].add(v2)
            adjacency[v2].add(v1)
    return adjacency

def find_chains(bonds, smc_bonds, adjacency):
    """Find all chains by following bonds in order, preserving sequence."""
    visited_bonds = set()
    chains = []
    
    # Build a map from vertex to list of bonds it's part of (for easier lookup)
    vertex_to_bonds = {}
    for i, (v1, v2) in enumerate(bonds):
        bond_pair = (min(v1, v2), max(v1, v2))
        if bond_pair not in smc_bonds:
            if v1 not in vertex_to_bonds:
                vertex_to_bonds[v1] = []
            if v2 not in vertex_to_bonds:
                vertex_to_bonds[v2] = []
            vertex_to_bonds[v1].append(i)
            vertex_to_bonds[v2].append(i)
    
    # Follow bonds in sequence to build chains
    for i, (v1, v2) in enumerate(bonds):
        bond_pair = (min(v1, v2), max(v1, v2))
        if bond_pair in smc_bonds or i in visited_bonds:
            continue
        
        # Start a new chain from this bond
        chain = []
        visited_bonds.add(i)
        
        # Build chain by following connected bonds
        # Start from v1 and go in direction of v2
        current_vertex = v2
        prev_vertex = v1
        chain.append(v1)
        chain.append(v2)
        
        # Follow the chain forward
        while True:
            next_vertex = None
            # Find the next bond that connects to current_vertex but not to prev_vertex
            if current_vertex in vertex_to_bonds:
                for bond_idx in vertex_to_bonds[current_vertex]:
                    if bond_idx in visited_bonds:
                        continue
                    bv1, bv2 = bonds[bond_idx]
                    if bv1 == current_vertex and bv2 != prev_vertex:
                        next_vertex = bv2
                        visited_bonds.add(bond_idx)
                        break
                    elif bv2 == current_vertex and bv1 != prev_vertex:
                        next_vertex = bv1
                        visited_bonds.add(bond_idx)
                        break
            
            if next_vertex is None:
                break
            
            chain.append(next_vertex)
            prev_vertex = current_vertex
            current_vertex = next_vertex
        
        # Check if chain is circular (first and last vertices are connected)
        if len(chain) >= 3:
            first_vertex = chain[0]
            last_vertex = chain[-1]
            # Check if chain loops back to start (first == last means circular)
            is_circular = (first_vertex == last_vertex)
            # Also check if first and last are connected (circular chain)
            if not is_circular and first_vertex in adjacency:
                is_circular = last_vertex in adjacency[first_vertex]
            
            # Also check if we have a bond connecting them in the original bonds list
            if not is_circular:
                for v1, v2 in bonds:
                    bond_pair = (min(v1, v2), max(v1, v2))
                    if bond_pair not in smc_bonds:
                        if (v1 == first_vertex and v2 == last_vertex) or (v1 == last_vertex and v2 == first_vertex):
                            is_circular = True
                            break
            
            # If chain loops back to start, remove the duplicate last vertex
            if is_circular and first_vertex == last_vertex and len(chain) > 1:
                chain = chain[:-1]  # Remove duplicate last vertex
            
            if is_circular:
                # For circular chains, we have all angles (one per vertex)
                chains.append((chain, True))
            else:
                # For linear chains, we need at least 3 vertices
                chains.append((chain, False))
    
    return chains

def generate_angles(chains, k_angle, theta0):
    """Generate angles for all triplets in chains."""
    angles = []
    
    for chain, is_circular in chains:
        if len(chain) < 3:
            continue
        
        if is_circular:
            # Circular chain: create angles for all triplets including wraparound
            for i in range(len(chain)):
                prev_idx = (i - 1) % len(chain)
                curr_idx = i
                next_idx = (i + 1) % len(chain)
                
                v1 = chain[prev_idx]
                v2 = chain[curr_idx]
                v3 = chain[next_idx]
                
                angles.append((v1, v2, v3, k_angle, theta0))
        else:
            # Linear chain: create angles for triplets (skip first and last)
            for i in range(1, len(chain) - 1):
                v1 = chain[i - 1]
                v2 = chain[i]
                v3 = chain[i + 1]
                
                angles.append((v1, v2, v3, k_angle, theta0))
    
    return angles

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_angles_from_bonds.py <bonds_file> [output_angles_file] [k_angle] [theta0]")
        sys.exit(1)
    
    bonds_file = sys.argv[1]
    
    # Determine output file
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        base_name = os.path.splitext(bonds_file)[0]
        output_file = base_name + "_angles.txt"
    
    # Get k_angle and theta0 from command line or use defaults from file
    bonds, smc_bonds, default_k_angle, default_theta0 = parse_bonds_file(bonds_file)
    
    k_angle = float(sys.argv[3]) if len(sys.argv) > 3 else default_k_angle
    theta0 = float(sys.argv[4]) if len(sys.argv) > 4 else default_theta0
    
    print(f"Reading bonds from: {bonds_file}")
    print(f"Found {len(bonds)} bonds, {len(smc_bonds)} SMC bonds")
    print(f"Using k_angle={k_angle}, theta0={theta0}")
    
    # Build adjacency and find chains
    adjacency = build_adjacency(bonds, smc_bonds)
    chains = find_chains(bonds, smc_bonds, adjacency)
    
    print(f"Found {len(chains)} chains")
    for i, (chain, is_circular) in enumerate(chains):
        chain_type = "circular" if is_circular else "linear"
        print(f"  Chain {i+1}: {len(chain)} vertices ({chain_type})")
    
    # Generate angles
    angles = generate_angles(chains, k_angle, theta0)
    
    print(f"Generated {len(angles)} angles")
    print(f"Writing to: {output_file}")
    
    # Write angles file
    with open(output_file, 'w') as f:
        f.write("; Angles file (generated from bonds file)\n")
        f.write(f"; k_angle = {k_angle}\n")
        f.write(f"; theta0 = {theta0}\n")
        f.write("; Format: v1 v2 v3 [k_angle] [theta0]\n")
        f.write("; (k_angle and theta0 are optional if specified in header)\n\n")
        
        for v1, v2, v3, k, t0 in angles:
            f.write(f"{v1} {v2} {v3}\n")
    
    print("Done!")

if __name__ == "__main__":
    main()
