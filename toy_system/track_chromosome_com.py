#!/usr/bin/env python3
"""
Track center of mass of two chromosomes from TSI trajectory files.

Chromosome 1: beads 0-1999 (2000 beads)
Chromosome 2: beads 2000-3999 (2000 beads)
"""

import os
import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

def parse_tsi_file(filename):
    """Parse a TSI file and return vertex positions as a dictionary {id: (x, y, z)}."""
    vertices = {}
    box = None
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('version'):
            i += 1
            continue
        elif line.startswith('box'):
            # Extract box dimensions
            parts = line[3:].split()
            if len(parts) >= 3:
                box = (float(parts[0]), float(parts[1]), float(parts[2]))
            i += 1
            continue
        elif line.startswith('vertex'):
            # Extract number of vertices from the same line
            parts = line[6:].strip().split()  # Skip "vertex" prefix
            if len(parts) > 0:
                n_vertices = int(parts[0])
            else:
                i += 1
                if i >= len(lines):
                    break
                n_vertices = int(lines[i].strip())
            i += 1
            
            # Read vertex positions
            for j in range(n_vertices):
                if i >= len(lines):
                    break
                parts = lines[i].strip().split()
                if len(parts) >= 4:
                    vid = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    vertices[vid] = (x, y, z)
                i += 1
            continue
        else:
            i += 1
    
    return vertices, box

def calculate_com(vertices, start_id, end_id):
    """Calculate center of mass for vertices with IDs from start_id to end_id (inclusive)."""
    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0
    count = 0
    
    for vid in range(start_id, end_id + 1):
        if vid in vertices:
            x, y, z = vertices[vid]
            x_sum += x
            y_sum += y
            z_sum += z
            count += 1
    
    if count == 0:
        return None
    
    return (x_sum / count, y_sum / count, z_sum / count)

def main():
    # Find all TSI files in TrajTSI directory
    tsi_dir = 'TrajTSI'
    if not os.path.exists(tsi_dir):
        print(f"Error: Directory {tsi_dir} not found!")
        sys.exit(1)
    
    # Get all TSI files and sort them numerically
    tsi_files = glob.glob(os.path.join(tsi_dir, 'dts*.tsi'))
    
    def extract_number(filename):
        match = re.search(r'dts(\d+)\.tsi', filename)
        return int(match.group(1)) if match else 0
    
    tsi_files.sort(key=extract_number)
    
    if len(tsi_files) == 0:
        print(f"Error: No TSI files found in {tsi_dir}!")
        sys.exit(1)
    
    print(f"Found {len(tsi_files)} TSI files")
    print("Tracking center of mass for:")
    print("  Chromosome 1: beads 0-1999")
    print("  Chromosome 2: beads 2000-3999")
    print()
    
    # Storage for results
    timesteps = []
    com1_x, com1_y, com1_z = [], [], []
    com2_x, com2_y, com2_z = [], [], []
    com1_initial = None
    com2_initial = None
    
    # Process each TSI file
    for tsi_file in tsi_files:
        step = extract_number(tsi_file)
        vertices, box = parse_tsi_file(tsi_file)
        
        if len(vertices) == 0:
            print(f"Warning: No vertices found in {tsi_file}")
            continue
        
        # Calculate center of mass for each chromosome
        com1 = calculate_com(vertices, 0, 1999)
        com2 = calculate_com(vertices, 2000, 3999)
        
        if com1 is None or com2 is None:
            print(f"Warning: Could not calculate COM for step {step}")
            continue
        
        timesteps.append(step)
        com1_x.append(com1[0])
        com1_y.append(com1[1])
        com1_z.append(com1[2])
        com2_x.append(com2[0])
        com2_y.append(com2[1])
        com2_z.append(com2[2])
        
        # Store initial COM for reference
        if com1_initial is None:
            com1_initial = com1
            com2_initial = com2
        
        # Calculate displacement from initial position
        disp1 = np.sqrt((com1[0] - com1_initial[0])**2 + 
                       (com1[1] - com1_initial[1])**2 + 
                       (com1[2] - com1_initial[2])**2)
        disp2 = np.sqrt((com2[0] - com2_initial[0])**2 + 
                       (com2[1] - com2_initial[1])**2 + 
                       (com2[2] - com2_initial[2])**2)
        
        if step % 10 == 0 or step == timesteps[0] or step == timesteps[-1]:
            print(f"Step {step:6d}: COM1 = ({com1[0]:8.3f}, {com1[1]:8.3f}, {com1[2]:8.3f}), "
                  f"displacement = {disp1:8.3f}")
            print(f"          COM2 = ({com2[0]:8.3f}, {com2[1]:8.3f}, {com2[2]:8.3f}), "
                  f"displacement = {disp2:8.3f}")
    
    # Convert to numpy arrays for easier analysis
    timesteps = np.array(timesteps)
    com1_x = np.array(com1_x)
    com1_y = np.array(com1_y)
    com1_z = np.array(com1_z)
    com2_x = np.array(com2_x)
    com2_y = np.array(com2_y)
    com2_z = np.array(com2_z)
    
    # Calculate displacements from initial positions
    disp1 = np.sqrt((com1_x - com1_initial[0])**2 + 
                   (com1_y - com1_initial[1])**2 + 
                   (com1_z - com1_initial[2])**2)
    disp2 = np.sqrt((com2_x - com2_initial[0])**2 + 
                   (com2_y - com2_initial[1])**2 + 
                   (com2_z - com2_initial[2])**2)
    
    # Print summary statistics
    print()
    print("=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    print(f"Initial COM1: ({com1_initial[0]:.6f}, {com1_initial[1]:.6f}, {com1_initial[2]:.6f})")
    print(f"Initial COM2: ({com2_initial[0]:.6f}, {com2_initial[1]:.6f}, {com2_initial[2]:.6f})")
    print()
    print(f"Final COM1:   ({com1_x[-1]:.6f}, {com1_y[-1]:.6f}, {com1_z[-1]:.6f})")
    print(f"Final COM2:   ({com2_x[-1]:.6f}, {com2_y[-1]:.6f}, {com2_z[-1]:.6f})")
    print()
    print(f"Max displacement COM1: {np.max(disp1):.6f}")
    print(f"Max displacement COM2: {np.max(disp2):.6f}")
    print(f"RMS displacement COM1: {np.sqrt(np.mean(disp1**2)):.6f}")
    print(f"RMS displacement COM2: {np.sqrt(np.mean(disp2**2)):.6f}")
    print()
    print(f"Final displacement COM1: {disp1[-1]:.6f}")
    print(f"Final displacement COM2: {disp2[-1]:.6f}")
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: COM positions over time (X component)
    axes[0, 0].plot(timesteps, com1_x, 'b-', label='Chromosome 1', linewidth=1.5)
    axes[0, 0].plot(timesteps, com2_x, 'r-', label='Chromosome 2', linewidth=1.5)
    axes[0, 0].axhline(com1_initial[0], color='b', linestyle='--', alpha=0.5, label='Initial COM1')
    axes[0, 0].axhline(com2_initial[0], color='r', linestyle='--', alpha=0.5, label='Initial COM2')
    axes[0, 0].set_xlabel('Step')
    axes[0, 0].set_ylabel('X position')
    axes[0, 0].set_title('Center of Mass X Position')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: COM positions over time (Y component)
    axes[0, 1].plot(timesteps, com1_y, 'b-', label='Chromosome 1', linewidth=1.5)
    axes[0, 1].plot(timesteps, com2_y, 'r-', label='Chromosome 2', linewidth=1.5)
    axes[0, 1].axhline(com1_initial[1], color='b', linestyle='--', alpha=0.5, label='Initial COM1')
    axes[0, 1].axhline(com2_initial[1], color='r', linestyle='--', alpha=0.5, label='Initial COM2')
    axes[0, 1].set_xlabel('Step')
    axes[0, 1].set_ylabel('Y position')
    axes[0, 1].set_title('Center of Mass Y Position')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: COM positions over time (Z component)
    axes[1, 0].plot(timesteps, com1_z, 'b-', label='Chromosome 1', linewidth=1.5)
    axes[1, 0].plot(timesteps, com2_z, 'r-', label='Chromosome 2', linewidth=1.5)
    axes[1, 0].axhline(com1_initial[2], color='b', linestyle='--', alpha=0.5, label='Initial COM1')
    axes[1, 0].axhline(com2_initial[2], color='r', linestyle='--', alpha=0.5, label='Initial COM2')
    axes[1, 0].set_xlabel('Step')
    axes[1, 0].set_ylabel('Z position')
    axes[1, 0].set_title('Center of Mass Z Position')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Displacement from initial position
    axes[1, 1].plot(timesteps, disp1, 'b-', label='Chromosome 1', linewidth=1.5)
    axes[1, 1].plot(timesteps, disp2, 'r-', label='Chromosome 2', linewidth=1.5)
    axes[1, 1].set_xlabel('Step')
    axes[1, 1].set_ylabel('Displacement from initial position')
    axes[1, 1].set_title('Center of Mass Displacement')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('chromosome_com_tracking.png', dpi=150)
    print()
    print(f"Plot saved to: chromosome_com_tracking.png")
    
    # Save data to text file
    output_file = 'chromosome_com_data.txt'
    with open(output_file, 'w') as f:
        f.write("# Step COM1_x COM1_y COM1_z COM2_x COM2_y COM2_z disp1 disp2\n")
        for i in range(len(timesteps)):
            f.write(f"{timesteps[i]} {com1_x[i]:.6f} {com1_y[i]:.6f} {com1_z[i]:.6f} "
                   f"{com2_x[i]:.6f} {com2_y[i]:.6f} {com2_z[i]:.6f} "
                   f"{disp1[i]:.6f} {disp2[i]:.6f}\n")
    print(f"Data saved to: {output_file}")

if __name__ == '__main__':
    main()
