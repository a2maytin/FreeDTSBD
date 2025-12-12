#!/bin/bash
# Fast script to generate dna_bonds.txt with global parameters
# Usage: ./generate_dna_bonds.sh [start_idx] [num_beads] [k] [l0] [k_angle] [theta0_deg]

START_IDX=${1:-1252}
NUM_BEADS=${2:-5000}
K=${3:-10}
L0=${4:-0.472}
K_ANGLE=${5:-0}
THETA0_DEG=${6:-180}

python3 << EOF
start_idx = $START_IDX
num_beads = $NUM_BEADS
k = $K
l0 = $L0
k_angle = $K_ANGLE
theta0_deg = $THETA0_DEG

print('; DNA bonds file - circular topology')
print('; Format: vertex_id1 vertex_id2 [k l0] (if k and l0 not specified, use global values)')
print(f'; DNA beads start at index {start_idx} in activeVertices')
print(f'; Circular chain: {num_beads} beads')
print(f'; Global bond parameters: k={k}, l0={l0} nm')
if k_angle > 0:
    print(f'; Global angle parameters: k_angle={k_angle}, theta0={theta0_deg} deg')
print()

for i in range(num_beads):
    id1 = start_idx + i
    id2 = start_idx + ((i + 1) % num_beads)
    print(f'{id1}  {id2}')
EOF

