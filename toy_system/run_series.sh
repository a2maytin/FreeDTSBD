#!/bin/bash
# Run a series of simulations with incrementing volume and area coupling parameters
# Each simulation restarts from the previous one's checkpoint (dts.res)

# Configuration
DTS_EXEC="../DTS"
INPUT_TEMPLATE="input.dts"
TOPOLOGY="topology.top"
NUM_STEPS=21
STEPS_PER_SIM=500000
TOTAL_STEPS=$((NUM_STEPS * STEPS_PER_SIM))

# Starting values
START_GAMMAV=0.7
START_AREA=0.34

# Ending values
END_GAMMAV=0.2
END_AREA=0.34

# Calculate increment per step 
GAMMAV_INCREMENT=$(awk "BEGIN {print ($END_GAMMAV - $START_GAMMAV) / ($NUM_STEPS - 1)}")
AREA_INCREMENT=$(awk "BEGIN {print ($END_AREA - $START_AREA) / ($NUM_STEPS - 1)}")

echo "=========================================="
echo "Running series of $NUM_STEPS simulations"
echo "Each simulation: $STEPS_PER_SIM steps"
echo "Total steps across all simulations: $TOTAL_STEPS"
echo "Starting: gammaV=$START_GAMMAV, area=$START_AREA"
echo "=========================================="
echo ""

# Check if DTS executable exists
if [ ! -f "$DTS_EXEC" ]; then
    echo "ERROR: DTS executable not found: $DTS_EXEC"
    exit 1
fi

# Check if input template exists
if [ ! -f "$INPUT_TEMPLATE" ]; then
    echo "ERROR: Input template not found: $INPUT_TEMPLATE"
    exit 1
fi

# Check if topology exists
if [ ! -f "$TOPOLOGY" ]; then
    echo "ERROR: Topology file not found: $TOPOLOGY"
    exit 1
fi

# Run each simulation
for step in $(seq 1 $NUM_STEPS); do
    # Calculate current parameter values using awk
    gammaV=$(awk "BEGIN {printf \"%.2f\", $START_GAMMAV + ($step - 1) * $GAMMAV_INCREMENT}")
    area=$(awk "BEGIN {printf \"%.2f\", $START_AREA + ($step - 1) * $AREA_INCREMENT}")
    
    echo "============================================================"
    echo "Step $step/$NUM_STEPS: gammaV=$gammaV, area=$area"
    echo "============================================================"
    
    # Create modified input file
    cp "$INPUT_TEMPLATE" input_current.dts
    
    # Calculate step range for this simulation
    # Step 1: 1 to 5000, Step 2: 5001 to 10000, Step 3: 10001 to 15000, etc.
    INITIAL_STEP=$(((step - 1) * STEPS_PER_SIM + 1))
    FINAL_STEP=$((step * STEPS_PER_SIM))
    
    # Update Set_Steps to have correct initial and final steps for this segment
    # Format: Set_Steps = <initial> <final>
    # When restarting, FreeDTS will continue from the restart step to FINAL_STEP
    sed -i "s/^Set_Steps = [0-9]* [0-9]*/Set_Steps = $INITIAL_STEP $FINAL_STEP/" input_current.dts
    
    # Update VolumeCoupling (line with "VolumeCoupling = SecondOrder ...")
    sed -i "s/^VolumeCoupling = SecondOrder [0-9.]* [0-9.]* [0-9.]*/VolumeCoupling = SecondOrder 0.0 10000 $gammaV/" input_current.dts
    
    # Update TotalAreaCoupling (line with "TotalAreaCoupling = HarmonicPotential ...")
    sed -i "s/^TotalAreaCoupling = HarmonicPotential [0-9.]* [0-9.]*/TotalAreaCoupling = HarmonicPotential 1000 $area/" input_current.dts
    
    # Check if restart file exists (from previous simulation)
    if [ -f "dts.res" ] && [ $step -gt 1 ]; then
        RESTART_ARG="-restart dts.res"
        echo "Restarting from dts.res (will run from previous checkpoint to step $FINAL_STEP)"
    else
        RESTART_ARG=""
        echo "Starting new simulation (steps $INITIAL_STEP to $FINAL_STEP)"
    fi
    
    # Run simulation - construct command properly
    CMD="$DTS_EXEC -in input_current.dts -top $TOPOLOGY"
    if [ -n "$RESTART_ARG" ]; then
        CMD="$CMD $RESTART_ARG"
    fi
    CMD="$CMD -nt 15"
    
    echo "Command: $CMD"
    $CMD
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "✓ Step $step completed successfully"
    else
        echo "✗ Step $step exited with code $EXIT_CODE"
        echo "Stopping series due to error"
        exit $EXIT_CODE
    fi
    
    echo ""
done

# Clean up temporary input file
rm -f input_current.dts

echo "=========================================="
echo "All $NUM_STEPS simulations completed!"
echo "=========================================="
