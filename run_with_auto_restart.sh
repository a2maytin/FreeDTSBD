#!/bin/bash
# Auto-restart wrapper for FreeDTS that restarts on segmentation fault
# Usage: ./run_with_auto_restart.sh [DTS arguments]

# Configuration
MAX_RESTARTS=10  # Maximum number of restart attempts
RESTART_DELAY=1  # Seconds to wait before restarting

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DTS_BINARY="${SCRIPT_DIR}/DTS"

# Default arguments if none provided
if [ $# -eq 0 ]; then
    ARGS="-in input.dts -top topology.top"
else
    ARGS="$@"
fi

# Extract restart file name from arguments or use default
RESTART_FILE=""
if [[ "$ARGS" == *"-restart"* ]]; then
    # Extract restart file from arguments
    RESTART_FILE=$(echo "$ARGS" | grep -oP '(?<=-restart\s)\S+' || echo "")
fi

# If no restart file specified in arguments, look for the most recent one
# This allows resuming from existing checkpoints automatically
if [ -z "$RESTART_FILE" ]; then
    # Look for restart files (typically .res extension or in current directory)
    RESTART_FILE=$(ls -t *.res 2>/dev/null | head -1)
    if [ -n "$RESTART_FILE" ]; then
        echo "Found existing restart file: $RESTART_FILE"
        echo "Will resume from checkpoint (use -restart <file> to override)"
        # Don't add to ARGS yet - we'll add it dynamically on each run to get latest
    fi
fi

# Change to the directory containing input files (if input.dts is specified)
if [[ "$ARGS" == *"-in"* ]]; then
    INPUT_FILE=$(echo "$ARGS" | grep -oP '(?<=-in\s)\S+' || echo "")
    if [ -n "$INPUT_FILE" ] && [ -f "$INPUT_FILE" ]; then
        cd "$(dirname "$(readlink -f "$INPUT_FILE")")" || cd "$(dirname "$INPUT_FILE")"
        INPUT_FILE=$(basename "$INPUT_FILE")
        ARGS=$(echo "$ARGS" | sed "s|-in [^ ]*|-in $INPUT_FILE|")
    fi
fi

# Track restart attempts
RESTART_COUNT=0

echo "=========================================="
echo "FreeDTS Auto-Restart Wrapper"
echo "=========================================="
echo "Command: $DTS_BINARY $ARGS"
echo "Max restarts: $MAX_RESTARTS"
echo "=========================================="

while [ $RESTART_COUNT -le $MAX_RESTARTS ]; do
    # Always check for latest restart file (use it if available)
    LATEST_RESTART=$(ls -t *.res 2>/dev/null | head -1)
    
    if [ $RESTART_COUNT -gt 0 ]; then
        echo ""
        echo "--- Restart attempt $RESTART_COUNT/$MAX_RESTARTS ---"
        
        # Update restart file to the most recent one
        if [ -n "$LATEST_RESTART" ]; then
            echo "Using restart file: $LATEST_RESTART"
            ARGS=$(echo "$ARGS" | sed "s|-restart [^ ]*||" | sed "s|$|-restart $LATEST_RESTART|")
        fi
        
        echo "Waiting $RESTART_DELAY seconds before restarting..."
        sleep $RESTART_DELAY
    else
        # First run - check if restart file exists
        if [ -n "$LATEST_RESTART" ] && [[ "$ARGS" != *"-restart"* ]]; then
            echo ""
            echo "Found existing restart file: $LATEST_RESTART"
            echo "Resuming from checkpoint (use -restart <file> to override, or delete .res files to start fresh)"
            ARGS="$ARGS -restart $LATEST_RESTART"
        elif [ -z "$LATEST_RESTART" ]; then
            echo ""
            echo "No restart file found. Starting new simulation..."
        fi
    fi
    
    # Run the simulation
    echo ""
    echo "Starting simulation..."
    $DTS_BINARY $ARGS
    EXIT_CODE=$?
    
    # Check exit code
    if [ $EXIT_CODE -eq 0 ]; then
        echo ""
        echo "=========================================="
        echo "Simulation completed successfully!"
        echo "=========================================="
        exit 0
    elif [ $EXIT_CODE -eq 139 ]; then
        # Segmentation fault (SIGSEGV)
        echo ""
        echo "ERROR: Segmentation fault detected (exit code 139)"
        RESTART_COUNT=$((RESTART_COUNT + 1))
        
        if [ $RESTART_COUNT -le $MAX_RESTARTS ]; then
            echo "Attempting automatic restart..."
        else
            echo "Maximum restart attempts ($MAX_RESTARTS) reached. Exiting."
            exit 139
        fi
    elif [ $EXIT_CODE -eq 134 ]; then
        # Abort (SIGABRT)
        echo ""
        echo "ERROR: Abort signal detected (exit code 134)"
        RESTART_COUNT=$((RESTART_COUNT + 1))
        
        if [ $RESTART_COUNT -le $MAX_RESTARTS ]; then
            echo "Attempting automatic restart..."
        else
            echo "Maximum restart attempts ($MAX_RESTARTS) reached. Exiting."
            exit 134
        fi
    else
        # Other error
        echo ""
        echo "ERROR: Simulation exited with code $EXIT_CODE"
        echo "Not restarting (only restarts on segfault/abort)"
        exit $EXIT_CODE
    fi
done

exit 1

