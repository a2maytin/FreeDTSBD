#!/bin/bash
# Simple one-liner style auto-restart wrapper
# Usage: ./run_with_auto_restart_simple.sh

MAX_RESTARTS=10
RESTART_COUNT=0

while [ $RESTART_COUNT -le $MAX_RESTARTS ]; do
    # Always check for latest restart file (use it if available)
    LATEST_RESTART=$(ls -t *.res 2>/dev/null | head -1)
    
    if [ $RESTART_COUNT -gt 0 ]; then
        echo "Restart attempt $RESTART_COUNT/$MAX_RESTARTS"
    else
        if [ -n "$LATEST_RESTART" ]; then
            echo "Found existing restart file: $LATEST_RESTART"
            echo "Resuming from checkpoint..."
        else
            echo "No restart file found. Starting new simulation..."
        fi
    fi
    
    # Use restart file if available, otherwise start fresh
    if [ -n "$LATEST_RESTART" ]; then
        ../DTS -in input.dts -top topology.top -restart "$LATEST_RESTART" -nt 30
    else
        ../DTS -in input.dts -top topology.top -nt 30
    fi
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "Simulation completed successfully!"
        exit 0
    elif [ $EXIT_CODE -eq 139 ] || [ $EXIT_CODE -eq 134 ]; then
        # Segfault (139) or Abort (134)
        RESTART_COUNT=$((RESTART_COUNT + 1))
        if [ $RESTART_COUNT -gt $MAX_RESTARTS ]; then
            echo "Max restarts reached. Exiting."
            exit $EXIT_CODE
        fi
        sleep 1
    else
        echo "Exited with code $EXIT_CODE (not a crash). Not restarting."
        exit $EXIT_CODE
    fi
done

