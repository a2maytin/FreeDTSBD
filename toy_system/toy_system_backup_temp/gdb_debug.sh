#!/bin/bash
# GDB script to debug the crash
cd /home/andrew/Desktop/Projects/FreeDTSBD/toy_system
export OMP_NUM_THREADS=1

gdb -batch -ex "set args -in input.dts -top topology.top" \
    -ex "run" \
    -ex "bt" \
    -ex "info registers" \
    -ex "list" \
    -ex "quit" \
    --args ../DTS 2>&1 | tee gdb_output.log

echo "=== GDB output saved to gdb_output.log ==="
cat gdb_output.log

