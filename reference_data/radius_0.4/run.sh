#!/bin/bash
export OPENMM_CPU_THREADS=1
pids=""
for w in Walker_*; do
  cd ${w}
  pid=$(python3 ../../simulation_files/run.py > screen.out & echo $!)
  echo "Running ${w} with PID ${pid}"
  pids="${pids}${pid} "
  cd ..
done
echo "To kill walkers call ' kill -9 ${pids} '"
