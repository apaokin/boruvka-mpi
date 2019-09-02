#!/bin/bash
# mpisubmit.pl ./gen_RMAT -- -s $1
# 3690
make -f Makefile_gdb
if [ ! -f rmat-$1 ]
then
  ./gen_RMAT -s $1
  # mpisubmit.pl -p 1  ./mst_reference_mpi -- -s $1

fi
echo '-------'
#echo mpirun -n $2 -print-all-exitcodes ./mst_reference_mpi -in rmat-$1   -nIters $3
#mpirun -n $2 -print-all-exitcodes ./mst_reference_mpi -in rmat-$1   -nIters $3

mpiexec -n $2  xterm  -e  gdb   --args   ./mst_reference_mpi -in rmat-$1 -nIters $3

# mpisubmit.pl -p $2 --stdout $1-$2-$3.txt ./mst_reference_mpi -- -in rmat-$1   -nIters $3

# mpisubmit.pl -p $2 --stdout $1-$2-$3.txt ./mst_reference_mpi -- -in rmat-$1   -nIters $3
