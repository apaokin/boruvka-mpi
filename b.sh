make
if [ ! -f rmat-$1 ]
then
  echo gen
  ./gen_RMAT -s $1
fi

if [ ! -f rmat-$1.vinfo ]
then
  echo gen_valid
  ./gen_valid_info -in rmat-$1
fi

# ./gen_valid_info -in rmat-10
echo '----------------------------'
rm rmat-$1.mst
# mpirun -n $2 ./mst_reference_mpi -in rmat-$1   -nIters $3 > res_$1_$2.txt
mpirun -n $2 ./mst_reference_mpi -in rmat-$1   -nIters $3

echo '----------------------------'
./validation -in_graph rmat-$1 -in_result rmat-$1.mst -in_valid rmat-$1.vinfo

# ./validation -in_graph rmat-10 -in_result rmat-10.mst -in_valid rmat-10.vinfo
