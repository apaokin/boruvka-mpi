make
# ./gen_RMAT -s $1
# ./gen_valid_info -in rmat-10
# ./gen_valid_info -in rmat-$1
echo '----------------------------'
rm rmat-$1.mst

	# ./mst_reference -in rmat-$1   -nIters 10
 ./mst_reference_mpi -in rmat-$1   -nIters 10

echo '----------------------------'
./validation -in_graph rmat-$1 -in_result rmat-$1.mst -in_valid rmat-$1.vinfo

# ./validation -in_graph rmat-10 -in_result rmat-10.mst -in_valid rmat-10.vinfo
