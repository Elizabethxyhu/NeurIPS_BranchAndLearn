for i in $(eval echo "{$1..$2}")
do
  nohup ./test "data/edge_POLSKA.txt" "BAL/BAL_POLSKA_100(${i}).txt" "result_${i}.txt" &
done
