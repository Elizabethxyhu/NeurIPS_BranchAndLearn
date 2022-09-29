for i in $(eval echo "{$1..$2}")
do
  nohup ./test "graph/graph_USA_100(${i}).txt" "BAL/BAL_USA_100(${i}).txt" "${i}_result.txt" 20 &
done
