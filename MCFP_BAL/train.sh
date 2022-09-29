for i in $(eval echo "{$1..$2}")
do
  nohup ./train "graph/graph_USA_100(${i}).txt" "training/train_USA_100(${i}).txt" "testing/test_USA_100(${i}).txt" "BAL/BAL_USA_100(${i}).txt" 20 3 &
done
