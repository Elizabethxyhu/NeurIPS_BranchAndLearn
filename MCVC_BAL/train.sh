for i in $(eval echo "{$1..$2}")
do
  nohup ./train "data/edge_POLSKA.txt" "data/train/train_POLSKA_100(${i}).txt" "data/test/test_POLSKA_100(${i}).txt" "BAL/BAL_POLSKA_100(${i}).txt" 1 "runtime/runtime_${i}.txt" &
done
