for i in $(eval echo "{$1..$2}")
do
  nohup ./train "training/train_jobNum10_100(${i}).txt" "testing/test_jobNum10_100(${i}).txt" "BAL/BAL_jobNum10_100(${i}).txt" 1 "runtime/runtime_${i}.txt" &
done