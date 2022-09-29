for i in $(eval echo "{$1..$2}")
do
  nohup ./test "BAL/BAL_jobNum10_100(${i}).txt" "result_${i}.txt" &
done
