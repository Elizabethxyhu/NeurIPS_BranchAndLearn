# NeurIPS_BranchAndLearn
This repository is the official implementation of the paper: Branch & Learn for Recursively and Iteratively Solvable Problems in Predict+Optimize

Download and extract the [datasets](https://mycuhk-my.sharepoint.com/:u:/g/personal/1155136882_link_cuhk_edu_hk/ERuXDuDI4lBDisCxvXPW4RIBgcBwf560ULJifXUkqnNxnw?e=Sd3h5P).

MCFP_BAL:
This is the package for the minimum cost flow problem with B&L. The training process and the test process are divided in two files. Please use the following commands:
1. compile the program:
g++ -std=c++11 train.cpp -o train / g++ -std=c++11 test.cpp -o test
2. conduct training:
./train [graph file] [train data] [test data] [result file] [flow] [iteration]
(for example: ./train data/graph_USA_100.txt data/train_USA_100.txt data/test_USA_100.txt data/BAL_USA_100.txt 20 3)
3. conduct testing:
./test [graph file] [test data] [result file] [flow]
(for example: ./test data/graph_USA_100.txt data/BAL_USA_100.txt result.txt 20)

MCVC_BAL:
This is the package for the minimum cost vertex cover problem with B&L. The training process and the test process are divided in two files.
Please use the following commands:
1. compile the program:
g++ -std=c++11 train.cpp -o train / g++ -std=c++11 test.cpp -o test
2. conduct training:
./train [graph file] [train data] [test data] [prediction file] [iteration] [runtime file]
(for example: ./train data/edge_POLSKA.txt data/train_POLSKA_100.txt data/test_POLSKA_100.txt BAL/BAL_POLSKA_100.txt 1 runtime/runtime.txt)
3. conduct testing:
./test [graph file] [predicted data] [result file]
(for example: ./test data/edge_POLSKA.txt BAL/BAL_POLSKA_100.txt result.txt)

MSS_BAL:
This is the package for the multi-stage scheduling problem with B&L. The training process and the test process are divided in two files.
Please use the following commands:
1. compile the program:
g++ -std=c++11 train.cpp -o train / g++ -std=c++11 test.cpp -o test
2. conduct training:
./train [train data] [test data] [prediction file] [iteration] [runtime file]
(for example: ./train data/train_jobNum10_100.txt data/test_jobNum10_100.txt BAL/BAL_jobNum10_100.txt 1 runtime/runtime.txt)
3. conduct testing:
./test [predicted data] [result file]
(for example: ./test BAL/BAL_jobNum10_100.txt result.txt)
