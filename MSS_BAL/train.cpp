#include <iostream>
#include <fstream>
#include <algorithm> 
#include <time.h>
#include <vector>
#include "PiecewiseLinearFunction.h"

using namespace std;

class timeRecord {
public:
	double time;
	int jobId;
	int machineId;
	void set(double timeTemp, int jobTemp, int machineTemp) {
		time = timeTemp;
		jobId = jobTemp;
		machineId = machineTemp;
	}
};

int jobNum = 10;
vector<vector<timeRecord>> timeReq(jobNum, vector<timeRecord>(2));
vector<vector<PiecewiseLinearFunction>> timeReqPara(jobNum, vector<PiecewiseLinearFunction>(2));
PiecewiseLinearFunction predictValue;

double Johnson(int cnt, vector<int> schedule, vector<int> remainJob) {
	double totalTime;
	if (cnt == jobNum) {
		double M1Time = 0;
		double M2Time = 0;
		for (int i = 0; i < jobNum; i++) {
			int currJob = schedule[i];
			M1Time = M1Time + timeReq[currJob][0].time;
			M2Time = M1Time < M2Time ? M2Time : M1Time;
			M2Time = M2Time + timeReq[currJob][1].time;
		}
		return M2Time;
	}
	else {
		vector<double> remainJobTime;
		vector<int> remainJobId;
		vector<int> remainMachineId;
		for (int i = 0; i < remainJob.size(); i++) {
			remainJobTime.push_back(timeReq[remainJob[i]][0].time);
			remainJobTime.push_back(timeReq[remainJob[i]][1].time);
			remainJobId.push_back(timeReq[remainJob[i]][0].jobId);
			remainJobId.push_back(timeReq[remainJob[i]][1].jobId);
			remainMachineId.push_back(timeReq[remainJob[i]][0].machineId);
			remainMachineId.push_back(timeReq[remainJob[i]][1].machineId);
		}
		int minElementIndex = std::min_element(remainJobTime.begin(), remainJobTime.end()) - remainJobTime.begin();
		int insertJobId = remainJobId[minElementIndex];
		int corrspondingMId = remainMachineId[minElementIndex];
		if (corrspondingMId == 0) {
			for (int i = 0; i < jobNum; i++) {
				if (schedule[i] == -1) {
					schedule[i] = insertJobId;
					break;
				}
			}
		}
		else {
			for (int i = jobNum - 1; i >= 0; i--) {
				if (schedule[i] == -1) {
					schedule[i] = insertJobId;
					break;
				}
			}
		}
		remainJob.erase(std::remove(remainJob.begin(), remainJob.end(), insertJobId), remainJob.end());
		cnt = cnt + 1;
		totalTime = Johnson(cnt, schedule, remainJob);
	}
	return totalTime;
}

double JohnsonTrain(int cnt, vector<int> schedule, vector<int> remainJob) {
	double totalTime = 0;
	if (cnt == jobNum) {
		double M1Time = 0;
		double M2Time = 0;
		for (int i = 0; i < jobNum; i++) {
			int currJob = schedule[i];
			M1Time = M1Time + timeReqPara[currJob][0].id[0];
			M2Time = M1Time < M2Time ? M2Time : M1Time;
			M2Time = M2Time + timeReqPara[currJob][1].id[0];
		}

		/*cout << "schedule: " << endl;
		for (int i = 0; i < jobNum; i++) {
			cout << schedule[i] << " ";
		}
		cout << endl;*/

		for (int j = 0; j < jobNum; j++) {
			remainJob.push_back(j);
			schedule[j] = -1;
		}
		totalTime = M2Time;

		predictValue.assign(timeReqPara[0][0].start_x[0], timeReqPara[0][0].end_x[0], 0, 0, totalTime, 0, 0);
		predictValue.merge();
		/*cout << "predictValue" << endl;
		predictValue.output();
		cout << endl;*/

		return totalTime;
	}
	else {
		vector<PiecewiseLinearFunction> remainJobTime;
		for (int i = 0; i < remainJob.size(); i++) {
			remainJobTime.push_back(timeReqPara[remainJob[i]][0]);
			remainJobTime.push_back(timeReqPara[remainJob[i]][1]);
		}
		PiecewiseLinearFunction minRemainJobTime = MinCompare(remainJobTime);
		minRemainJobTime.merge();
		/*cout << "minRemainJobTime" << endl;
		minRemainJobTime.output();
		cout << endl;*/

		for (int i = 0; i < minRemainJobTime.a.size(); i++) {
			for (int j = 0; j < jobNum; j++) {
				timeReqPara[j][0].start_x[0] = minRemainJobTime.start_x[i];
				timeReqPara[j][1].start_x[0] = minRemainJobTime.start_x[i];
				timeReqPara[j][0].end_x[0] = minRemainJobTime.end_x[i];
				timeReqPara[j][1].end_x[0] = minRemainJobTime.end_x[i];
			}
			vector<int> scheduleTemp;
			vector<int> remainJobTemp;
			scheduleTemp = schedule;
			remainJobTemp = remainJob;
			int cntTemp = cnt;

			int insertJobId = minRemainJobTime.jobId[i];
			int corrspondingMId = minRemainJobTime.machineId[i];
			if (corrspondingMId == 0) {
				for (int i = 0; i < jobNum; i++) {
					if (schedule[i] == -1) {
						schedule[i] = insertJobId;
						break;
					}
				}
			}
			else {
				for (int i = jobNum - 1; i >= 0; i--) {
					if (schedule[i] == -1) {
						schedule[i] = insertJobId;
						break;
					}
				}
			}
			remainJob.erase(std::remove(remainJob.begin(), remainJob.end(), insertJobId), remainJob.end());
			cnt = cnt + 1;
			totalTime = JohnsonTrain(cnt, schedule, remainJob);
			schedule = scheduleTemp;
			remainJob = remainJobTemp;
			cnt = cntTemp;
		}

	}
	return totalTime;
}

int main(int argc, char* argv[]) {

  if (argc != 6) {
		std::cerr << "usage: ./train [train data] [test data] [prediction file] [iteration] [runtime file]" << std::endl;
		exit(0);
	}

	int featureNum = 8;
	vector<double> alpha(featureNum, 1);
	vector<double> regretRec(featureNum, 0);
	/*for(int i = 0; i < featureNum; i++){
	  alpha[i] = (rand() % (MAX-MIN+1))+MIN;
	}*/
	//double alpha[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
	//double regretRec[8] = { 0,0,0,0,0,0,0,0 };
	double regret = INF;
	vector<PiecewiseLinearFunction> ans;
	PiecewiseLinearFunction finalRes;
	int trainNum = 70;
	int testNum = 30;

	int k = 0;

	ifstream infile;
	//ifstream flowfile;
	double benchmarkId;

	std::cout << "==================================================== training ====================================================" << endl;
	time_t Tstart, Tend;
	Tstart = time(NULL);

	int iteration = atoi(argv[4]);
	//int iteration = 1;
	for (int T = 0; T < iteration; T++) {
		for (int k = 0; k < featureNum; k++) {

			infile.open(argv[1]);
			//infile.open("train_GEANT_100.txt");

			for (int benN = 0; benN < trainNum; benN++) {

				vector<vector<double>> trainFeature(jobNum * 2, vector<double>(featureNum, 0));
				vector<double> trainTimeReq(jobNum * 2, 0);

				for (int i = 0; i < jobNum * 2; i++) {
					infile >> benchmarkId;
					for (int j = 0; j < featureNum; j++) {
						infile >> trainFeature[i][j];
					}
					infile >> trainTimeReq[i];
				}

				vector<double> A(jobNum * 2, 0);
				vector<double> B(jobNum * 2, 0);
				for (int i = 0; i < jobNum * 2; i++) {
					A[i] = trainFeature[i][k];
					B[i] = 0;
					for (int d = 0; d < featureNum; d++) {
						B[i] = B[i] + trainFeature[i][d] * alpha[d];
					}
					B[i] = B[i] - trainFeature[i][k] * alpha[k];
				}

				int countNum = 0;
				for (int i = 0; i < jobNum; i++) {
					PiecewiseLinearFunction l0, l1;
					l0.assign(-INF, INF, A[countNum], B[countNum], trainTimeReq[countNum], i, 0);
					timeReqPara[i][0] = l0;
					timeReq[i][0].set(trainTimeReq[countNum], i, 0);
					countNum++;
					l1.assign(-INF, INF, A[countNum], B[countNum], trainTimeReq[countNum], i, 1);
					timeReqPara[i][1] = l1;
					timeReq[i][1].set(trainTimeReq[countNum], i, 1);
					countNum++;
					//cost[i].assign(-INF, INF, A[i], B[i], trainCost[i]);
					//cost[i].output();
				}

				/*for (int i = 0; i < jobNum; i++) {
					cout << timeReq[i][0].time << " " << timeReq[i][1].time << endl;
				}*/

				vector<int> scheduleForTrue(jobNum, -1);
				vector<int> remainJobForTrue(jobNum);
				for (int i = 0; i < jobNum; i++) {
					remainJobForTrue[i] = i;
				}
				double realValue = Johnson(0, scheduleForTrue, remainJobForTrue);
				//std::cout << realValue << " ";

				vector<int> scheduleForPre(jobNum, -1);
				vector<int> remainJobForPre(jobNum);
				for (int i = 0; i < jobNum; i++) {
					remainJobForPre[i] = i;
				}
				double preValue = JohnsonTrain(0, scheduleForPre, remainJobForPre);

				predictValue.computeRegret(realValue);
				//predictValue.output();
				//predictValue.merge();
				ans.push_back(predictValue);

				predictValue.start_x.clear();
				predictValue.end_x.clear();
				predictValue.a.clear();
				predictValue.b.clear();
				predictValue.id.clear();
				predictValue.jobId.clear();
				predictValue.machineId.clear();
			}

			PiecewiseLinearFunction res;
			res.assign(-INF, INF, 0, 0, 0, 0, 0);
			int g = 0;
			while (g < ans.size()) {
				res = Plus(res, ans[g]);
				g++;
			}
			//std::cout << endl;
			std::cout << "=============================== the " << k << "th alpha done================================" << endl;
			//res.output();
			res.mergeSimple();

			vector<double>::iterator smallest = min_element(begin(res.id), end(res.id));
			int index = distance(begin(res.id), smallest);
			double Rtemp = *smallest;
			std::cout << "minimum regret = " << Rtemp << " at position " << index << endl;

			if (Rtemp < regret) {
				regret = Rtemp;
				if (index == 0) {
					alpha[k] = res.end_x[index] - 1;
				}
				else if (index == res.start_x.size() - 1) {
					alpha[k] = res.start_x[index] + 1;
				}
				else {
					alpha[k] = (res.start_x[index] + res.end_x[index]) / 2;
				}
				std::cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			else {
				std::cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			//std::cout << "minimum regret = " << regret << endl;
			regretRec[k] = regret;

			if (k == (featureNum - 1)) {
				finalRes = res;
			}

			ans.clear();
			infile.close();
			//flowfile.close();
		}
		if (regretRec[featureNum - 1] == regretRec[0])
			break;
	}
	std::cout << "================================================= training done ==================================================" << endl;
	vector<double>::iterator finalSmallest = min_element(begin(finalRes.id), end(finalRes.id));
	double totalRegret = *finalSmallest / trainNum;
	//cout << "trainRegret = " << totalRegret << endl;
	Tend = time(NULL);
	double diff = difftime(Tend, Tstart);
	std::cout << "Time = " << diff << endl;

	ofstream ofile1;
	ofile1.open(argv[5]);
	//ofile1.open("runtime.txt");
	ofile1 << diff << endl;
	ofile1.close();

	// read test data 
	// int testbenchmarkNum, testfeatureNum; 
	infile.open(argv[2]);
	//infile.open("test_GEANT_100.txt");
	double testFeatures[testNum][jobNum * 2][featureNum];
	double testRealTime[testNum][jobNum * 2];
	ofstream ofile;
	ofile.open(argv[3]);
	//ofile.open("BAL_GEANT.txt");
	for (int i = 0; i < testNum; i++) {
		for (int j = 0; j < jobNum * 2; j++) {
			double benchmarkId;
			infile >> benchmarkId;
			double predictedTime = 0.0;
			for (int k = 0; k < featureNum; k++) {
				infile >> testFeatures[i][j][k];
				predictedTime += testFeatures[i][j][k] * alpha[k];
			}
			// std::cout << std::endl;
			infile >> testRealTime[i][j];
			ofile << i << " " << testRealTime[i][j] << " " << predictedTime << std::endl;
		}
	}
	//ofile << "Time = " << diff << endl;
	ofile.close();
	infile.close();
}