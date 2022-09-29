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
vector<vector<timeRecord>> trueTimeReq(jobNum, vector<timeRecord>(2));
vector<vector<timeRecord>> preTimeReq(jobNum, vector<timeRecord>(2));
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

double JohnsonTest(int cnt, vector<int> schedule, vector<int> remainJob) {
	double totalTime;
	if (cnt == jobNum) {
		double M1Time = 0;
		double M2Time = 0;
		for (int i = 0; i < jobNum; i++) {
			int currJob = schedule[i];
			M1Time = M1Time + trueTimeReq[currJob][0].time;
			M2Time = M1Time < M2Time ? M2Time : M1Time;
			M2Time = M2Time + trueTimeReq[currJob][1].time;
		}
		return M2Time;
	}
	else {
		vector<double> remainJobTime;
		vector<int> remainJobId;
		vector<int> remainMachineId;
		for (int i = 0; i < remainJob.size(); i++) {
			remainJobTime.push_back(preTimeReq[remainJob[i]][0].time);
			remainJobTime.push_back(preTimeReq[remainJob[i]][1].time);
			remainJobId.push_back(preTimeReq[remainJob[i]][0].jobId);
			remainJobId.push_back(preTimeReq[remainJob[i]][1].jobId);
			remainMachineId.push_back(preTimeReq[remainJob[i]][0].machineId);
			remainMachineId.push_back(preTimeReq[remainJob[i]][1].machineId);
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
		totalTime = JohnsonTest(cnt, schedule, remainJob);
	}
	return totalTime;
}

int main(int argc, char* argv[]) {
  
  if (argc != 3) {
		std::cerr << "usage: ./test [predicted data] [result file]" << std::endl;
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

	ifstream infile;
	double benchmarkId;

	//infile.open("LR_GEANT.txt");
	infile.open(argv[1]);

	double totalRegret = 0;
	double totalSum = 0;
	for (int benN = 0; benN < testNum; benN++) {

		vector<double> trueTimeTemp(jobNum * 2, 0);
		vector<double> preTimeTemp(jobNum * 2, 0);
		for (int i = 0; i < jobNum * 2; i++) {
			infile >> benchmarkId;
			infile >> trueTimeTemp[i];
			infile >> preTimeTemp[i];
		}
		int countNum = 0;
		for (int i = 0; i < jobNum; i++) {
			timeReq[i][0].set(trueTimeTemp[countNum], i, 0);
			trueTimeReq[i][0].set(trueTimeTemp[countNum], i, 0);
			preTimeReq[i][0].set(preTimeTemp[countNum], i, 0);
			countNum++;
			timeReq[i][1].set(trueTimeTemp[countNum], i, 1);
			trueTimeReq[i][1].set(trueTimeTemp[countNum], i, 1);
			preTimeReq[i][1].set(preTimeTemp[countNum], i, 1);
			countNum++;
		}
		
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
		double preValue = JohnsonTest(0, scheduleForPre, remainJobForPre);
		//std::cout << preValue << endl;

		double regret = abs(preValue - realValue);
		totalRegret = totalRegret + regret;
		totalSum = totalSum + realValue;
	}
	totalRegret = totalRegret / testNum;
	totalSum = totalSum / testNum;
	cout << "Regret: " << totalRegret << " ";
	cout << "TOV: " << totalSum << endl;
	ofstream ofile;
	ofile.open(argv[2]);
	//ofile << ifile << " ";
	ofile << totalRegret << " " << totalSum << endl;
	//ofile << totalRegret << endl;
	ofile.close();
	infile.close();
}