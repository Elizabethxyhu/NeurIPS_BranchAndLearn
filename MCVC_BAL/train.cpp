// MinWeightVertexCover.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <time.h>
#include <stdlib.h>
#include "PiecewiseLinearFunction.h"

using namespace std;

#define MIN -10
#define MAX 10

class MWVC {
public:
	int edgeNum, nodeNum;
	vector<double> realCost;
	vector<PiecewiseLinearFunction> cost;
	vector<vector<bool>> edge;

	MWVC(int eN, int nN, vector<double> rC, vector<PiecewiseLinearFunction> c, vector<vector<bool>> e) {
		edgeNum = eN;
		nodeNum = nN;
		realCost = rC;
		cost = c;
		edge = e;
	}

};

class Result {
public:
	double totalWeight;
	vector<int> finalSelected;
};

Result De_MWVC(MWVC G, int cnt, vector<int> selected, int haveToCover) {

	Result totalResult;
	cnt = cnt + 1;
	if (cnt == G.nodeNum || haveToCover == 0) {
		double totalWeight = 0;

		/*for (int i = 0; i < G.nodeNum; i++) {
			for (int j = 0; j < G.nodeNum; j++) {
				if (G.edge[i][j] == 1) {
					totalWeight = INFINITY;
				}
			}
		}*/
		if (haveToCover > 0) {
			totalWeight = INFINITY;
		}

		if (totalWeight == 0) {
			for (int i = 0; i < selected.size(); i++) {
				totalWeight = totalWeight + G.realCost[selected[i]];
			}
		}

		totalResult.totalWeight = totalWeight;
		totalResult.finalSelected = selected;

		/*for (int i = 0; i < selected.size(); i++) {
			cout << selected[i] << " ";
		}
		cout << ": " << totalWeight << endl;*/

		return totalResult;
	}

	else {
		vector<int> selected1 = selected;
		selected1.push_back(cnt);
		MWVC G1 = G;
		int haveToCover1 = haveToCover;
		for (int i = 0; i < G1.nodeNum; i++) {
			if (G1.edge[i][cnt] == 1) {
				G1.edge[i][cnt] = 0;
				G1.edge[cnt][i] = 0;
				haveToCover1 = haveToCover1 - 1;
			}
		}

		vector<int> selected2 = selected;
		MWVC G2 = G;
		int haveToCover2 = haveToCover;

		Result result1 = De_MWVC(G1, cnt, selected1, haveToCover1);
		Result result2 = De_MWVC(G2, cnt, selected2, haveToCover2);
		Result result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(De_MWVC(G1, cnt, selected1), De_MWVC(G2, cnt, selected2));
		if (result1.totalWeight > result2.totalWeight) {
			result = result2;
		}
		else {
			result = result1;
		}

		return result;
	}
}

PiecewiseLinearFunction CostTrain_MWVC(MWVC G, int cnt, vector<int> selected) {

	PiecewiseLinearFunction totalResult;
	cnt = cnt + 1;
	if (cnt == G.nodeNum) {
		PiecewiseLinearFunction totalWeight;
		totalWeight.assign(-INF, INF, 0, 0, 0);

		for (int i = 0; i < G.nodeNum; i++) {
			for (int j = 0; j < G.nodeNum; j++) {
				if (G.edge[i][j] == 1) {
					totalWeight.a[0] = INF;
					totalWeight.b[0] = INF;
					totalWeight.id[0] = INF;
				}
			}
		}

		if (totalWeight.a[0] == 0) {
			for (int i = 0; i < selected.size(); i++) {
				totalWeight = Plus(totalWeight, G.cost[selected[i]]);
			}
		}

		totalResult = totalWeight;

		/*for (int i = 0; i < selected.size(); i++) {
			cout << selected[i] << " ";
		}
		cout << ": " << totalWeight.id[0] << endl;*/
		//totalWeight.output();

		return totalResult;
	}

	else {
		vector<int> selected1 = selected;
		selected1.push_back(cnt);
		MWVC G1 = G;
		for (int i = 0; i < G1.nodeNum; i++) {
			G1.edge[i][cnt] = 0;
			G1.edge[cnt][i] = 0;
		}

		vector<int> selected2 = selected;
		MWVC G2 = G;

		PiecewiseLinearFunction result1 = CostTrain_MWVC(G1, cnt, selected1);
		PiecewiseLinearFunction result2 = CostTrain_MWVC(G2, cnt, selected2);
		PiecewiseLinearFunction result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(CostTrain_MWVC(G1, cnt, selected1), CostTrain_MWVC(G2, cnt, selected2));
		if (result1.id[0] != INF && result2.id[0] != INF) {
			result = InfCombination(result1, result2);
			result.merge();
			/*cout << "result" << endl;
			result.output();
			cout << endl;*/
		}
		else if (result1.id[0] == INF) {
			result = result2;
		}
		else {
			result = result1;
		}

		return result;
	}
}

double CostTest_MWVC(vector<int> selected, vector<double> realCost) {
	double totalWeight = 0;
	for (int i = 0; i < selected.size(); i++) {
		totalWeight = totalWeight + realCost[selected[i]];
	}
	return totalWeight;
}

int main(int argc, char* argv[])
{
	if (argc != 7) {
		std::cerr << "usage: ./train [graph file] [train data] [test data] [prediction file] [iteration] [runtime file]" << std::endl;
		exit(0);
	}

	int nodeNum = 12, edgeNum = 18;
	int featureNum = 8;
	vector<double> alpha(featureNum, -1);
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

	infile.open(argv[1]);
	vector<vector<int>> edgeTemp(edgeNum, vector<int>(2, 0));
	for (int i = 0; i < edgeNum; i++) {
		infile >> edgeTemp[i][0];
		infile >> edgeTemp[i][1];
	}
	infile.close();

	cout << "==================================================== training ====================================================" << endl;
	time_t Tstart, Tend;
	Tstart = time(NULL);

	int iteration = atoi(argv[5]);
	for (int T = 0; T < iteration; T++) {
		for (int k = 0; k < featureNum; k++) {

			infile.open(argv[2]);
			for (int benN = 0; benN < trainNum; benN++) {

				vector<PiecewiseLinearFunction> cost;
				cost.resize(nodeNum);
				vector<vector<double>> trainFeature(nodeNum, vector<double>(featureNum, 0));
				vector<double> trainCost(nodeNum, 0);

				for (int i = 0; i < nodeNum; i++) {
					infile >> benchmarkId;
					for (int j = 0; j < featureNum; j++) {
						infile >> trainFeature[i][j];
					}
					infile >> trainCost[i];
				}

				vector<double> A(nodeNum, 0);
				vector<double> B(nodeNum, 0);
				for (int i = 0; i < nodeNum; i++) {
					A[i] = trainFeature[i][k];
					B[i] = 0;
					for (int d = 0; d < featureNum; d++) {
						B[i] = B[i] + trainFeature[i][d] * alpha[d];
					}
					B[i] = B[i] - trainFeature[i][k] * alpha[k];
				}

				for (int i = 0; i < nodeNum; i++) {
					cost[i].assign(-INF, INF, A[i], B[i], trainCost[i]);
					//cost[i].output();
				}
        
        /*double flowID;
      	vector<double> flow(edgeNum, 0);
      	for (int i = 0; i < edgeNum; i++) {
          flowfile >> flowID;
      		flowfile >> flow[i];
          //cout << flow[i] << " ";
      	}
      	vector<double>::iterator smallest = min_element(begin(flow), end(flow));
      	int realIndex = distance(begin(flow), smallest);*/
      
      	vector<vector<bool>> edge(nodeNum, vector<bool>(nodeNum, 0));
      	for (int i = 0; i < edgeNum; i++) {
      		//if (i != realIndex) {
      			int src = edgeTemp[i][0];
      			int des = edgeTemp[i][1];
      			edge[src][des] = 1;
      			edge[des][src] = 1;
      		//}
      	}

				MWVC G(edgeNum, nodeNum, trainCost, cost, edge);

				vector<int> S;
        int HTC = edgeNum;
				Result realValue = De_MWVC(G, -1, S, HTC);

				S.clear();
				PiecewiseLinearFunction predictValue = CostTrain_MWVC(G, -1, S);

				//cout << "realWeight: " << realValue.totalWeight << endl;
				//cout << realValue.totalWeight << " ";
				/*cout << "finalSelected: ";
				for (int i = 0; i < realWeight.finalSelected.size(); i++) {
					cout << realWeight.finalSelected[i] << " ";
				}
				cout << endl;
				preWeight.output();
				cout << endl;*/
				predictValue.computeRegret(realValue.totalWeight);
				//predictValue.output();
				predictValue.merge();
				ans.push_back(predictValue);
			}

			PiecewiseLinearFunction res;
			res.assign(-INF, INF, 0, 0, 0);
			int g = 0;
			while (g < ans.size()) {
				res = Plus(res, ans[g]);
				g++;
			}
			//cout << endl;
			cout << "=============================== the " << k << "th alpha done================================" << endl;
			//res.output();
			res.merge();

			vector<double>::iterator smallest = min_element(begin(res.id), end(res.id));
			int index = distance(begin(res.id), smallest);
			double Rtemp = *smallest;
			cout << "minimum regret = " << Rtemp << " at position " << index << endl;

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
				cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			else {
				cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			cout << "minimum regret = " << regret << endl;
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
	cout << "================================================= training done ==================================================" << endl;
	vector<double>::iterator finalSmallest = min_element(begin(finalRes.id), end(finalRes.id));
	double totalRegret = *finalSmallest / trainNum;
	//cout << "trainRegret = " << totalRegret << endl;
	Tend = time(NULL);
	double diff = difftime(Tend, Tstart);
	cout << "Time = " << diff << endl;
  ofstream ofile1;
  ofile1.open(argv[6]);
  ofile1 << diff << endl;
  ofile1.close();


	// read test data 
	// int testbenchmarkNum, testfeatureNum; 
	infile.open(argv[3]);
	double testFeatures[testNum][nodeNum][featureNum];
	double testRealCost[testNum][nodeNum];
	ofstream ofile;
	ofile.open(argv[4]);
	for (int i = 0; i < testNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			double benchmarkId;
			infile >> benchmarkId;
			double predictedCost = 0.0;
			for (int k = 0; k < featureNum; k++) {
				infile >> testFeatures[i][j][k];
				predictedCost += testFeatures[i][j][k] * alpha[k];
			}
			// std::cout << std::endl;
			infile >> testRealCost[i][j];
			ofile << i << " " << testRealCost[i][j] << " " << predictedCost << std::endl;
		}
	}
	//ofile << "Time = " << diff << endl;
	ofile.close();
	infile.close();


}

