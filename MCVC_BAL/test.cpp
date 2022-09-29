// MinWeightVertexCover.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <time.h>
#include "PiecewiseLinearFunction.h"

using namespace std;

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
	if (argc != 4) {
		std::cerr << "usage: ./test [graph file] [predicted data] [result file]" << std::endl;
		exit(0);
	}

	int nodeNum = 12, edgeNum = 18;
	int trainNum = 70;
	int testNum = 30;

	int k = 0;

	ifstream infile;
	double benchmarkId;

	infile.open(argv[1]);
	vector<vector<int>> edgeTemp(edgeNum, vector<int>(2, 0));
	for (int i = 0; i < edgeNum; i++) {
		infile >> edgeTemp[i][0];
		infile >> edgeTemp[i][1];
	}
	infile.close();


	// read test data 
	infile.open(argv[2]);
  string ifile = argv[3];

	double totalRegret = 0;
  double totalSum = 0;
	for (int benN = 0; benN < testNum; benN++) {

		vector<double> predictedCost(nodeNum);
		vector<double> realCost(nodeNum);

		for (int j = 0; j < nodeNum; j++) {
			infile >> benchmarkId;
			infile >> realCost[j];
			infile >> predictedCost[j];
		}

		double preValue;

		vector<PiecewiseLinearFunction> cost;
		cost.resize(nodeNum);
		for (int i = 0; i < nodeNum; i++) {
			cost[i].assign(-INF, INF, 0, 0, 0);
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

		MWVC GTrue(edgeNum, nodeNum, realCost, cost, edge);
		vector<int> S;
		int HTC = edgeNum;
		Result realValue = De_MWVC(GTrue, -1, S, HTC);
		//cout << "realCost: " << realValue.totalWeight << endl;

		MWVC GPre(edgeNum, nodeNum, predictedCost, cost, edge);
		S.clear();
		Result preSelected = De_MWVC(GPre, -1, S, HTC);
		preValue = CostTest_MWVC(preSelected.finalSelected, realCost);
		//cout << "preCost: " << preValue << endl;

		double regret = abs(preValue - realValue.totalWeight);
		totalRegret = totalRegret + regret;
    totalSum = totalSum + realValue.totalWeight;
	}
	totalRegret = totalRegret / testNum;
  totalSum = totalSum / testNum;
	cout << "avgRegret: " << totalRegret << " ";
  cout << "avgTOV: " << totalSum << endl;
	ofstream ofile;
	ofile.open(argv[3]);
  //ofile << ifile << " ";
	ofile << totalRegret << " " << totalSum << endl;
  //ofile << totalRegret << endl;
	ofile.close();
  infile.close();

}

