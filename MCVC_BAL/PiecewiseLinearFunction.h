#pragma once
#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stack>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

#define INF 100000
#define EPS 0.0001

struct PiecewiseLinearFunction
{
	vector<double> a;
	vector<double> b;
	vector<double> start_x;
	vector<double> end_x;
	vector<double> id;
	vector<double> pv;
	double value(double x);
	void assign(double x, double y, double z, double k, double id);
	void output();
	void merge();
	void mergeFinal();
	bool equal(const double& a, const double& b);
	//void predict(vector<double> v);
	void computeRegret(double v);
};

void PiecewiseLinearFunction::computeRegret(double v) {
	for (int i = 0; i < id.size(); i++) {
		id[i] = abs(v - id[i]);
	}
}

/*void PiecewiseLinearFunction::predict(vector<double> v) {
	int i = 0;
	while (i < a.size()) {
		int temp = 0;
		double res = 0;
		char* p = strtok(const_cast<char*>(id[i].c_str()), " ");
		while (p != NULL)
		{
			temp = atoi(p);
			p = strtok(NULL, " ");
			if (temp > 0) {
				res = res + v[temp - 1];
			}
		}
		pv.push_back(res);
		i++;
	}
}*/

bool PiecewiseLinearFunction::equal(const double& a, const double& b) {

	if (abs(a - b) <= EPS)
		return true;

	return false;
}

void PiecewiseLinearFunction::merge() {
	int i = 0;
	vector<double>::iterator a_iter = this->a.begin();
	vector<double>::iterator b_iter = this->b.begin();
	vector<double>::iterator start_iter = this->start_x.begin();
	vector<double>::iterator end_iter = this->end_x.begin();
	vector<double>::iterator id_iter = this->id.begin();
	while (a_iter != this->a.end()) {
		double a_num = *a_iter;
		double b_num = *b_iter;
		double id_num = *id_iter;
		a_iter++;
		b_iter++;
		start_iter++;
		id_iter++;
		//end_iter++;
		while (a_iter != this->a.end() && equal(*a_iter, a_num) && equal(*b_iter, b_num) && equal(*id_iter, id_num)) {
			//while (a_iter != this->a.end() && equal(*id_iter, id_num)) {
			a_iter = this->a.erase(a_iter);
			b_iter = this->b.erase(b_iter);
			start_iter = this->start_x.erase(start_iter);
			end_iter = this->end_x.erase(end_iter);
			id_iter = this->id.erase(id_iter);
		}
		end_iter++;
	}
}

void PiecewiseLinearFunction::mergeFinal() {
	int i = 0;
	vector<double>::iterator a_iter = this->a.begin();
	vector<double>::iterator b_iter = this->b.begin();
	vector<double>::iterator start_iter = this->start_x.begin();
	vector<double>::iterator end_iter = this->end_x.begin();
	vector<double>::iterator id_iter = this->id.begin();
	while (a_iter != this->a.end()) {
		double a_num = *a_iter;
		double b_num = *b_iter;
		double id_num = *id_iter;
		a_iter++;
		b_iter++;
		start_iter++;
		id_iter++;
		//end_iter++;
		//while (a_iter != this->a.end() && equal(*a_iter, a_num) && equal(*b_iter, b_num) && equal(*id_iter, id_num)) {
		while (a_iter != this->a.end() && equal(*id_iter, id_num)) {
			a_iter = this->a.erase(a_iter);
			b_iter = this->b.erase(b_iter);
			start_iter = this->start_x.erase(start_iter);
			end_iter = this->end_x.erase(end_iter);
			id_iter = this->id.erase(id_iter);
		}
		end_iter++;
	}
}

void PiecewiseLinearFunction::output() {
	int i = 0;
	while (i != start_x.size()) {
		std::cout << "[" << start_x.at(i) << ", " << end_x.at(i) << "]: ";
		std::cout << " {" << id.at(i) << "} ";
		std::cout << a.at(i) << "x + " << b.at(i) << endl;
		i++;
	}
	//cout << endl;
}

double PiecewiseLinearFunction::value(double x) {
	int i = 0;
	while (i != start_x.size()) {
		if (x >= start_x[i] && x < end_x[i]) {
			return a[i] * x + b[i];
			break;
		}
		i++;
	}
	return -1;
}

void PiecewiseLinearFunction::assign(double x, double y, double z, double k, double l) {
	start_x.push_back(x);
	end_x.push_back(y);
	a.push_back(z);
	b.push_back(k);
	id.push_back(l);
}

PiecewiseLinearFunction DivisionSpecial(PiecewiseLinearFunction fun1, double k) {
	PiecewiseLinearFunction ds;
	double s = fun1.id[0] / k;
	ds.assign(fun1.start_x[0], fun1.end_x[0], fun1.a[0] / k, fun1.b[0] / k, s);
	return ds;
}

PiecewiseLinearFunction MultipleSpecial(PiecewiseLinearFunction fun1, double k) {
	PiecewiseLinearFunction ds;
	double s = fun1.id[0] * k;
	ds.assign(fun1.start_x[0], fun1.end_x[0], fun1.a[0] * k, fun1.b[0] * k, s);
	return ds;
}

PiecewiseLinearFunction Plus(PiecewiseLinearFunction fun1, PiecewiseLinearFunction fun2) {
	int i = 0;
	int j = 0;
	double s;
	double cur_pos = min(fun1.start_x[i], fun2.start_x[j]);
	PiecewiseLinearFunction comb_fun;
	while (i < fun1.start_x.size() || j < fun2.start_x.size()) {
		if (cur_pos >= fun1.start_x[i] && cur_pos < fun1.end_x[i] && cur_pos >= fun2.start_x[j] && cur_pos < fun2.end_x[j]) {
			comb_fun.a.push_back(fun1.a[i] + fun2.a[j]);
			comb_fun.b.push_back(fun1.b[i] + fun2.b[j]);
			comb_fun.start_x.push_back(cur_pos);
			comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
			s = fun1.id[i] + fun2.id[j];
			comb_fun.id.push_back(s);
		}
		else {
			if (cur_pos >= fun1.start_x[i] && cur_pos < fun1.end_x[i]) {
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(fun1.end_x[i]);
				s = fun1.id[i];
				comb_fun.id.push_back(s);
			}
			if (cur_pos >= fun2.start_x[j] && cur_pos < fun2.end_x[j]) {
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(fun2.end_x[j]);
				s = fun2.id[j];
				comb_fun.id.push_back(s);
			}
		}


		if (fun1.pv.size() && fun2.pv.size()) {
			comb_fun.pv.push_back(fun1.pv[i] + fun2.pv[j]);
		}
		else if (!fun1.pv.size() && fun2.pv.size()) {
			comb_fun.pv.push_back(fun2.pv[j]);
		}
		cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
		if (cur_pos >= fun1.end_x[i])
			i++;
		if (cur_pos >= fun2.end_x[j])
			j++;
	}
	return comb_fun;
}

PiecewiseLinearFunction PlusSpecial(PiecewiseLinearFunction fun1, PiecewiseLinearFunction fun2, int p) {
	PiecewiseLinearFunction ps;
	double s = fun1.id[0] + fun2.id[p];
	ps.assign(fun2.start_x[p], fun2.end_x[p], fun1.a[0] + fun2.a[p], fun1.b[0] + fun2.b[p], s);
	return ps;
}

PiecewiseLinearFunction Minus(PiecewiseLinearFunction fun1, PiecewiseLinearFunction fun2) {
	int i = 0;
	int j = 0;
	double s;
	double cur_pos = min(fun1.start_x[i], fun2.start_x[j]);
	PiecewiseLinearFunction comb_fun;
	while (i < fun1.start_x.size() || j < fun2.start_x.size()) {
		if (cur_pos >= fun1.start_x[i] && cur_pos < fun1.end_x[i] && cur_pos >= fun2.start_x[j] && cur_pos < fun2.end_x[j]) {
			comb_fun.a.push_back(fun1.a[i] - fun2.a[j]);
			comb_fun.b.push_back(fun1.b[i] - fun2.b[j]);
			comb_fun.start_x.push_back(cur_pos);
			comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
			s = fun1.id[i] - fun2.id[j];
			comb_fun.id.push_back(s);
		}
		else {
			if (cur_pos >= fun1.start_x[i] && cur_pos < fun1.end_x[i]) {
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(fun1.end_x[i]);
				s = fun1.id[i];
				comb_fun.id.push_back(s);
			}
			if (cur_pos >= fun2.start_x[j] && cur_pos < fun2.end_x[j]) {
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(fun2.end_x[j]);
				s = fun2.id[j];
				comb_fun.id.push_back(s);
			}
		}


		if (fun1.pv.size() && fun2.pv.size()) {
			comb_fun.pv.push_back(fun1.pv[i] - fun2.pv[j]);
		}
		else if (!fun1.pv.size() && fun2.pv.size()) {
			comb_fun.pv.push_back(fun2.pv[j]);
		}
		cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
		if (cur_pos >= fun1.end_x[i])
			i++;
		if (cur_pos >= fun2.end_x[j])
			j++;
	}
	return comb_fun;
}

PiecewiseLinearFunction MinusSpecial(PiecewiseLinearFunction fun1, PiecewiseLinearFunction fun2, int p) {
	PiecewiseLinearFunction ms;
	double s = fun1.id[0] - fun2.id[p];
	ms.assign(fun2.start_x[p], fun2.end_x[p], fun1.a[0] - fun2.a[p], fun1.b[0] - fun2.b[p], s);
	return ms;
}

PiecewiseLinearFunction SupCombination(PiecewiseLinearFunction fun1, PiecewiseLinearFunction fun2) {
	int i = 0;
	int j = 0;
	double cur_pos = min(fun1.start_x[i], fun2.start_x[j]);
	PiecewiseLinearFunction comb_fun;
	while (i < fun1.start_x.size() || j < fun2.start_x.size()) {
		double x = (fun2.b[j] - fun1.b[i]) / (fun1.a[i] - fun2.a[j]);
		double comp_d = max(fun1.start_x[i], fun2.start_x[j]);
		if (x >= fun1.start_x[i] && x < fun1.end_x[i] && x >= fun2.start_x[j] && x < fun2.end_x[j]) {
			if (fun1.value(comp_d) > fun2.value(comp_d) + EPS) {
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(x);
				comb_fun.id.push_back(fun1.id[i]);
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.start_x.push_back(x);
				comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
				comb_fun.id.push_back(fun2.id[j]);
				cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
			}
			else if (fun1.value(comp_d) + EPS < fun2.value(comp_d)) {
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(x);
				comb_fun.id.push_back(fun2.id[j]);
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.start_x.push_back(x);
				comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
				comb_fun.id.push_back(fun1.id[i]);
				cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
			}
			else if (abs(fun1.value(comp_d) - fun2.value(comp_d)) < EPS) {
				double comp_e = comp_d + EPS;
				if (fun1.value(comp_e) > fun2.value(comp_e)) {
					comb_fun.a.push_back(fun1.a[i]);
					comb_fun.b.push_back(fun1.b[i]);
					comb_fun.start_x.push_back(cur_pos);
					comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
					cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
					comb_fun.id.push_back(fun1.id[i]);
				}
				if (fun1.value(comp_e) < fun2.value(comp_e)) {
					comb_fun.a.push_back(fun2.a[j]);
					comb_fun.b.push_back(fun2.b[j]);
					comb_fun.start_x.push_back(cur_pos);
					comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
					cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
					comb_fun.id.push_back(fun2.id[j]);
				}
			}
		}
		else {
			if (fun1.value(comp_d) >= fun2.value(comp_d)) {
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.id.push_back(fun1.id[i]);
			}
			if (fun1.value(comp_d) < fun2.value(comp_d)) {
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.id.push_back(fun2.id[j]);
			}

			comb_fun.start_x.push_back(cur_pos);
			comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
			cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
		}
		if (cur_pos >= fun1.end_x[i])
			i++;
		//cout << "i is " << i << endl;
		if (cur_pos >= fun2.end_x[j])
			j++;
		//cout << comb_fun.a.size() << "  ";
		//cout << "j is " << j << endl;
		//cout << "fun 1 size is "<<fun1.a.size()<<endl;
		//cout << "fun 2 size is " << fun2.a.size() << endl;
		//bool ans1 = i < fun1.a.size();
		//bool ans2 = j < fun2.a.size();
		//cout << "ans 1 is "<<ans1<<endl;
		//cout << "ans 2 is "<<ans2 << endl;
	}
	return comb_fun;
}

PiecewiseLinearFunction InfCombination(PiecewiseLinearFunction fun1, PiecewiseLinearFunction fun2) {
	int i = 0;
	int j = 0;
	double cur_pos = min(fun1.start_x[i], fun2.start_x[j]);
	PiecewiseLinearFunction comb_fun;
	while (i < fun1.start_x.size() || j < fun2.start_x.size()) {
		double x = (fun2.b[j] - fun1.b[i]) / (fun1.a[i] - fun2.a[j]);
		double comp_d = max(fun1.start_x[i], fun2.start_x[j]);
		if (x >= fun1.start_x[i] && x < fun1.end_x[i] && x >= fun2.start_x[j] && x < fun2.end_x[j]) {
			if (fun1.value(comp_d) + EPS < fun2.value(comp_d)) {
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(x);
				comb_fun.id.push_back(fun1.id[i]);
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.start_x.push_back(x);
				comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
				comb_fun.id.push_back(fun2.id[j]);
				cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
			}
			else if (fun1.value(comp_d) > fun2.value(comp_d) + EPS) {
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.start_x.push_back(cur_pos);
				comb_fun.end_x.push_back(x);
				comb_fun.id.push_back(fun2.id[j]);
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.start_x.push_back(x);
				comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
				comb_fun.id.push_back(fun1.id[i]);
				cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
			}
			else if (abs(fun1.value(comp_d) - fun2.value(comp_d)) < EPS) {
				double comp_e = comp_d + EPS;
				if (fun1.value(comp_e) <= fun2.value(comp_e)) {
					comb_fun.a.push_back(fun1.a[i]);
					comb_fun.b.push_back(fun1.b[i]);
					comb_fun.start_x.push_back(cur_pos);
					comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
					cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
					comb_fun.id.push_back(fun1.id[i]);
				}
				if (fun1.value(comp_e) > fun2.value(comp_e)) {
					comb_fun.a.push_back(fun2.a[j]);
					comb_fun.b.push_back(fun2.b[j]);
					comb_fun.start_x.push_back(cur_pos);
					comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
					cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
					comb_fun.id.push_back(fun2.id[j]);
				}
			}
		}
		else {
			if (fun1.value(comp_d) <= fun2.value(comp_d)) {
				comb_fun.a.push_back(fun1.a[i]);
				comb_fun.b.push_back(fun1.b[i]);
				comb_fun.id.push_back(fun1.id[i]);
			}
			if (fun1.value(comp_d) > fun2.value(comp_d)) {
				comb_fun.a.push_back(fun2.a[j]);
				comb_fun.b.push_back(fun2.b[j]);
				comb_fun.id.push_back(fun2.id[j]);
			}

			comb_fun.start_x.push_back(cur_pos);
			comb_fun.end_x.push_back(min(fun1.end_x[i], fun2.end_x[j]));
			cur_pos = min(fun1.end_x[i], fun2.end_x[j]);
		}
		if (cur_pos >= fun1.end_x[i])
			i++;
		if (cur_pos >= fun2.end_x[j])
			j++;
	}
	return comb_fun;
}

PiecewiseLinearFunction MinCompare(vector<PiecewiseLinearFunction> p) {
	while (p[0].start_x.size() == 0) {
		p.erase(p.begin());
	}
	PiecewiseLinearFunction minCap = p[0];
	for (int i = 1; i < p.size(); i++) {
		if (p[i].start_x.size() != 0) {
			minCap = InfCombination(minCap, p[i]);
		}
	}

	return minCap;
}

PiecewiseLinearFunction MaxCompare(vector<PiecewiseLinearFunction> p) {
	while (p[0].start_x.size() == 0) {
		p.erase(p.begin());
	}
	PiecewiseLinearFunction maxCap = p[0];
	for (int i = 1; i < p.size(); i++) {
		if (p[i].start_x.size() != 0) {
			maxCap = SupCombination(maxCap, p[i]);
		}
	}

	return maxCap;
}
/*
int main() {
	PiecewiseLinearFunction l1, l2, l3;
	l1.assign(-INF, INF, 1, 4, "1 ");
	l2.assign(-INF, INF, 3, -2, "2 ");
	l3 = InfCombination(l1, l2);
	l3.output();
	vector<PiecewiseLinearFunction> cap;
	PiecewiseLinearFunction l1, l2, l3;
	l1.assign(-INF, INF, 1, 2, "1 ");
	l2.assign(-INF, INF, 2, -1, "2 ");
	l3.assign(-INF, INF, 4, -3, "3 ");
	cap.push_back(l1);
	cap.push_back(l2);
	cap.push_back(l3);
	PiecewiseLinearFunction l4 = MinCompare(cap);
	l4.output();
	PiecewiseLinearFunction l4 = Minus(l1, l2);
	l4.output();
}*/
