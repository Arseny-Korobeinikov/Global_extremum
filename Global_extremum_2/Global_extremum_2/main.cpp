#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <chrono>
using namespace std;
#define _USE_MATH_DEFINES

double f_R(double m, double a, double b, double f_a, double f_b) {
	return (m * (b - a) + ((f_b - f_a) * (f_b - f_a)) / (m * (b - a)) - 2 * (f_b + f_a));
}

struct Result {
	double x;
	double y;
	int k;
	double m;
	double time = 0.0;
	Result(double x_, double y_, int k_, double m_) {
		x = x_;
		y = y_;
		k = k_;
		m = m_;
	}
	const Result& operator=(const Result& obj) {
		if (this == &obj) return (*this);
		x = obj.x;
		y = obj.y;
		k = obj.k;
		m = obj.m;

		return (*this);
	}
	friend ostream& operator<<(ostream& os, const Result& obj)
	{
		os << "arg_extremum = " << obj.x << "\nobj_func_extremum = " << obj.y
			<< "\ncount operation to search = " << obj.k
			<< "\ntime = " << obj.time
			<< "\nConst L = " << obj.m << endl << endl;
		return os;
	}
};

class Segment {
public:
	double a;
	double b;
	double R;

	Segment() : a(), b() {}
	Segment(double a_, double b_, double R_) {
		a = a_;
		b = b_;
		R = R_;
	}

	bool operator == (const Segment& obj) {
		return (a == obj.a && b == obj.b);
	}

};

//class SegmentForQueue: public Segment{
//public:
//	double R;
//	Segment seg_vec;
//
//	SegmentForQueue(double a_, double b_, double R_, Segment& seg_vec_) {
//		seg_vec = seg_vec_;
//		a = a_;
//		b = b_;
//		R = R_;	
//	}
//	bool operator > (const SegmentForQueue& obj) {
//		return R > obj.R;
//	}
//	bool operator < (const SegmentForQueue& obj) {
//		return R < obj.R;
//	}
//	
//};

struct Comp {
	bool operator () (const Segment* obj1, const Segment* obj2) {
		return obj1->R < obj2->R;
	}
};

Result algorithm_for_searching(double a, double b, double (*f)(double x), double r, int max_count_operation, double eps) {
	int count_operation = 0;
	vector<Segment*>  seg_v;
	priority_queue<Segment*, vector<Segment*>, Comp> seg_q;
	double res_arg;
	if (f(a) > f(b)) {
		res_arg = b;
	}
	else {
		res_arg = a;
	}

	double max_M = abs((f(b) - f(a)) / (b - a)), m = 0;
	if (max_M == 0)
		m = 1;
	else
		m = r * max_M;


	double x_k = 0.5 * (b + a) - (f(b) - f(a)) / (2 * m);
	if (f(res_arg) > f(x_k)) {
		res_arg = x_k;
	}

	double M1 = abs((f(x_k) - f(a)) / (x_k - a));
	double M2 = abs((f(b) - f(x_k)) / (b - x_k));
	if ((M1 * r > m) || (M2 * r > m)) {
		m = r * max(M1, M2);
	}
	double R = f_R(m, a, x_k, f(a), f(x_k));
	Segment* tmp_segment = new Segment(a, x_k, R);
	seg_v.push_back(tmp_segment);


	R = f_R(m, x_k, b, f(x_k), f(b));
	tmp_segment = new Segment(x_k, b, R);
	seg_v.push_back(tmp_segment);


	seg_q.push(seg_v[0]);
	seg_q.push(seg_v[1]);
	while (count_operation < max_count_operation) {
		Segment* tmp = seg_q.top();
		seg_q.pop();
		double a_ = tmp->a, b_ = tmp->b;

		double x_k = 0.5 * (tmp->b + tmp->a) - (f(tmp->b) - f(tmp->a)) / (2 * m);
		if (f(res_arg) > f(x_k)) {
			res_arg = x_k;
		}
		M1 = abs((f(x_k) - f(tmp->a)) / (x_k - tmp->a));
		M2 = abs((f(tmp->b) - f(x_k)) / (tmp->b - x_k));

		R = f_R(m, tmp->a, x_k, f(tmp->a), f(x_k));
		tmp->b = x_k;
		tmp->R = R;

		R = f_R(m, x_k, b_, f(x_k), f(b_));
		tmp_segment = new Segment(x_k, b_, R);
		seg_v.push_back(tmp_segment);

		if (M1 * r > m || M2 * r > m) {
			m = r * max(M1, M2);
			int s = seg_q.size();
			for (int i = 0; i < s; i++) {
				seg_q.pop();
			}

			for (int i = 1; i < seg_v.size(); i++) {
				R = f_R(m, seg_v[i]->a, seg_v[i]->b, f(seg_v[i]->a), f(seg_v[i]->b));
				seg_v[i]->R = R;
				seg_q.push(seg_v[i]);
			}
		}
		else {
			seg_q.push(tmp);
			seg_q.push(tmp_segment);
		}
		count_operation++;
		if (b_ - a_ < eps) {
			break;
		}
	}

	Result res = { res_arg, f(res_arg), count_operation, m };

	for (int i = 0; i < seg_v.size(); i++) {
		delete seg_v[i];
	}
	return res;
}

double f_sin(double x) {
	return -sin(x);
}

double f_x(double x) {
	return (x - 10) * (x - 10) + 3;
}

double f_1(double x) {
	return (sin(x) + sin(10 * x / 3));
}

double f_2(double x) {
	double res = 0;
	for (int k = 1; k < 6; k++) {
		res += k * sin((k + 1) * x + k);
	}
	return -res;
}

double f_3(double x) {
	return ((3 * x - 1.4) * sin(18 * x));
}

double f_4(double x) {
	return -((x + sin(x)) * exp(-(x * x)));
}

double f_5(double x) {
	return sin(x) + sin(10 * x / 3) + log(x) - 0.84 * x + 3;
}

double f_6(double x) {
	//const double PI = acos(-1.0);
	const double PI = 3.141592653589793;
	return -sin(2 * PI * x) * exp(-x);
}

double f_7(double x) {
	return (x * x - 5 * x + 6) / (x * x + 1);
}

double f_8(double x) {
	return (-x + sin(3 * x) - 1);
}

double f_9(double x) {
	return (2 * (x - 3) * (x - 3) + exp(x * x / 2));
}

int main() {
	double r = 1.5, error_of_arg = 0.001, max_count_operation = 100000;
	auto start = chrono::high_resolution_clock::now();
	Result res_f1 = algorithm_for_searching(2.7, 7.5, f_1, r, max_count_operation, error_of_arg);
	auto end = chrono::high_resolution_clock::now();
	double time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f1.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f2 = algorithm_for_searching(0.0, 10.0, f_2, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f2.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f3 = algorithm_for_searching(0.0, 1.2, f_3, 2, max_count_operation, error_of_arg); // not true prog: 0.629343	 true: 0.96609
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f3.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f4 = algorithm_for_searching(-10.0, 10.0, f_4, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f4.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f5 = algorithm_for_searching(2.7, 7.5, f_5, 5, max_count_operation, error_of_arg); // not true prog: 7.06785		true: 5.19978
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f5.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f6 = algorithm_for_searching(0.0, 4.0, f_6, 1.8, max_count_operation, error_of_arg); // not prog: 2.22485		true: 0.222485   !!
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f6.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f7 = algorithm_for_searching(-5.0, 5.0, f_7, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f7.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f8 = algorithm_for_searching(0.0, 6.5, f_8, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f8.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f9 = algorithm_for_searching(-3.0, 3.0, f_9, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f9.time = time_taken;

	cout << "f1:" << endl << res_f1
		<< "f2:" << endl << res_f2
		<< "f3:" << endl << res_f3
		<< "f4:" << endl << res_f4
		<< "f5:" << endl << res_f5
		<< "f6:" << endl << res_f6
		<< "f7:" << endl << res_f7
		<< "f8:" << endl << res_f8
		<< "f9:" << endl << res_f9;

	return 0;
}
