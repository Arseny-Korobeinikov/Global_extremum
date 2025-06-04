#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <chrono>
#include <omp.h>
using namespace std;

#define _USE_MATH_DEFINES

double f_R(double m, double a, double b, double f_a, double f_b) {
    return (m * (b - a) + ((f_b - f_a) * (f_b - f_a)) / (m * (b - a)) - 2 * (f_b + f_a));
}

struct Result {
    double x, y, m, time;
    int k;
    Result() = default;
    Result(double x_, double y_, int k_, double m_, double time_)
        : x(x_), y(y_), k(k_), m(m_), time(time_) {
    }
    Result(double x_, double y_, int k_, double m_) {
        x = x_;
        y = y_;
        k = k_;
        m = m_;
    }
    friend ostream& operator<<(ostream& os, const Result& r) {
        os << "arg_extremum = " << r.x
            << "\nobj_func_extremum = " << r.y
            << "\ncount operations = " << r.k
            << "\nConst L = " << r.m
            << "\ntime = " << r.time << "s\n\n";
        return os;
    }
};

struct Segment {
    double a, b, R;
    Segment(double a_, double b_, double R_) : a(a_), b(b_), R(R_) {}
};

struct Comp {
    bool operator()(const Segment* s1, const Segment* s2) const {
        return s1->R < s2->R;
    }
};

Result algorithm_for_searching(double a, double b, double (*f)(double), double r, int max_ops, double eps) {

    vector<Segment*> segs;
    priority_queue<Segment*, vector<Segment*>, Comp> pq;
    double bestX = (f(a) < f(b) ? a : b);
    double M = fabs((f(b) - f(a)) / (b - a));
    double m = (M == 0 ? 1 : r * M);

    {
        double fa = f(a), fb = f(b);
        double xk = 0.5 * (a + b) - (fb - fa) / (2 * m);
        if (f(xk) < f(bestX)) bestX = xk;

        double R1 = f_R(m, a, xk, fa, f(xk));
        segs.push_back(new Segment(a, xk, R1));

        double R2 = f_R(m, xk, b, f(xk), fb);
        segs.push_back(new Segment(xk, b, R2));

        pq.push(segs[0]);
        pq.push(segs[1]);
    }

    vector<Segment*> local_maxR;
    vector<double> nextX, nextF;

    int iter = 0;
    bool done = false;

#pragma omp parallel
    {
        int nthr = omp_get_num_threads();
#pragma omp single
        {
            local_maxR.resize(nthr, nullptr);
            nextX.resize(nthr, 0.0);
            nextF.resize(nthr, 0.0);
        }

        while (!done && iter < max_ops) {
#pragma omp single
            {
                double maxM = 0.0;
                for (auto& s : segs) {
                    double localM = fabs((f(s->b) - f(s->a)) / (s->b - s->a));
                    if (localM > maxM) maxM = localM;
                }
                m = (maxM == 0 ? 1 : r * maxM);

                while (!pq.empty()) pq.pop();
                for (auto& s : segs) {
                    s->R = f_R(m, s->a, s->b, f(s->a), f(s->b));
                    pq.push(s);
                }

                int cnt = min((int)segs.size(), nthr);
                for (int i = 0; i < cnt; ++i) {
                    local_maxR[i] = pq.top();
                    pq.pop();
                }

                for (int i = cnt; i < nthr; ++i) {
                    local_maxR[i] = nullptr;
                }

                ++iter;
            }

#pragma omp for schedule(static)
            for (int i = 0; i < nthr; ++i) {
                Segment* seg = local_maxR[i];
                if (!seg) {
                    continue;
                }
                double xk = 0.5 * (seg->a + seg->b) - (f(seg->b) - f(seg->a)) / (2 * m);
                nextX[i] = xk;
                nextF[i] = f(xk);
            }

#pragma omp single
            {
                for (int i = 0; i < nthr; ++i) {
                    Segment* seg = local_maxR[i];
                    if (!seg) continue;

                    if (nextF[i] < f(bestX)) {
                        bestX = nextX[i];
                    }

                    if (seg->b - seg->a < eps) {
                        done = true;
                        break;
                    }
                    segs.push_back(new Segment(nextX[i], seg->b, 0.0));
                    seg->b = nextX[i];
                }
            }
        }
    }


    double bestY = f(bestX);
    Result res(bestX, bestY, iter, m);

    for (auto p : segs) {
        delete p;
    }
    return res;
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
    Result res_f3 = algorithm_for_searching(0.0, 1.2, f_3, 2, max_count_operation, error_of_arg);
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
    Result res_f5 = algorithm_for_searching(2.7, 7.5, f_5, 5, max_count_operation, error_of_arg);
    end = chrono::high_resolution_clock::now();
    time_taken =
        chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    res_f5.time = time_taken;

    start = chrono::high_resolution_clock::now();
    Result res_f6 = algorithm_for_searching(0.0, 4.0, f_6, 1.8, max_count_operation, error_of_arg);
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