#include <iostream>
#include <vector>
#include <cmath>
#include <queue>
#include <chrono>
#include <algorithm>
#include <functional>

using namespace std;


double f_R(double m, double a, double b, double f_a, double f_b) {
    return (m * (b - a)
        + (f_b - f_a) * (f_b - f_a) / (m * (b - a))
        - 2.0 * (f_b + f_a));
}

struct Segment {
    double a, b, R;
    Segment(double _a, double _b, double _R) : a(_a), b(_b), R(_R) {}
};

struct Comp {
    bool operator()(const Segment* s1, const Segment* s2) {
        return s1->R < s2->R;
    }
};

struct Result {
    double t_min;
    double f_min;
    int iters;
    double m;
};

Result algorithm_for_searching(
    double a, double b,
    function<double(double)> f,
    double r,
    int max_iters,
    double eps)
{
    priority_queue<Segment*, vector<Segment*>, Comp> pq;
    vector<Segment*> all;

    double fa = f(a), fb = f(b);
    double best_t = (fa < fb ? a : b);
    double best_f = min(fa, fb);

    double M = abs((fb - fa) / (b - a));
    double m = (M == 0 ? 1.0 : r * M);

    double x = 0.5 * (a + b) - (fb - fa) / (2 * m);
    double fx = f(x);

    if (fx < best_f) {
        best_f = fx;
        best_t = x;
    }

    pq.push(new Segment(a, x, f_R(m, a, x, fa, fx)));
    pq.push(new Segment(x, b, f_R(m, x, b, fx, fb)));

    int iters = 0;

    while (!pq.empty() && iters < max_iters) {
        Segment* s = pq.top(); pq.pop();

        if (s->b - s->a < eps)
            break;

        double a_ = s->a, b_ = s->b;
        double fa_ = f(a_), fb_ = f(b_);

        double x_ = 0.5 * (a_ + b_) - (fb_ - fa_) / (2 * m);
        double fx_ = f(x_);

        if (fx_ < best_f) {
            best_f = fx_;
            best_t = x_;
        }

        double M1 = abs((fx_ - fa_) / (x_ - a_));
        double M2 = abs((fb_ - fx_) / (b_ - x_));

        if (max(M1, M2) * r > m)
            m = r * max(M1, M2);

        pq.push(new Segment(a_, x_, f_R(m, a_, x_, fa_, fx_)));
        pq.push(new Segment(x_, b_, f_R(m, x_, b_, fx_, fb_)));

        iters++;
    }

    return { best_t, best_f, iters, m };
}

vector<double> linear_map(
    double t,
    int n,
    const vector<double>& a,
    const vector<double>& b,
    int m)
{
    int K = 1 << m;
    int idx = int(t * pow(K, n));

    vector<double> x(n);

    for (int i = 0; i < n; i++) {
        int c = idx % K;
        idx /= K;
        x[i] = a[i] + (b[i] - a[i]) * double(c) / (K - 1);
    }

    return x;
}

double rastrigin(const vector<double>& x) {
    double res = 10.0 * x.size();
    for (double xi : x)
        res += xi * xi - 10.0 * cos(2 * M_PI * xi);
    return res;
}

struct ND_Function {
    int n;
    vector<double> a, b;
    int peano_level;

    double operator()(double t) const {
        vector<double> x = linear_map(t, n, a, b, peano_level);
        return rastrigin(x);
    }
};

int main() {

    ND_Function F;
    F.n = 2;
    F.a = { -5.12, -5.12 };
    F.b = {  5.12,  5.12 };
    F.peano_level = 15; //9

    double r = 1.655;
    int max_iters = 20000;
    double eps = 1e-10;

    auto start = chrono::high_resolution_clock::now();

    Result res = algorithm_for_searching(
        0.0, 1.0,
        [&](double t) { return F(t); },
        r, max_iters, eps
    );

    auto end = chrono::high_resolution_clock::now();
    double time =
        chrono::duration_cast<chrono::nanoseconds>(end - start).count();

    cout << "===== Multidimensional AGP + Peano =====\n";
    cout << "t* = " << res.t_min << "\n";
    cout << "f(t*) = " << res.f_min << "\n";
    cout << "iterations = " << res.iters << "\n";
    cout << "time (sec) = " << time*1e-9 << "\n";

    vector<double> x_min = linear_map(
        res.t_min, F.n, F.a, F.b, F.peano_level);

    cout << "x* = (";
    for (int i = 0; i < x_min.size(); i++) {
        cout << x_min[i];
        if (i + 1 < x_min.size()) cout << ", ";
    }
    cout << ")\n";

    return 0;
}
