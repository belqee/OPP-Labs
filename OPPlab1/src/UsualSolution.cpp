#include "Solution.h"

const double UsualSolution::ti_plus = 0.00001;
const double UsualSolution::ti_minus = -0.0001;
double UsualSolution::ti = ti_plus;

static int n = 0;

void UsualSolution::print_v(double *v) const {
    for (int i = 0; i < N; ++i) {
        cout << v[i] << ' ';
    }
    cout << endl;
}

double UsualSolution::multiply_v(const double *a, const double *b) const {
    double result = 0.0;
    for (size_t i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

UsualSolution::UsualSolution(int N) : A(N) {
    this->N = N;
    b = new double[N];
    x = new double[N];
    for (size_t i = 0; i < N; ++i) {
        b[i] = N + 1;
    }
    for (size_t i = 0; i < N; ++i) {
        x[i] = 0;
    }
    norm_denominator = find_norm(b, N);
}

double UsualSolution::find_norm(const double *v, int size) const {
    double result = 0.0;
    for (int i = 0; i < size; i++) {
        result += v[i] * v[i];
    }
    result = pow(result, 0.5);
    return result;
}

void UsualSolution::run(double epsilon) {
    while (!accuracy_check(epsilon)) {
        proximity_function();
    }
}

void UsualSolution::proximity_function() {
    double *new_x = new double[N];
    for (int i = 0; i < N; i++) {
        new_x[i] = multiply_v(A[i], x);
    }

    for (int i = 0; i < N; i++) {
        x[i] = x[i] - ti * (new_x[i] - b[i]);
    }
    delete[] new_x;
}

bool UsualSolution::accuracy_check(double epsilon) {
    double *result = new double[N];
    for (int i = 0; i < N; i++) {
        result[i] = multiply_v(A[i], x);
    }

    for (int i = 0; i < N; i++) {
        result[i] = (result[i] - b[i]);
    }
    double norm_numerator = find_norm(result, N);
    delete[] result;
    return norm_numerator / norm_denominator < epsilon;
}

void UsualSolution::print_result() {
    for (int i = 0; i < N; i++) {
        cout << x[i] << ' ';
    }
    cout << endl;
}

UsualSolution::UsualSolution() {
    cout << "I dont exist";
}

UsualSolution::~UsualSolution() {
    delete[] x;
    delete[] b;
}
