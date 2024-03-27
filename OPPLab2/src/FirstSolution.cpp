#include "Solution.h"


static int n = 0;

double FirstSolution::multiply_v(const double *a, const double *b) const {
    double result = 0.0;
    #pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

FirstSolution::FirstSolution(int N) : UsualSolution(N) {
}

double FirstSolution::find_norm(const double *v, int size) const {
    double result = 0.0;
#pragma omp parallel for reduction(+:result)
    for (int i = 0; i < size; i++) {
        result += v[i] * v[i];
    }
    result = pow(result, 0.5);
    return result;
}

void FirstSolution::run(double epsilon) {
    while (!accuracy_check(epsilon)) {
        proximity_function();
    }
}

void FirstSolution::proximity_function() {
    double *new_x = new double[N];
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        new_x[i] = multiply_v(A[i], x);
    }
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        x[i] = x[i] - ti * (new_x[i] - b[i]);
    }
    delete[] new_x;
}

bool FirstSolution::accuracy_check(double epsilon) {
    double *result = new double[N];
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        result[i] = multiply_v(A[i], x);
    }
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        result[i] = (result[i] - b[i]);
    }
    double norm_numerator = find_norm(result, N);
    delete[] result;
    return norm_numerator / norm_denominator < epsilon;
}

void FirstSolution::print_result() {
    for (int i = 0; i < N; i++) {
        cout << x[i] << ' ';
    }
    cout << endl;
}


FirstSolution::~FirstSolution() {
    delete[] x;
    delete[] b;
}
