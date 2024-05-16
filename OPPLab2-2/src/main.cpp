#include <iostream>
#include<chrono>
#include "Matrix.h"
#include <cmath>
#include <omp.h>
#include <vector>
#include <chrono>

Matrix::Matrix(int N) : N(N) {
    matrix = new double *[N];

#pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        matrix[i] = new double[N];

        for (size_t j = 0; j < N; ++j) {
            if (i == j) {
                matrix[i][j] = 2.0;
            } else {
                matrix[i][j] = 1.0;
            }
        }
    }
}

Matrix::~Matrix() {
    for (size_t i = 0; i < N; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

int Matrix::get_size() const {
    return N;
}

double *Matrix::operator[](size_t index) {
    if (index < N) return matrix[index];
    printf("Out of matrix size: %zu >= N\n", index, N);
    return nullptr;
}

const double *Matrix::operator[](size_t index) const {
    if (index < N) return matrix[index];
    printf("Out of matrix size: %zu >= N\n", index, N);
    return nullptr;
}

Matrix::Matrix() {
    N = 0;
    matrix = nullptr;
}

void Matrix::show() {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << matrix[i][j];
        }
        cout << endl;
    }
}

int N = 9600;
const double t_plus = 0.00001;
const double epsilon = 0.00001;

double multiply_v(const double *a, const std::vector<double> &b) {
    double result = 0.0;
    for (size_t i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void proximity() {
    Matrix A(N);
    std::vector<double> b(N, N + 1);
    std::vector<double> result(N);
    std::vector<double> x(N, 0.0);
    std::vector<double> new_x(N);
    double b_norm = 0;
    double res = 0.0;


#pragma omp parallel shared(N, b, x, A, new_x, b_norm, res)
    {
        do {
#pragma omp single
            {
                res = 0.0;
                b_norm = 0;
            }

#pragma omp for //schedule(dynamic, atoi(argv[1])) reduction(+:new_x)
            for (int i = 0; i < N; i++) {
                new_x[i] = multiply_v(A[i], x);
            }
#pragma omp for
            for (int i = 0; i < N; i++) {
                new_x[i] = x[i] - t_plus * (new_x[i] - b[i]);
            }
#pragma omp single
            {
                x = new_x;
            }
#pragma omp for
            for (int i = 0; i < N; i++) {
                result[i] = multiply_v(A[i], x) - b[i];
            }
#pragma omp for reduction(+:res)
            for (int i = 0; i < N; i++) {
                res += result[i] * result[i];
            }
#pragma omp for reduction(+:b_norm)
            for (int i = 0; i < b.size(); i++) {
                b_norm += b[i] * b[i];
            }

        } while (res / b_norm >= epsilon * epsilon);
    }
}

int main(int argc, char **argv) {
    omp_set_num_threads(4);
    auto start_time = std::chrono::high_resolution_clock::now();
    proximity();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;
    return 0;
}