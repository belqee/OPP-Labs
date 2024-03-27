#include "Matrix.h"

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

void Matrix::show(){
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << matrix[i][j];
        }
        cout << endl;
    }
}
