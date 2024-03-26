#pragma once

#include <iostream>
#include <vector>

using namespace std;

class Matrix {
private:
    int N;
    vector<vector<double>> matrix;
public:
    explicit Matrix(int N);
    Matrix();
    ~Matrix();

    int get_size() const;

    vector<double>& operator[](size_t index);

    const vector<double>& operator[](size_t index) const;

    void show();
};