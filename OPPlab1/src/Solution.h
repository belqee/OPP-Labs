#pragma once

#include "Matrix.h"
#include <iostream>
#include <cmath>

#define MAIN_THREAD 0

class AllSolution {
private:
    AllSolution();
protected:
    static const double ti_plus;
    static const double ti_minus;
    static double ti;
    int N;
    double *b;
    double *x;
    Matrix A;
    double norm_denominator = 0;

    virtual bool accuracy_check(double epsilon);

    double multiply_v(const double *a, const double *b) const;

    double find_norm(const double *v, int size) const;

    virtual void proximity_function();


public:

    void print_v(double* v) const;

    ~AllSolution();

    explicit AllSolution(int N);


    virtual void run(double epsilon);


    void print_result();
};


