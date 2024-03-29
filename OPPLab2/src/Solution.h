#pragma once

#include "Matrix.h"
#include <iostream>
#include <cmath>
#include <omp.h>
#include <vector>
#include <chrono>

#define MAIN_THREAD 0

class UsualSolution {
private:
    UsualSolution();

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

    virtual double multiply_v(const double *a, const double *b) const;

    virtual double find_norm(const double *v, int size) const;

    virtual void proximity_function();


public:

    void print_v(double *v) const;

    virtual ~UsualSolution();

    explicit UsualSolution(int N);

    virtual void run(double epsilon);

    virtual void print_result();
};

class FirstSolution : virtual public UsualSolution {
private:
    int size, rank, block_begin, block_end, block_size;

    void proximity_function() override;

    bool accuracy_check(double epsilon) override;

    double *block_result = nullptr;
    double *result = nullptr;

    double multiply_v(const double *a, const double *b) const override;

    double find_norm(const double *v, int size) const override;

public:
    FirstSolution(int N);

    ~FirstSolution();

    void print_result() override;

    void run(double epsilon) override;

};

