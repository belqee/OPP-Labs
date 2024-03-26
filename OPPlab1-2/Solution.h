#pragma once

#include "Matrix.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>

#define MAIN_THREAD 0

class UsualSolution {
private:
    UsualSolution();

protected:
    static const double ti_plus;
    static const double ti_minus;
    static double ti;
    int N;
    std::vector<double> b;
    std::vector<double> x;
    Matrix A;
    double norm_denominator = 0;

    virtual bool accuracy_check(double epsilon);

    virtual double multiply_v(const std::vector<double>& a, const std::vector<double>& b) const;

    double find_norm(const std::vector<double>& v, int size) const;

    virtual void proximity_function();


public:

    void print_v(const std::vector<double>& v) const;

    virtual ~UsualSolution();

    explicit UsualSolution(int N);

    virtual void run(double epsilon);

    virtual void print_result();
};

class FirstSolution : public UsualSolution {
private:
    int size, rank, block_begin, block_end, block_size;

    void proximity_function() override;

    bool accuracy_check(double epsilon) override;

    std::vector<double> block_result;
    std::vector<double> result;
public:
    explicit FirstSolution(int N);

    ~FirstSolution();

    void print_result() override;

    void run(double epsilon) override;
};

class SecondSolution : public UsualSolution {
private:
    int size, rank, block_begin, block_end, block_size, destination, sender;

    void proximity_function() override;

    bool accuracy_check(double epsilon) override;

    double multiply_v(const std::vector<double>& a, const std::vector<double>& b, int offset, int count) const;

    double find_norm_b();

    std::vector<double> block_result;
    std::vector<double> result;
    std::vector<double> tmp;
public:
    explicit SecondSolution(int N);

    ~SecondSolution();

    void print_result() override;

    void run(double epsilon) override;
};
