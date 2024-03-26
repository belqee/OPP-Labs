#include<mpi.h>
#include<vector>
#include<iostream>
#include<chrono>
#include<cmath>

class Matrix
{
private:
    std::vector< std::vector<double> > matrix;
    int N;
public:
    Matrix(int N);

    int size() {
        return N;
    }
    const std::vector<double>& operator[](int index) const;
    std::vector<double>& operator[](int index);
};

Matrix::Matrix(int size) : N(size)
{
    matrix.resize(N, std::vector<double>(N, 1.0));
    for (size_t i = 0; i < matrix.size(); i++)
    {
        for (size_t j = 0; j < matrix.size(); j++)
        {
            if (i == j)
            {
                matrix[i][j] = 2.0;
            }
        }
    }
}

std::vector<double>& Matrix::operator[](int index) {
    return matrix[index];
}

const std::vector<double>& Matrix::operator[](int index) const {
    return matrix[index];
}

extern std::chrono::microseconds durationGLOBAL;

class UsualSolution
{
protected:
    double ti = ti_plus;
    static const double ti_plus;
    static const double ti_minus;
    Matrix A;
    int N;
    std::vector<double> x;
    std::vector<double> b;
    double norm_denominator = 0;

    virtual void proximity_function() = 0;
    virtual bool accuracy_check(double epsilon) = 0;
    double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const;
    double find_norm(const std::vector<double>& row) const;
public:
    UsualSolution(Matrix A) : A(A), N(A.size()) {};
    UsualSolution(Matrix A, std::vector<double> x, std::vector<double> b);
    virtual ~UsualSolution() = default;

    virtual void print_result() = 0;
    virtual void execute(double epsilon);
};

class SecondSolution : public UsualSolution
{
private:
    std::vector<double> x_process;
    std::vector<double> result;
    std::vector<double> tmp;

    int size, rank, ibeg, iend, count_for_process, destination, sender;

    int find_norm_b();
    void proximity_function() override;
    bool accuracy_check(double epsilon) override;
    double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column, int offset, int count) const;
public:
    SecondSolution(Matrix A, std::vector<double> x, std::vector<double> b);

    void print_result() override;
};

const double UsualSolution::ti_plus = 0.00001;
const double UsualSolution::ti_minus = -0.00001;

UsualSolution::UsualSolution(Matrix A, std::vector<double> x, std::vector<double> b)
        : A(A), N(A.size()), x(x), b(b)
{
    norm_denominator = find_norm(b);
}

double UsualSolution::multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const
{
    double result = 0.0;
    for (int i = 0; i < N; ++i) {
        result += row[i] * column[i];
    }
    return result;
}

double UsualSolution::find_norm(const std::vector<double>& row) const
{
    double result = 0.0;
    for (size_t i = 0; i < row.size(); i++)
    {
        result += row[i] * row[i];
    }
    result = std::pow(result, 0.5);
    return result;
}

void UsualSolution::execute(double epsilon)
{
    while (accuracy_check(epsilon) == false)
    {
        proximity_function();
    }
}

#define FIRST_THREAD 0

SecondSolution::SecondSolution(Matrix A, std::vector<double> x, std::vector<double> b)
        : UsualSolution(A)
{
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (N % size != 0)
    {
        throw std::runtime_error("Error: N % size is not 0");
    }

    count_for_process = ceil((double)N / size);
    ibeg = count_for_process * rank;
    iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);
    count_for_process = iend - ibeg;

    destination = (rank + 1) % size;
    sender = (rank - 1 + size) % size;

    this->x = x;
    this->b = b;
    x_process.resize(count_for_process, 0.0);
    result.resize(count_for_process, 0.0);
    tmp.resize(count_for_process);
    norm_denominator = find_norm_b();
};

int SecondSolution::find_norm_b()
{
    double sum_part = 0;
    double sum = 0;
    for (int i = 0; i < count_for_process; i++)
    {
        sum_part += b[i] * b[i];
    }
    MPI_Allreduce(&sum_part, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

double SecondSolution::multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column, int offset, int count) const
{
    double result = 0.0;
    for (int i = offset, j = 0; i < offset + count; i++, j++) {
        result += row[i] * column[j];
    }
    return result;
}

void SecondSolution::proximity_function()
{
    x_process.assign(x_process.size(), 0);
    for (int k = 0; k < size; k++)
    {
        int block = (rank + k) % size;

        for (int i = ibeg, j = 0; i < iend; i++, j++)
        {
            x_process[j] += multiply_row_by_column(A[i], x, block * count_for_process, count_for_process);
        }

        MPI_Sendrecv(&x[0], count_for_process, MPI_DOUBLE, destination, 0,
                     &tmp[0], count_for_process, MPI_DOUBLE, sender, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        x = tmp;
        //MPI_Sendrecv_replace(&x[0], count_for_process, MPI_DOUBLE, destination, 0, sender, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (int i = 0; i < count_for_process; i++)
    {
        x[i] = x[i] - ti * (x_process[i] - b[i]);
    }
}

bool SecondSolution::accuracy_check(double epsilon)
{
    result.assign(result.size(), 0);
    for (int k = 0; k < size; k++)
    {
        int block = (rank + k) % size;

        for (int i = ibeg, j = 0; i < iend; i++, j++)
        {
            result[j] += multiply_row_by_column(A[i], x, block * count_for_process, count_for_process);
        }
        MPI_Sendrecv(x.data(), count_for_process, MPI_DOUBLE, destination, 0,
                     tmp.data(), count_for_process, MPI_DOUBLE, sender, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        x = tmp;

    }

    double part_norm = 0;

    for (int i = 0; i < count_for_process; i++)
    {
        part_norm += (result[i] - b[i]) * (result[i] - b[i]);
    }

    double norm_numerator;
    MPI_Allreduce(&part_norm, &norm_numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return norm_numerator / norm_denominator < epsilon * epsilon ? true : false;
}

void SecondSolution::print_result()
{
    int global_x_size = 0;
    int part_x_size = x.size();
    MPI_Allreduce(&part_x_size, &global_x_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::cout << "Count of elements X = " << global_x_size << std::endl;
    }
    for (int i = 0; i < size; i++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i)
        {
            for (size_t i = 0; i < x.size(); i++)
            {
                std::cout << x[i] << ' ';
            }
        }
    }
    if (rank == 0)
    {
        std::cout << std::endl;
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int N = 1000;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count_for_process = ceil((double)N / size);
    int ibeg = count_for_process * rank;
    int iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);
    count_for_process = iend - ibeg;

    Matrix A(N);
    std::vector<double> x(count_for_process, 0);
    std::vector<double> b(count_for_process, N + 1);
    try {
        auto start_time = std::chrono::high_resolution_clock::now();

        SecondSolution solver(A, x, b);
        solver.execute(0.0001);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        if (rank == 0)
        {
            std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;
        }
        solver.print_result();
    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }

    MPI_Finalize();
}
