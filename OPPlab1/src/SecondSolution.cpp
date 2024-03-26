
#include "Solution.h"

#define FIRST_THREAD 0

SecondSolution::SecondSolution(int N) : UsualSolution(N) {
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size == 1) {
        MPI_Finalize();
        throw runtime_error("Error: size is 1");
    }
    block_size = (N + size - 1) / size;
    block_begin = block_size * rank;
    block_end = block_size * (rank + 1) > N ? N : block_size * (rank + 1);
    block_size = block_end - block_begin;
    destination = (rank + 1) % size;
    sender = (rank - 1 + size) % size;
    delete[] x;
    delete[] b;
    x = new double[block_size];
    b = new double[block_size];
    block_result = new double[block_size];
    tmp = new double[block_size];
    result = new double[block_size];
    for (int i = block_begin, j = 0; i < block_end; i++, j++) {
        this->x[j] = x[i];
        this->b[j] = b[i];
    }
    norm_denominator = find_norm_b();
};

SecondSolution::~SecondSolution() {
    MPI_Barrier(MPI_COMM_WORLD);
}

double SecondSolution::find_norm_b() {
    double sum_part = 0;
    double sum = 0;
    for (int i = 0; i < block_size; i++) {
        sum_part += b[i] * b[i];
    }
    MPI_Allreduce(&sum_part, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

double SecondSolution::multiply_v(const double *a, const double *b, int offset, int count) const {
    double sum = 0.0;
    for (int i = offset, j = 0; i < offset + count; i++, j++) {
        sum += a[i] * b[j];
    }
    return sum;
}

void SecondSolution::run(double epsilon) {
    int count = 0;
    while (!accuracy_check(epsilon)) {
        proximity_function();
        count++;
    }
    if (rank == 0) cout << count << " cycles" << endl;
    delete[] result;
    delete[] block_result;
    delete[] tmp;
    result = nullptr;
    block_result = nullptr;
    tmp = nullptr;
}

void SecondSolution::proximity_function() {
    for (int i = 0; i < block_size; ++i) {
        block_result[i] = 0;
    }
    for (int k = 0; k < size; k++) {
        int block = (rank + k) % size;
        for (int i = block_begin, j = 0; i < block_end; i++, j++) {
            block_result[j] += multiply_v(A[i], x, block * block_size, block_size);
        }
        MPI_Sendrecv(&x[0], block_size, MPI_DOUBLE, destination, 0,
                     &tmp[0], block_size, MPI_DOUBLE, sender, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < block_size; ++i) {
            x[i] = tmp[i];
        }
    }
    for (int i = 0; i < block_size; i++) {
        x[i] = x[i] - ti * (block_result[i] - b[i]);
    }
}

bool SecondSolution::accuracy_check(double epsilon) {
    for (int i = 0; i < block_size; ++i) {
        result[i] = 0;
    }
    for (int k = 0; k < size; k++) {
        int block = (rank + k) % size;
        for (int i = block_begin, j = 0; i < block_end; i++, j++) {
            result[j] += multiply_v(A[i], x, block * block_size, block_size);
        }
        MPI_Sendrecv(x, block_size, MPI_DOUBLE, destination, 0,
                     tmp, block_size, MPI_DOUBLE, sender, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < block_size; ++i) {
            x[i] = tmp[i];
        }
    }
    double part_norm = 0;
    for (int i = 0; i < block_size; i++) {
        part_norm += (result[i] - b[i]) * (result[i] - b[i]);
    }
    double norm_numerator;
    MPI_Allreduce(&part_norm, &norm_numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return norm_numerator / norm_denominator < epsilon * epsilon;
}

void SecondSolution::print_result() {
    int global_x_size = 0;
    int part_x_size = block_size;
    MPI_Allreduce(&part_x_size, &global_x_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        cout << "Count of elements X = " << global_x_size << endl;
    }
    for (int i = 0; i < size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i) {
            for (int i = 0; i < block_size; i++) {
                cout << x[i] << ' ';
            }
        }
    }
    if (rank == 0) {
        cout << endl;
    }
}
