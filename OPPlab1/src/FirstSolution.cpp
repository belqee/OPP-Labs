#include "Solution.h"

void FirstSolution::proximity_function() {
    double *new_x = new double[N];
    for (int i = block_begin; i <= block_end; i++) {
        new_x[i] = x[i] - ti * (result[i]);
    }
    MPI_Gather(&new_x[block_begin], block_size, MPI_DOUBLE, x, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    delete[] new_x;
}


bool FirstSolution::accuracy_check(double epsilon) {
    for (int i = block_begin, j = 0; j < block_size; i++, j++) {
        block_result[j] = multiply_v(A[i], x) - b[i];
    }
    MPI_Gather(block_result, block_size, MPI_DOUBLE, result, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(result, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double norm_numerator = find_norm(result, N);
    return norm_numerator / norm_denominator < epsilon;
}

FirstSolution::FirstSolution(int N) : UsualSolution(N) {

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    block_size = (N + size - 1) / size;

    block_begin = rank * block_size;
    block_end = size - 1 == rank ? N - 1 : block_size * (rank + 1) - 1;
    block_size = block_end - block_begin + 1;
    block_result = new double[block_size];
    result = new double[N];
    cout << "i'm : " << rank << " my range is: " << block_begin << "-" << block_end << " : " << block_size
         << endl;
}

FirstSolution::~FirstSolution() {
    MPI_Barrier(MPI_COMM_WORLD);
}

void FirstSolution::print_result() {
    if (rank == 0) {
        for (int i = 0; i < 1; i++) {
            cout << x[i] << ' ';
        }
        cout << endl;
    }
}

void FirstSolution::run(double epsilon) {
    int count = 0;
    while (!accuracy_check(epsilon)) {
        proximity_function();
        count++;
    }
    if (rank == 0) cout << count << " cycles" << endl;
    delete[] result;
    delete[] block_result;
    result = nullptr;
    block_result = nullptr;
}
