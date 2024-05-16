#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <windows.h>

using namespace std;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < size; i++) {

        MPI_Barrier(MPI_COMM_WORLD);
        if (i == rank) {
            printf("%d\n", i);
            //cout << i << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}