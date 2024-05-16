#include <mpi.h>
#include <iostream>
#include <cstring>

struct Matrix {
    double* matrix;
    int columns;
    int lines;
    Matrix() : columns(0), lines(0), matrix(nullptr) {}
    Matrix(const int& lines, const int& columns) : columns(columns),
                                                   lines(lines),
                                                   matrix(new double[lines * columns]()) {}

    ~Matrix() {
        delete[] matrix;
    }

    void operator=(const Matrix& other) {
        columns = other.columns;
        lines = other.lines;
        delete[] matrix;
        matrix = new double[columns * lines];
        memcpy(matrix, other.matrix, columns * lines * sizeof(double));
    }

    void fill() {
        for (int i = 0; i < lines * columns; i++) {
            matrix[i] = i;
        }
    }
};

void send_pieces_A(Matrix& matrix, const int& size, int ranks[], const int& lines, const int& columns, const int& p1, MPI_Comm Vertical_Comm) {
    int rows = lines / p1; // кол-во строк у процесса
    int count = rows * columns;             // кол-во элементов

    double* recv_buf = new double[count]();

    if (ranks[1] == 0) {
        MPI_Scatter(matrix.matrix, count, MPI_DOUBLE, recv_buf, count, MPI_DOUBLE, 0, Vertical_Comm);
        matrix = Matrix(rows, columns);
        memcpy(matrix.matrix, recv_buf, count * sizeof(double));
    }
    delete[] recv_buf;
}


double scalar_mul(double* A, double* B, int len) {
    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += A[i] * B[i];
    }
    return sum;
}

double* mul_minor(Matrix& A, Matrix& B) {
    double* C = new double[A.lines * B.columns]();
    for (int i = 0; i < A.lines; i++) {
        for (int j = 0; j < B.columns; j++) {
            for (int k = 0; k < A.columns; k++) {
                C[j + i * B.columns] += *(A.matrix + i * A.columns + k) * *(B.matrix + j + k * B.columns);
            }
        }
    }
    return C;
}

void send_pieces_B(Matrix& B, const int& lines, const int& columns, const int& size, const int ranks[], const int& p2, MPI_Datatype COLUMN, MPI_Comm Horizontal_Comm) {
    int count = columns / p2;
    double* my_column = new double[lines * count];

    if (ranks[1] == 0) {
        MPI_Scatter(B.matrix, 1, COLUMN, my_column, lines * count, MPI_DOUBLE, 0, Horizontal_Comm);
        B = Matrix(lines, count);
        memcpy(B.matrix, my_column, lines * count * sizeof(double));
    }

    delete[] my_column;
}

void broadcast_A(Matrix& A, const int& lines, const int& columns, const int& p1, MPI_Comm Horizontal_Comm) {
    int count = lines / p1 * columns;
    MPI_Bcast(A.matrix, count, MPI_DOUBLE, 0, Horizontal_Comm);
}

void broadcast_B(Matrix& B, const int& lines, const int& columns, const int& p2, MPI_Comm Vertical_Comm) {
    int count = columns / p2 * lines;
    MPI_Bcast(B.matrix, count, MPI_DOUBLE, 0, Vertical_Comm);
}

double* gather_full_matrix_at_null_process(const double* C, const int& n1, const int& n3, const int& p1, const int& p2, const int ranks[2], MPI_Comm Horizontal_Comm, MPI_Comm Vertical_Comm) {
    double* vertical_result = new double[n1 / p1 * n3]();

    for (int i = 0; i < n1 / p1; i++) { // цикл по строчкам
        MPI_Gather(C + i * (n3 / p2), n3 / p2, MPI_DOUBLE, vertical_result + i * n3, n3 / p2, MPI_DOUBLE, 0, Horizontal_Comm);
    }

    double* horizontal_result = new double[n1 * n3]();

    MPI_Gather(vertical_result, n1 / p1 * n3, MPI_DOUBLE, horizontal_result, n1 / p1 * n3, MPI_DOUBLE, 0, Vertical_Comm);

    delete[] vertical_result;
    return horizontal_result;
}


MPI_Comm create_grid_of_process(const int& size, const int& rank, int& p1, int& p2) {
    int arr_dims[2]{ p1, p2 };
    int arr_period[2]{ 0,0 };
    MPI_Comm My_Comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, arr_dims, arr_period, 0, &My_Comm);
    return My_Comm;
}

MPI_Datatype get_new_datatype(const int& lines, const int& columns, const int& p2) {
    MPI_Datatype COLUMN_NOT_RESIZED;
    MPI_Datatype COLUMN_RESIZED;
    MPI_Type_vector(lines, columns / p2, columns, MPI_DOUBLE, &COLUMN_NOT_RESIZED);
    MPI_Type_commit(&COLUMN_NOT_RESIZED);

    MPI_Type_create_resized(COLUMN_NOT_RESIZED, 0, (columns / p2) * sizeof(double), &COLUMN_RESIZED);
    MPI_Type_commit(&COLUMN_RESIZED);
    return COLUMN_RESIZED;
}

MPI_Comm create_horizontal_comm(MPI_Comm& Old_Comm) {
    MPI_Comm New_Comm;
    int remain_dims[2]{ false, true };
    MPI_Cart_sub(Old_Comm, remain_dims, &New_Comm);
    return New_Comm;
}

MPI_Comm create_vertical_comm(MPI_Comm& Old_Comm) {
    MPI_Comm New_Comm;
    int remain_dims[2]{ true, false };
    MPI_Cart_sub(Old_Comm, remain_dims, &New_Comm);
    return New_Comm;
}



void print_result(const double* result, const int& n1, const int& n3, const int& rank) {
    if (rank == 0) {
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n3; j++) {
                std::cout << result[j + i * n3] << " ";
            }
            std::cout << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int size, rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    int n3 = atoi(argv[3]);
    int p1 = atoi(argv[4]);
    int p2 = atoi(argv[5]);

    Matrix A;
    Matrix B;

    if (rank == 0) {
        A = Matrix(n1, n2);
        A.fill();
        B = Matrix(n2, n3);
        B.fill();
    }

    int my_ranks[2];
    MPI_Comm Grid_Comm = create_grid_of_process(size, rank, p1, p2);

    MPI_Cart_coords(Grid_Comm, rank, 2, my_ranks);

    MPI_Comm Horizontal_Comm = create_horizontal_comm(Grid_Comm);
    MPI_Comm Vertiacal_Comm = create_vertical_comm(Grid_Comm);

    double start = MPI_Wtime();

    send_pieces_A(A, size, my_ranks, n1, n2, p1, Vertiacal_Comm);

    MPI_Datatype COLUMN = get_new_datatype(n2, n3, p2);

    send_pieces_B(B, n2, n3, size, my_ranks, p2, COLUMN, Horizontal_Comm);

    broadcast_A(A, n1, n2, p1, Horizontal_Comm);
    broadcast_B(B, n2, n3, p2, Vertiacal_Comm);

    double* C = mul_minor(A, B);

    double* result = gather_full_matrix_at_null_process(C, n1, n3, p1, p2, my_ranks, Horizontal_Comm, Vertiacal_Comm);
    if (rank == 0)
        std::cout << MPI_Wtime() - start << std::endl;
    delete result;
    delete C;
    MPI_Finalize();
}