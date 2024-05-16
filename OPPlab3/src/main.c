#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int size;

typedef struct {
    int rows;
    int cols;
    double *data;
} Matrix;

typedef struct {
    MPI_Comm commGrib;
    MPI_Comm commCol;
    MPI_Comm commRow;
    MPI_Datatype rowType;
    MPI_Datatype colType;
    int *dims;
} MPIData;

int CreateMatrix(Matrix *Matrix) {
    Matrix->data = (double *) malloc(sizeof(double) * Matrix->rows * Matrix->cols);
    if (!Matrix->data) {
        return ENOMEM;
    }
    return 0;
}

void FreeMatrix(Matrix *Matrix) {
    free(Matrix->data);
}

void PrintMatrix(Matrix matrix) {
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.cols; j++) {
            printf("%f ", matrix.data[i * matrix.cols + j]);
        }
        putchar('\n');
    }
}

void MatrixInitialization(Matrix *Matrix) {
    for (int i = 0; i < Matrix->rows; i++) {
        for (int j = 0; j < Matrix->cols; j++) {
            Matrix->data[i * Matrix->cols + j] = i * Matrix->cols + j + 1;
        }
    }
}

void MPIDataFree(MPIData *MPIdata) {
    MPI_Type_free(&MPIdata->rowType);
    MPI_Type_free(&MPIdata->colType);
    free(MPIdata->dims);
}

void colmRowInitialization(MPIData *MPIdata) {
    int remains_dims[2];
    remains_dims[0] = 0;
    remains_dims[1] = 1;
    MPI_Cart_sub(MPIdata->commGrib, remains_dims, &MPIdata->commCol);

    remains_dims[0] = 1;
    remains_dims[1] = 0;
    MPI_Cart_sub(MPIdata->commGrib, remains_dims, &MPIdata->commRow);
}

void MakeGrib(MPIData *MPIdata, Matrix *matrix1, Matrix *matrix2) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int *dims = calloc(2, sizeof(int));
    MPI_Dims_create(size, 2, dims);
    int periods[2] = {0, 0}; //массив флагов. 0 - нет цикла
    int reorder = 0; //флаг переупоряд. процессов
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &MPIdata->commGrib);
    MPIdata->dims = dims;
    MPI_Type_contiguous(matrix1->cols, MPI_DOUBLE, &MPIdata->rowType);
    MPI_Datatype colTmp;
    MPI_Type_vector(matrix2->rows, 1, matrix2->cols, MPI_DOUBLE, &colTmp);
    MPI_Type_commit(&colTmp);
    MPI_Type_create_resized(colTmp, 0, sizeof(double), &MPIdata->colType);
    MPI_Type_commit(&MPIdata->rowType);
    MPI_Type_commit(&MPIdata->colType);
    colmRowInitialization(MPIdata);
}

int StepOne(MPIData *MPIdata, Matrix *matrix, Matrix *proccesMatrixRow) {
    int rank, size;
    MPI_Comm_rank(MPIdata->commRow, &rank);
    MPI_Comm_size(MPIdata->commRow, &size);

    int *sendCounts = (int *) malloc(sizeof(int) * size);
    if (!sendCounts) {
        perror("Memory error");
        return ENOMEM;
    }
    int *displs = (int *) malloc(sizeof(int) * size);
    if (!displs) {
        free(sendCounts);
        perror("Memory error");
        return ENOMEM;
    }
    int countForProcess = matrix->rows / size;

    for (int i = 0; i < size; ++i) {
        sendCounts[i] = countForProcess;
        displs[i] = (matrix->rows / size) * i * matrix->cols;
    }

    MPI_Scatterv(matrix->data, sendCounts, displs, MPIdata->rowType,
                 proccesMatrixRow->data, countForProcess, MPIdata->rowType, 0, MPIdata->commRow);

    free(sendCounts);
    free(displs);
    return 0;
}

int StepTwo(MPIData *MPIdata, Matrix *matrix, Matrix *proccesMatrixCol) {
    int rank, size;
    MPI_Comm_rank(MPIdata->commCol, &rank);
    MPI_Comm_size(MPIdata->commCol, &size);

    int *sendCounts = (int *) malloc(sizeof(int) * size);
    if (!sendCounts) {
        perror("Memory error");
        return ENOMEM;
    }
    int *displs = (int *) malloc(sizeof(int) * size);
    if (!sendCounts) {
        free(sendCounts);
        perror("Memory error");
        return ENOMEM;
    }
    int countForProcess = matrix->cols / size;

    for (int i = 0; i < size; i++) {
        sendCounts[i] = countForProcess;
        displs[i] = i;
    }

    MPI_Scatterv(matrix->data, sendCounts, displs, MPIdata->colType,
                 proccesMatrixCol->data, countForProcess, MPIdata->colType, 0, MPIdata->commCol);

    free(sendCounts);
    free(displs);
    return 0;
}

void forwardStep1(MPIData *MPIdata, Matrix *matrix, Matrix *proccesMatrixRow) {
    int coords[2];
    int rank, size, countForProcess;
    MPI_Comm_rank(MPIdata->commGrib, &rank);
    MPI_Comm_size(MPIdata->commRow, &size);
    countForProcess = matrix->rows / size;
    MPI_Cart_coords(MPIdata->commGrib, rank, 2, coords);

    if (coords[1] == 0) {
        MPI_Scatter(matrix->data, countForProcess, MPIdata->rowType,
                    proccesMatrixRow->data, countForProcess * matrix->cols, MPI_DOUBLE, 0, MPIdata->commRow);
    }
}

void forwardStep2(MPIData *MPIdata, Matrix *matrix, Matrix *proccesMatrixCol) {
    int coords[2];
    int rank, size, countForProcess;
    MPI_Comm_rank(MPIdata->commGrib, &rank);
    MPI_Comm_size(MPIdata->commCol, &size);
    countForProcess = matrix->cols / size;
    MPI_Cart_coords(MPIdata->commGrib, rank, 2, coords);
    MPI_Datatype colTmp;
    MPI_Type_vector(proccesMatrixCol->rows, 1, proccesMatrixCol->cols, MPI_DOUBLE, &colTmp);
    MPI_Type_commit(&colTmp);
    MPI_Datatype col;
    MPI_Type_create_resized(colTmp, 0, 1 * sizeof(double), &col);
    MPI_Type_commit(&col);
    if (coords[0] == 0) {
        //MPI_Scatter(matrix->data, countForProcess, MPIdata->colType,
        //    proccesMatrixCol->data, countForProcess * matrix->rows, MPI_DOUBLE, 0, MPIdata->commCol);*/
        MPI_Scatter(matrix->data, countForProcess, MPIdata->colType,
                    proccesMatrixCol->data, countForProcess, col, 0, MPIdata->commCol);
    }
}

int forwardEveryone(MPIData *MPIdata, Matrix *matrix1, Matrix *matrix2, Matrix *proccesMatrixRow,
                    Matrix *proccesMatrixCol) {
    forwardStep1(MPIdata, matrix1, proccesMatrixRow);
    forwardStep2(MPIdata, matrix2, proccesMatrixCol);
    //3
    int size, countForProcess;
    MPI_Comm_size(MPIdata->commCol, &size);
    countForProcess = matrix1->rows / size;
    MPI_Bcast(proccesMatrixRow->data, countForProcess * proccesMatrixRow->cols, MPI_DOUBLE, 0, MPIdata->commCol);
    //4
    MPI_Datatype colTmp;
    MPI_Datatype col;
    MPI_Type_vector(proccesMatrixCol->rows, 1, proccesMatrixCol->cols, MPI_DOUBLE, &colTmp);
    MPI_Type_commit(&colTmp);
    MPI_Type_create_resized(colTmp, 0, 1 * sizeof(double), &col);
    MPI_Type_commit(&col);

    MPI_Bcast(proccesMatrixCol->data, countForProcess, col, 0, MPIdata->commRow);

    return 0;
}

void ProccesMatrixInit(MPIData *MPIdata, Matrix *proccesMatrixRow, Matrix *proccesMatrixCol, Matrix *proccesResult,
                       int row1, int col1, int row2, int col2) {
    int rank, size;

    MPI_Comm_rank(MPIdata->commRow, &rank);
    MPI_Comm_size(MPIdata->commRow, &size);
    proccesMatrixRow->cols = col1;
    proccesMatrixRow->rows = row1 / size;

    MPI_Comm_rank(MPIdata->commCol, &rank);
    MPI_Comm_size(MPIdata->commCol, &size);
    proccesMatrixCol->cols = col2 / size;
    proccesMatrixCol->rows = row2;

    proccesResult->rows = proccesMatrixRow->rows;
    proccesResult->cols = proccesMatrixCol->cols;

    CreateMatrix(proccesMatrixRow);
    CreateMatrix(proccesMatrixCol);
    CreateMatrix(proccesResult);
}

void multiplication(Matrix matrixStr, Matrix matrixColm, Matrix result) {
    for (int i = 0; i < matrixStr.rows; ++i) {
        for (int j = 0; j < matrixStr.cols; ++j) {
            for (int k = 0; k < matrixColm.cols; ++k) {
                result.data[i * result.cols + k] +=
                        matrixStr.data[i * matrixStr.cols + j] * matrixColm.data[j * matrixColm.cols + k];
            }
        }
    }
}

void collectAll(Matrix *proccesresult, Matrix *result, MPIData *MPIdata) {
    MPI_Datatype dataMatrix;
    MPI_Type_vector(proccesresult->rows, proccesresult->cols, result->cols, MPI_DOUBLE, &dataMatrix);
    MPI_Type_commit(&dataMatrix);
    MPI_Datatype dataMatrixResized;
    MPI_Type_create_resized(dataMatrix, 0, proccesresult->cols * sizeof(double), &dataMatrixResized);
    MPI_Type_commit(&dataMatrixResized);
    int size;
    MPI_Comm_size(MPIdata->commGrib, &size);
    int *recvCount = malloc(sizeof(int) * size);
    int *displs = malloc(sizeof(int) * size);;
    for (int i = 0; i < MPIdata->dims[0]; i++) // 0 - ãîðèçîíò
    {
        for (int j = 0; j < MPIdata->dims[1]; j++) {
            recvCount[j * MPIdata->dims[0] + i] = 1;
            displs[j * MPIdata->dims[0] + i] = j * MPIdata->dims[0] * proccesresult->rows + i;
        }
    }
    MPI_Gatherv(proccesresult->data, proccesresult->cols * proccesresult->rows, MPI_DOUBLE, result->data, recvCount,
                displs, dataMatrixResized, 0, MPIdata->commGrib);
    free(recvCount);
    free(displs);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    clock_t start, end;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0){
        double cpu_time_used;
        start = clock();
    }
    Matrix matrix1;
    Matrix matrix2;
    Matrix result;
    int row1, col1, row2, col2;
    row1 = 2000;
    col1 = row2 = 2000;
    col2 = 2000;
    matrix1.rows = row1;
    matrix1.cols = col1;
    matrix2.rows = row2;
    matrix2.cols = col2;
    result.rows = row1;
    result.cols = col2;
    CreateMatrix(&matrix1);
    CreateMatrix(&matrix2);
    CreateMatrix(&result);

    if (rank == 0) {
        MatrixInitialization(&matrix1);
        MatrixInitialization(&matrix2);
    }

    MPIData MPIdata;

    Matrix processMatrix1;
    Matrix processMatrix2;
    Matrix processResult;

    MakeGrib(&MPIdata, &matrix1, &matrix2);

    if (row1 % MPIdata.dims[0] != 0 || col2 % MPIdata.dims[1] != 0) {
        if (rank == 0) {
            perror("Error :(");
        }
        FreeMatrix(&matrix1);
        FreeMatrix(&matrix2);
        FreeMatrix(&result);
        MPIDataFree(&MPIdata);
        MPI_Finalize();
        return 0;
    }

    ProccesMatrixInit(&MPIdata, &processMatrix1, &processMatrix2, &processResult, row1, col1, row2, col2);
    forwardEveryone(&MPIdata, &matrix1, &matrix2, &processMatrix1, &processMatrix2);
    multiplication(processMatrix1, processMatrix2, processResult);

    //for (int i = 0; i < 6; i++)
    //{
    //	MPI_Barrier(MPI_COMM_WORLD);
    //	if (rank == i)
    //	{
    //		printf("Process id: %d\n", i);
    //		PrintMatrix(processResult);
    //	}
    //}


    collectAll(&processResult, &result, &MPIdata);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

        printf("Прошло времени: %f секунд\n", cpu_time_used);
        printf("Matrix %d x %d x %d", row1, col1, col2);
        printf("size: %d\n", size);
        //PrintMatrix(result);
    }

    FreeMatrix(&processMatrix1);
    FreeMatrix(&processMatrix2);
    FreeMatrix(&processResult);
    FreeMatrix(&matrix1);
    FreeMatrix(&matrix2);
    FreeMatrix(&result);
    MPIDataFree(&MPIdata);
    MPI_Finalize();
}