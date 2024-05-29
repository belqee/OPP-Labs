#include <stdio.h>
#include <mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


#define WORKS_SIZE 20
#define COUNT_THREAD 2

typedef struct {
    int* localWorks;              // Локальные задачи для каждого процесса
    int* works;                   // Все задачи, которые будут распределяться
    pthread_mutex_t Mutex;        // Мьютекс для синхронизации потоков
    int LocalWorksCount;          // Количество локальных задач
    int counter;                  // Счетчик выполненных задач
    int refusalCounter;           // Счетчик отказов на запрос задач
    int totalSleepTime;           // Общее время выполнения задач
} SharedState;

void* Getter(void* arg) {
    SharedState* state = (SharedState*)arg;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    while (1) {
        pthread_mutex_lock(&state->Mutex);      // Захватываем мьютекс
        if (state->refusalCounter >= size - 1) { // Если количество отказов достигло размера, выходим из цикла
            pthread_mutex_unlock(&state->Mutex); // Освобождаем мьютекс
            break;
        }
        pthread_mutex_unlock(&state->Mutex);

        int wantWork = -1;                      // Инициализируем запрос на работу
        MPI_Status status;
        MPI_Recv(&wantWork, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status); // Получаем запрос от любого процесса
        int sendRank = status.MPI_SOURCE;
        pthread_mutex_lock(&state->Mutex);
        if (state->LocalWorksCount >= 2) {
            int sharedWork = state->LocalWorksCount / 2;
            MPI_Send(&sharedWork, 1, MPI_INT, sendRank, 0, MPI_COMM_WORLD); // Отправляем количество задач
            int* startShared = state->localWorks + state->counter + sharedWork + state->LocalWorksCount % 2;
            MPI_Send(startShared, sharedWork, MPI_INT, sendRank, 0, MPI_COMM_WORLD); // Отправляем сами задачи
            state->LocalWorksCount -= sharedWork; // Уменьшаем количество локальных задач
        } else {
            wantWork = -1;
            MPI_Send(&wantWork, 1, MPI_INT, sendRank, 0, MPI_COMM_WORLD); // Отправляем сигнал отказа
            state->refusalCounter++;
        }
        pthread_mutex_unlock(&state->Mutex);
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Getter with rank %d closed \n", rank);
    return NULL;
}

void GetWork(int* visited, SharedState* state) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < size; ++i) {
        if (visited[i] != 1) {     // не посещен
            MPI_Send(&rank, 1, MPI_INT, i, 1, MPI_COMM_WORLD); // Отправляем запрос на получение задач
            int countNewWork;
            MPI_Recv(&countNewWork, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Получаем количество задач
            if (countNewWork == -1) {    // получен отказ
                printf("No work for %d from %d \n", rank, i);
                visited[i] = 1;
                continue;
            }
            printf("Work from %d to %d count %d\n", i, rank, countNewWork);
            pthread_mutex_lock(&state->Mutex);
            int* buff = (int*)malloc(countNewWork * sizeof(int)); // память под новые задачи
            pthread_mutex_unlock(&state->Mutex);

            MPI_Recv(buff, countNewWork, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Получаем задачи

            pthread_mutex_lock(&state->Mutex);
            memcpy(state->localWorks, buff, countNewWork); //  копируем задачи в локальный массив
            state->LocalWorksCount = countNewWork;
            pthread_mutex_unlock(&state->Mutex);
            break;
        }
    }
}

void* Work(void* arg) {
    SharedState* state = (SharedState*)arg;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int* visited = (int*)calloc(size, sizeof(int)); // Выделяем память под массив посещенных процессов
    visited[rank] = 1;
    while (1) {
        if (state->LocalWorksCount > 0) {
            pthread_mutex_lock(&state->Mutex);
            printf("Thread %d is doing work\n", rank);
            int workDuration = state->localWorks[state->counter]; // Получаем длительность текущей задачи
            state->LocalWorksCount--;           // Уменьшаем количество локальных задач
            state->counter++;                   // Увеличиваем счетчик выполненных задач
            state->totalSleepTime += workDuration; // Увеличиваем общее время выполнения задач
            pthread_mutex_unlock(&state->Mutex);

            sleep(workDuration);
            printf("process with rank %d: time: %d sec, numb: %d\n", rank, workDuration, state->counter);
        } else {
            pthread_mutex_lock(&state->Mutex);
            state->LocalWorksCount = -1;
            state->counter = 0;
            pthread_mutex_unlock(&state->Mutex);
            GetWork(visited, state);         // Запрашиваем новые задачи
        }
        if (state->LocalWorksCount == -1) {
            break;
        }
    }
    printf("Worker with rank %d is closed \n", rank); //
    free(visited);
    return NULL;
}

void InitWorks(SharedState* state) {
    srand(time(NULL));
    for (int i = 0; i < WORKS_SIZE; ++i) {
        state->works[i] = i + 1;
    }
}

int GetCount(int num, int count, int size) {
    int startIndex = size * num / count;
    int endIndex = size * (num + 1) / count;
    int countRow = endIndex - startIndex;
    return countRow;
}

void ShareWork(SharedState* state) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int* counts = (int*)malloc(size * sizeof(int));
    for (int i = 0; i < size; ++i) {
        counts[i] = GetCount(i, size, WORKS_SIZE);
    }
    int* shift = (int*)malloc(size * sizeof(int));
    shift[0] = 0;
    for (int i = 1; i < size; ++i) {
        shift[i] = shift[i - 1] + counts[i - 1];   //смещение для каждого процесса
    }
    pthread_mutex_lock(&state->Mutex);             // Захватываем мьютекс
    MPI_Scatterv(state->works, counts, shift, MPI_INT, state->localWorks, counts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    pthread_mutex_unlock(&state->Mutex);
    free(shift);
    free(counts);
}

int main(int argc, char* argv[]) {
    int rank, size, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        printf("Error: Required thread level not provided\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    SharedState* state = (SharedState*)malloc(sizeof(SharedState)); // Выделяем память под структуру состояния
    pthread_mutex_init(&state->Mutex, NULL);
    state->LocalWorksCount = 0;
    state->counter = 0;
    state->refusalCounter = 0;
    state->totalSleepTime = 0;
    double start_time, end_time;
    start_time = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    state->works = NULL;
    state->localWorks = NULL;

    int totalWorksSum = 0;

    if (rank == 0) {
        state->works = (int*)malloc(WORKS_SIZE * sizeof(int));
        InitWorks(state);
        for (int i = 0; i < WORKS_SIZE; ++i) {
            totalWorksSum += state->works[i];
        }
        printf("Total sum of works: %d\n", totalWorksSum);
    }

    for (int iteration = 0; iteration < 2; ++iteration) {
        state->LocalWorksCount = GetCount(rank, size, WORKS_SIZE);
        state->localWorks = (int*)malloc(state->LocalWorksCount * sizeof(int));
        state->counter = 0;
        ShareWork(state);                   // Распределяем
        pthread_t threads[COUNT_THREAD];    // Массив потоков

        pthread_create(&threads[0], NULL, Getter, (void*)state); //  поток геттера
        pthread_create(&threads[1], NULL, Work, (void*)state);   // поток рабочего

        for (int i = 0; i < COUNT_THREAD; ++i) {
            pthread_join(threads[i], NULL); // Ожидаем завершения потоков
        }

        free(state->localWorks);
    }

    if (rank == 0) {
        free(state->works);
    }

    end_time = MPI_Wtime();
    printf("Execution time of process %d = %f\n", rank, end_time - start_time);

    pthread_mutex_destroy(&state->Mutex);
    free(state);
    MPI_Finalize();
    return 0;
}
