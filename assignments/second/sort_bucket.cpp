#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define INF (-1)

/*
bool isSorted(float *data, int size, bool isAscending) {
    int i;
    for (i = 0; i < size; i++) {
        if (isAscending) {
            if (i && data[i] < data[i - 1])
                return false;
        } else if (i && data[i] > data[i - 1]) {
            return false;
        }
    }
    return true;
}
 */

int check(float *data, int nitems) {
    double sum = 0;
    int sorted = 1;
    int i;

    for (i = 0; i < nitems; i++) {
        sum += data[i];
        if (i && data[i] < data[i - 1]) sorted = 0;
    }
    printf("sum=%f, sorted=%d\n", sum, sorted);
    return sorted;
}

int floatComparator(const void *x1, const void *x2) {
    float *f1 = (float *) x1;
    float *f2 = (float *) x2;
    float diff = *f1 - *f2;

    return (diff < 0) ? -1 : 1;
}

int main(int argc, char **argv) {
    int numberProcessor, processorId, N, i;
    float *dSend, *dRecv;
    const float xMin = 1.0;

    MPI_Init(&argc, &argv
    );
    MPI_Comm_size(MPI_COMM_WORLD, &numberProcessor);
    MPI_Comm_rank(MPI_COMM_WORLD, &processorId);

    N = atoi(argv[1]);
    const float xMax = N * 10;
    const int MASTER_RANK = 0; // the master's rank

    int nrecv = N / numberProcessor;
    dSend = (float *) malloc(N * sizeof(float));
    dRecv = (float *) malloc(nrecv * sizeof(float));

    if (processorId == MASTER_RANK) {
        fprintf(stdout,
                "Generating %d numbers to be sorted on %d processors\n", N, numberProcessor);
        for (i = 0; i < N; i++) {
            dSend[i] = drand48() * (xMax - xMin - 1) + xMin;
        }
    }

    double total_s = MPI_Wtime();
    MPI_Scatter(dSend, nrecv, MPI_FLOAT, dRecv, nrecv, MPI_FLOAT, 0, MPI_COMM_WORLD);

    double bucketing_s = MPI_Wtime();
    int bucketCount = numberProcessor;
    float *bucket = (float *)malloc(bucketCount * nrecv * sizeof(float));
    for(i = 0; i < N; i++)
        bucket[i] = INF;

    int *itemCount = (int *)malloc(bucketCount * sizeof(int));
    for (i = 0; i < bucketCount; ++i)
        itemCount[i] = 0;

    float step = (xMax - xMin) / bucketCount;
    for (i = 0; i < N / numberProcessor; i++) {
        int bktno = (int) ((dRecv[i] - xMin) / step);
        int index = bktno * nrecv + itemCount[bktno];
        bucket[index] = dRecv[i];
        ++itemCount[bktno];
    }

    double bucketing_t = MPI_Wtime() - bucketing_s;

    int *recvCount = (int *) calloc(bucketCount, sizeof(int));
    MPI_Alltoall(itemCount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD);

    int *sdispls = (int *) calloc(bucketCount, sizeof(int));
    int *rdispls = (int *) calloc(bucketCount, sizeof(int));
    for (i = 1;  i < bucketCount; i++) {
        sdispls[i] = i * nrecv;
        rdispls[i] = rdispls[i - 1] + recvCount[i - 1];
    }

    float *big_bucket = (float*)malloc(N * sizeof(float));
    MPI_Alltoallv(bucket, itemCount, sdispls, MPI_FLOAT, big_bucket, recvCount, rdispls, MPI_FLOAT, MPI_COMM_WORLD);

    int totalCount = 0;
    for (i = 0; i < bucketCount; i++)
        totalCount += recvCount[i];
    double sorting_s = MPI_Wtime();

    qsort(big_bucket, totalCount, sizeof(float), floatComparator);
    double sorting_t = MPI_Wtime() - sorting_s;

//    memset(recvCount, 0, bucketCount * sizeof(int));

    for (i = 0; i < bucketCount; i ++){
        recvCount[i] = 0;
    }

    MPI_Gather(&totalCount, 1, MPI_INT, recvCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    rdispls[0] = 0;
    for (i = 1; i < bucketCount; i++) {
        rdispls[i] = rdispls[i - 1] + recvCount[i - 1];
    }

    MPI_Gatherv(big_bucket, totalCount, MPI_FLOAT, dSend, recvCount, rdispls, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if (processorId == MASTER_RANK && check(dSend, N)) {
        fprintf(stdout, "total time: %f, parallel: %f\n", MPI_Wtime() - total_s, bucketing_t + sorting_t);
        fprintf(stdout, "The sorted data array from %f to %f\n", dSend[0], dSend[N - 1]);
    }

    MPI_Finalize();
    return 0;
}
