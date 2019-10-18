//
// Created by guowenmi on 2019-10-01.
/**
 * parallel bucket sort
 *
 * how to do?
 * it will use master/slave mode.
 *
 * 1. First, Master generates random numbers needed to be sorted. (only in Master)
 * 2. Second, Master scatter numbers to each processors evenly.
 * 3. Each processor divides its own numbers into proper small buckets.
 * 4. Each processor scatters its own small buckets to all other processors
 *    and gather its own proper numbers from all other processors. (by call alltoall() and alltoallv() functions)
 * 5. Each processor sorts its own numbers
 * 6. Master gather sorted numbers from all other processes (large bucket).
 * 7. Master combines these large bucket in proper order.
 *
 * the steps, 3, 4 and 5, will be paralleled.
 *
 */
//first.cpp Adding numbers using two nodes C++ version
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdint.h>

#include <mpi.h>

typedef unsigned long long int U_LL_INT;

int MASTER_ID = 0;

U_LL_INT getRandom ();// use to generate random.

// for generate random
U_LL_INT a = 1664525;
U_LL_INT m = 4294967296;//2 << 32;
U_LL_INT c = 1013904223;
U_LL_INT n_prev = 12345; // seed, it will be used as n[0]
U_LL_INT currRandom = 0;
U_LL_INT n_next;



U_LL_INT sumInCircle = 0;
bool isInCircle = false;
double pi;


void getRandom(U_LL_INT sampleNumber);
bool pointIsInCircle (U_LL_INT randomNumber);
U_LL_INT getRandomInLeapfrog(U_LL_INT random0);
U_LL_INT countInCircleNumber (U_LL_INT seed, int processorid, int processorNumber, U_LL_INT loopNumber);

U_LL_INT A = 1;
U_LL_INT C = 1;

// 5 * 8
U_LL_INT randomArray [40];

//always use argc and argv, as mpirun will pass the appropriate parms.
int main(int argc,char* argv[]) {

    std::cout << __LINE__ << ", hello, This is at the begin of main " << std::endl;

    double starttime, endtime;

    MPI::Init(argc, argv);

    // What is my ID and how many processes are in this pool?
    int myid = MPI::COMM_WORLD.Get_rank(); // get my id
    int num_processor = MPI::COMM_WORLD.Get_size(); // get the number of processors

    std::cout << __LINE__ << ", This is id " << myid << " out of " << num_processor << std::endl;

    // Get the number the user wants
    U_LL_INT num_needtosort = atoll(argv[1]);   // the number of numbers
    U_LL_INT send_count = num_needtosort / num_processor;   //

    U_LL_INT array_random_number[num_needtosort];   // create an array with length = num_needtosort
    int *recvbuf = new int[send_count];     // N / P numbers in per processor

    // if it is Master, it needs to generate random numbers array, and then scatter to all processors.
    if (myid == 0) {
        // fill the array using random numbers
        for (int i = 0; i < num_needtosort; i++)
            array_random_number[i] = getRandom();
    }

    // starting time calculation of the sort
    starttime = MPI_Wtime();

    // every processs need to call this function
    // (send_buf, send_count, send_type, recv_buf, rev_count, recv_buf, recv_type, root_id)
    MPI::COMM_WORLD.Scatter(array_random_number, send_num, MPI_UNSIGNED_LONG_LONG, recvbuf, send_num, MPI_UNSIGNED_LONG_LONG, MASTER_ID);





    ///////////////////// split line

    getRandom (numproc);// generate the random as seeds
    U_LL_INT loopNumber = N ; // numproc ;

    time1 = MPI_Wtime();
    std::cout << __LINE__ << ", getRandom cost time = " << (time1 - time0)<< std::endl;

    if (myid == 0) { // master. need to distribute and gather
        // send seeds to all slaves
        for (int i = 1; i < numproc; i ++){
            // MPI_Send(void* data, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator)
            MPI_Send(&randomArray[i], 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
        }

        double startTime = MPI_Wtime();

        // do master's task
        sumInCircle = countInCircleNumber(randomArray[0], myid, numproc, loopNumber);

        std::cout << __LINE__ << ", myid = 0, sumInCircle = " << sumInCircle << ", N = " << N << std::endl;

        // gather all slaves' result and deal with
        for (int i = 1; i < numproc; i ++){
            //MPI_Recv(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator, MPI_Status* status)
            U_LL_INT slaveInCircleSum = 0;
            MPI_Recv(&slaveInCircleSum, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sumInCircle = sumInCircle + slaveInCircleSum;

            std::cout << __LINE__ << ", MPI_Recv from = " << i << ", slaveInCircleSum = " << slaveInCircleSum << ", sumInCircle = " << sumInCircle << std::endl;
        }

        std::cout << __LINE__ << ", sumInCircle = " << sumInCircle << ", N = " << N << std::endl;
        pi = 4.0 * sumInCircle / N ;

        double endTime = MPI_Wtime ();
        std::cout << __LINE__ << ", The pi is " << pi << ", cost time is " << (endTime - startTime) << std::endl;
    } else /* if (myid == 1) */ {
        // slaves. are processes of all slaves same? need to test
        U_LL_INT sumInCircleSlave = 0;

        //MPI_Recv(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator, MPI_Status* status)
        MPI_Recv (&currRandom, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        U_LL_INT result = countInCircleNumber (currRandom, myid, numproc, loopNumber);

        // MPI_Send(void* data, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator)
        MPI_Send(&result, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);

        std::cout << __LINE__ << " I am " << myid << ", have sent to Master" << std::endl;
    }

    MPI::Finalize();

    time1 = MPI_Wtime();
    std::cout << __LINE__ << ", The pi is " << pi << ", cost time is " << (time1 - time00) << std::endl;
}


U_LL_INT countInCircleNumber (U_LL_INT randomSeed, int processorId, int numproc, U_LL_INT loopNumber){

    std::cout << __LINE__ << ", processorId = " << processorId << ", randomSeed = " << randomSeed << ", numproc = " << numproc << ", loopNumber = " << loopNumber << std::endl;

    U_LL_INT sum = 0;
    for (U_LL_INT index = processorId; index < loopNumber;){

        //    std::cout << __LINE__ << ", index = " << index << ", currRandom = " << currRandom << std::endl;
        double time1 = MPI_Wtime();
        // Partial result for node 0
        if (index < numproc) {
            currRandom = randomSeed; //randomArray[processorId];
        } else {
            currRandom = getRandomInLeapfrog (currRandom);
        }

        double time0 = MPI_Wtime();
        //    std::cout << __LINE__ << ", currRandom = " << currRandom << ", cost time = " << (time1 - time0) << std::endl;

        bool isInCircleTmp = pointIsInCircle(currRandom);
        if (isInCircleTmp)
            sum = sum + 1;

        time1 = MPI_Wtime();
        //    std::cout << __LINE__ << ", isInCircleTmp = " << isInCircleTmp << ", sumInCircle = " << sum << ", pointIsInCircle cost time = " << (time1 - time0) << std::endl;

        index = index + numproc;
    }

    std::cout << __LINE__ << ", processorId = " << processorId << ", sumInCircle = " << sum << std::endl;

    return sum;
}

// generate one random number
void getRandom() {
        n_next = (a * n_prev + c) % m;
        //    std::cout << __LINE__ << ", n_prev = " << n_prev << ", n_next = " << n_next << std::endl;
        n_prev = n_next;
        return n_next;
}

// generate randomNumber random
void getRandom(U_LL_INT randomNumber) {

    for (int i = 0; i < randomNumber; i++) {
        n_next = (a * n_prev + c) % m;

        //    std::cout << __LINE__ << ", n_prev = " << n_prev << ", n_next = " << n_next << std::endl;

        randomArray[i] = n_next;
        n_prev = n_next;
    }
}

/**
 *  Formula: x[i+k] = (Ax[i] + C) mod m
 *  k = the number of processors
 * @param random0
 * @param k
 * @return
 */
U_LL_INT getRandomInLeapfrog(U_LL_INT random_prev) {
//    std::cout << __LINE__ << ", A = " << A << ", random_prev = " << random_prev << ", C = " << C << std::endl;

    return (A * random_prev + C) % m;
//    return result;
}

U_LL_INT radius = 65536; //2 << 16;

bool pointIsInCircle (U_LL_INT randomNumber){

    U_LL_INT x = randomNumber % radius;
    U_LL_INT y = randomNumber / radius;

//    std::cout << __LINE__ << ", randomNumber = " << randomNumber << ", radius = " << radius << ", x = " << x << ", y = " << y << std::endl;

//    if ((x / radius) ** 2 + (y / radius) ** 2 <= 1
    if (x * x + y * y < (radius + 0.5) * (radius + 0.5))
        return true;

    return false;
}
