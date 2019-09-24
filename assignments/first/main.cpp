//
// Created by guowenmi on 2019-08-14.
/**
 * calculate pi in parellel environment
 *
 * how to do?
 * it will use master/slave mode.
 *
 * 1. Firstly, master divides different tasks to slaves.
 *      master computes the first Random number and communicates the number and seed to each slave.
 * 2. Each slave computes its random number respectively
 * 3. Each slave generator another corresponding number to compose a point,
 *      and then judge if this point in a unit circle
 * 4. Slaves feedback these results to Master
 * 5. Master count the number in circle and compute the value of pi
 *
 * the steps, 2, 3 and 4, will be paralleled.
 *
 *  It will use Leapfrog method to generate the random number in different processor.
 *  Generating a point by using
 */
//first.cpp Adding numbers using two nodes C++ version
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdint.h>

#include <mpi.h>

typedef unsigned long long int U_LL_INT;
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

    double time00 = MPI_Wtime();
    double time0 = MPI_Wtime();
    double time1 = MPI_Wtime();

    MPI::Init(argc, argv);

    // What is my ID and how many processes are in this pool?
    int myid = MPI::COMM_WORLD.Get_rank(); // get my id
    int numproc = MPI::COMM_WORLD.Get_size(); // get the number of processors

    std::cout << __LINE__ << ", This is id " << myid << " out of " << numproc << std::endl;
    // Get the number the user wants
    U_LL_INT N = atoll(argv[1]); // is there a question due to different types?

    // calculating A and C
    for (int i = 0; i < numproc; i++)
        A = A * a;
    A = A % m;

    U_LL_INT tmpA = 1;
    for (int i = 1; i < numproc; i++) {
        tmpA = tmpA * a;
        C = C + tmpA;
    }
    C = (c * C) % m;

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
    if (x * x + y * y < radius * radius)
        return true;

    return false;
}
