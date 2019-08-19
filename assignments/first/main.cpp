//
// Created by guowenmi on 2019-08-14.
/**
 * calculate pi in parellel environment
 *
 * how to do?
 * this will use master/slave mode.
 *
 * 1. Firstly, master divides different tasks to slaves.
 *      master computes the first Random number and communicates the number and seed to each slave.
 * 2. Each slave computes its random number respectively
 * 3. Each slave generator another corresponding number to compose a point,
 *      and then decide if this point in a unit circle
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

void getRandom(int sampleNumber);
bool pointIsInCircle (U_LL_INT randomNumber);
U_LL_INT getRandomInLeapfrog(U_LL_INT random0);

U_LL_INT A = a;
U_LL_INT C = 1;

// 5 * 8
U_LL_INT randomArray [40];

//always use argc and argv, as mpirun will pass the appropriate parms.
int main(int argc,char* argv[]) {
    MPI::Init(argc, argv);

    // What is my ID and how many processes are in this pool?
    int myid = MPI::COMM_WORLD.Get_rank(); // get my id
    int numproc = MPI::COMM_WORLD.Get_size(); // get the number of processors

    std::cout << "numproc = " << numproc << std::endl;

    std::cout << "This is id " << myid << " out of " << numproc << std::endl;

    // Get the number the user wants
    // is also the number of sample
    U_LL_INT N = atoi(argv[1]); // is there a question due to different types?

    std::cout << __LINE__ << ", a = " << a << ", m = " << m << std::endl;
    double time0 = MPI_Wtime();

    // calculating A and C
    for (int i = 1; i < numproc; i++)
        A = A * a;
    A = A % m;
    double time1 = MPI_Wtime();
    std::cout << __LINE__ << ", A = " << A << ", cost time = " << (time1 - time0) << std::endl;


    U_LL_INT tmpA = 1;
    for (int i = 1; i < numproc; i++) {
        tmpA = tmpA * a;
        C = C + tmpA;
    }
    C = (c * C) % m;

    time0 = MPI_Wtime();
    std::cout << __LINE__ << ", C = " << C << ", cost time = " << (time1 - time0)<< std::endl;

    getRandom (numproc);// generate the random as seeds

    time1 = MPI_Wtime();
    std::cout << __LINE__ << ", getRandom =  cost time = " << (time1 - time0)<< std::endl;

    if (myid == 0) { // master. need to distribute and gather

        double startTime = MPI_Wtime();
        std::cout << __LINE__ << std::endl;

        for (unsigned long int index = 0; index < N / numproc; index ++){
            // Master sends currRandom to slaves
            // MPI_Send(void* data, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator)
            for (int i = 1; i < numproc; i++){
//                MPI::COMM_WORLD.Send(&randomArray, 1, MPI::MPI_UNSIGNED_LONG_LONG, i, 0);
                MPI_Send(&randomArray, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
            }

            std::cout << __LINE__ << ", currRandom = " << currRandom << std::endl;
            time1 = MPI_Wtime();
            // Partial result for node 0
            if (index < numproc) {
                currRandom = randomArray[0];
            }
            else {
                currRandom = getRandomInLeapfrog (currRandom);
            }

            time0 = MPI_Wtime();
            std::cout << __LINE__ << ", currRandom = " << currRandom << ", cost time = " << (time1 - time0) << std::endl;

            if (pointIsInCircle(currRandom))
                sumInCircle = sumInCircle + 1;
            std::cout << __LINE__ << std::endl;


            time1 = MPI_Wtime();
            std::cout << __LINE__ << ", sumInCircle = " << sumInCircle << ", pointIsInCircle cost time = " << (time1 - time0) << std::endl;

            //Master waits to receive 'sum1' from slave
            //MPI_Recv(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator, MPI_Status* status)
            for (int i = 1; i < numproc; i++) {
                MPI::COMM_WORLD.Recv(&isInCircle, 1, MPI::INT, i, 0);
                if (isInCircle)
                    sumInCircle = sumInCircle + 1;
            }
        }
        std::cout << __LINE__ << ", N = " << N << std::endl;

        pi = sumInCircle / N ;

        double endTime = MPI_Wtime ();

        std::cout << "The pi is " << pi << ", cost time is " << (endTime - startTime) << std::endl;

    } else /* if (myid == 1) */ {
        // slaves. are processes of all slaves same? need to test

        for (unsigned long int index = 0; index < N / numproc; index ++) {
            // Slave waits to receive 'currRandom' from master
            //MPI_Recv(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator, MPI_Status* status)
            MPI::COMM_WORLD.Recv(&randomArray, 1, MPI::INT, 0, 0); // MPI::COMM_WORLD.Recv

            if (index < numproc) {
                currRandom = randomArray[myid];
            } else {
                currRandom = getRandomInLeapfrog(currRandom);
            }
            U_LL_INT slaveRandom = getRandomInLeapfrog(currRandom);

            time0 = MPI_Wtime();

            bool isInCircle = pointIsInCircle(slaveRandom);

            time1 = MPI_Wtime();
            std::cout << __LINE__ << ", index = " << index << ", isInCircle = " << isInCircle << ", pointIsInCircle cost time = " << (time1 - time0) << std::endl;

            // Slave sends 'isInCircle' to master
            MPI_Send(&isInCircle, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    MPI::Finalize();
}


void getRandom(int randomNumber) {

    for (int i = 0; i < randomNumber; i++) {
        n_next = (a * n_prev + c) % m;
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
    return (A * random_prev + C) % m;
//    return result;
}

U_LL_INT radius = 2 << 16;

bool pointIsInCircle (U_LL_INT randomNumber){

    U_LL_INT x = randomNumber % radius;
    U_LL_INT y = randomNumber / radius;
//    if ((x / radius) ** 2 + (y / radius) ** 2 <= 1
    if (x * x + y * y <= radius * radius)
        return true;

    return false;
}


