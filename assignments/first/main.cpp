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
#include <stdatomic.h>
#include <stdint.h>

#include "mpi.h"


//always use argc and argv, as mpirun will pass the appropriate parms.
int main(int argc,char* argv[]) {
    return 0;
}




