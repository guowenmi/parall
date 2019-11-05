/*                                                
*/
	
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  <tgmath.h> 

#include "mpi.h"





/*****************************************************************************/
// Routines used in the sequential implementation of the bucket sort
float* create_buckets(int nbuckets, int nitems);
void bucket_sort(float *data, int ndata, float x1, float x2, int nbuckets,
		 float* bucket);

int check(float *data,int nitems) {
  double sum=0;
  int sorted=1;
  int i;

//	printf("The numbers are:\n");
	
  for(i=0;i<nitems;i++) {
 // 	printf("%f \n", data[i]);
	
     sum+=data[i];
     if(i && data[i]<data[i-1]) sorted=0;
  }
  printf("\nsum=%f, sorted=%d\n",sum,sorted);

  return sorted;
}

int compare(const void* x1, const void* x2);
/*****************************************************************************/

// Sequential implementation of the bucket sort routine. The full
// range x1 to x2 will be divided into a number of equally spaced
// subranges according to the number of buckets. All the buckets are
// contained in the single one dimensional array "bucket".
void bucket_sort(float *data, int ndata, float x1, float x2, int nbuckets,
		 float *bucket) 
{
  int i, count;

  // The range covered by one bucket
  float stepsize = (x2 - x1) / nbuckets;

  // The number of items thrown into each bucket. We would expect each
  // bucket to have a similar number of items, but they won't be
  // exactly the same. So we keep track of their numbers here.
  int* nitems = (int*)malloc(nbuckets * sizeof(int));
  for (i = 0; i < nbuckets; ++i) nitems[i] = 0;

  // Toss the data items into the correct bucket
  for (i = 0; i < ndata; ++i) {

    // What bucket does this data value belong to?
    int bktno = (int)floor((data[i] - x1) / stepsize);
    int idx = bktno * ndata + nitems[bktno];

    //printf("DATA %d %f %d %d\n", i, data[i], bktno, idx);

    // Put the data value into this bucket
    bucket[idx] = data[i];
    ++nitems[bktno];
  }
  
  // Sort each bucket using the standard library qsort routine. Note
  // that we need to input the correct number of items in each bucket
  count = 0;
  for (i = 0; i < nbuckets; ++i) {
    if(nitems[i]) {
      qsort(&bucket[i*ndata], nitems[i], sizeof(float), compare);
      memcpy(data,&bucket[i*ndata],nitems[i]*sizeof(float));
      data+=nitems[i];
    }
  }

  // Don't need the number of items anymore
  free(nitems);
  
}

/*****************************************************************************/

// Create a data array to hold the given number of buckets for the
// given number of total data items. All buckets are held contiguously
// in the
float* create_buckets(int nbuckets, int nitems)
{
  int i;

  int ntotal = nbuckets * nitems;

  // Pointer to an array of more pointers to each bucket
  float* bucket = (float*)calloc(ntotal, sizeof(float*));
  for (i=0; i<ntotal; ++i) bucket[i] = 0;

  // return the address of the array of pointers to float arrays
  return bucket;
}

/*****************************************************************************/

// The comparison function to use with the library qsort
// function. This will tell qsort to sort the numbers in ascending
// order.
int compare(const void* x1, const void* x2) {
  const float* f1 = (const float*)x1;
  const float* f2 = (const float*)x2;
  float diff = *f1 - *f2;

  return (diff < 0) ? -1 : 1;
}

/*****************************************************************************/

/*****************************************************************************/


void distribute_among_small_buckets(float *data, int ndata, float x1, float x2, int nbuckets,
		 float *bucket, int *numOfItemsInEachSmallBucket) 
{
  int i, count;

  // The range covered by one bucket
  float stepsize = (x2 - x1) / nbuckets;

  for (i = 0; i < nbuckets; ++i) numOfItemsInEachSmallBucket[i] = 0;

  // Toss the data items into the correct bucket
  for (i = 0; i < ndata; ++i) {

    // What bucket does this data value belong to?
    int bktno = (int)floor((data[i] - x1) / stepsize);
    int idx = bktno * ndata + numOfItemsInEachSmallBucket[bktno];

    //printf("DATA %d %f %d %d\n", i, data[i], bktno, idx);

    // Put the data value into this bucket
    bucket[idx] = data[i];
    ++numOfItemsInEachSmallBucket[bktno];
  }  
}


int main(int argc,char *argv[]) {


  // Specify here the full range that the data numbers can cover. We
  // may assume that all the numbers are positive definite
  const float xmin = 10.0;
  const float xmax = 250000;
  int i;

  if(argc != 2) 
  {
  	printf("The command line should be like this:\n %s {number of random numbers to be generated} \n", argv[0]);
	return 0;
  }

  int nitems=atoi(argv[1]);
  
  MPI_Init(&argc, &argv);
  
  // What is my ID and how many processes are in this pool?
  int myid;
  int numproc;
  MPI_Comm_size(MPI_COMM_WORLD, &numproc); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); 

	if (nitems % numproc != 0)
	{
		printf("The number of random numbers to be generated must be divisible by that of processors\n");
		MPI_Finalize();
		return 0;
	}
	
  double t_total_start, t_total_end, t_parallel_start1, t_parallel_end1, t_parallel_start2, t_parallel_end2;
  
  t_total_start =	MPI_Wtime();
  

	// "data" is for storing the original random numbers.
  float *data = NULL;

  if (myid == 0)
  {
  	data = (float*)malloc(nitems*sizeof(float));
	
	  for(i=0;i<nitems;i++)
	    data[i]=drand48()*(xmax-xmin-1)+xmin;

	printf("Before sorting, calculate the sum of the random numbers:");
	check(data,nitems);
  }
  
  // A largeBucket will be used for accepting initial numbers distributed by the master using MPI_Scatter,
  // as well as collecting all the numbers whose values are within a specific range later
	float *largeBucket = (float*)malloc(nitems*sizeof(float));	

	// Scatter random numbers to each processor's large bucket, including the master's own large bucket.
	// Each processor receives nitems/numproc numbers
	  MPI_Scatter(data, nitems/numproc, MPI_FLOAT, largeBucket, nitems/numproc, MPI_FLOAT, 0, MPI_COMM_WORLD);

	t_parallel_start1 = MPI_Wtime();
	
	// Each processor has numproc small buckets. each of which can hold nitems/numproc numbers at most.
	float *smallBuckets = (float*)malloc(nitems*sizeof(float));	

	int *numOfItemsInEachSmallBucket = (int*)malloc(numproc*sizeof(int));

	// distribute from my large bucket among my own small buckets, but not sort them
	distribute_among_small_buckets(largeBucket, nitems/numproc, xmin, xmax, numproc, smallBuckets, numOfItemsInEachSmallBucket);

	int *sdispls = (int*) malloc(numproc * sizeof(int));
	int *recvcounts = (int*) malloc(numproc * sizeof(int));
	int *rdispls = (int*) malloc(numproc * sizeof(int));

	t_parallel_end1 = MPI_Wtime();
	
	// notify other processors (including my own processor) how many numbers their large buckets need to be prepared to receive from me
	MPI_Alltoall(numOfItemsInEachSmallBucket, 1, MPI_INT, recvcounts, 1, MPI_INT, 
			 MPI_COMM_WORLD);

	int totalReceivedInLargeBucket = 0;
	
	for (int i=0; i<numproc; i++)
	{
		sdispls[i] = i*(nitems/numproc);
		rdispls[i] = totalReceivedInLargeBucket;
		totalReceivedInLargeBucket += recvcounts[i];
	}

	// transmit the numbers in my small buckets to correspoding large buckets. i.e.
	// Those in my 1st small buckets go to the large bucket of processor 0 
	// Those in my 2nd small buckets go to the large bucket of processor 1
	// Those in my 3rd small buckets go to the large bucket of processor 2
	// ....
	MPI_Alltoallv(smallBuckets, numOfItemsInEachSmallBucket,
					  sdispls, MPI_FLOAT, largeBucket,
					  recvcounts, rdispls, MPI_FLOAT,
					  MPI_COMM_WORLD);

	t_parallel_start2 =  MPI_Wtime();
	
	// sort the numbers in my large bucket
	qsort(largeBucket, totalReceivedInLargeBucket, sizeof(float), compare);

	t_parallel_end2 = MPI_Wtime();
	
	// notify the master how many numbers there are in my large bucket
	MPI_Gather(&totalReceivedInLargeBucket, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (int i=1; i<numproc; i++)
	{		
		rdispls[i] = rdispls[i-1] + recvcounts[i-1];
	}
	
	// The master collects all the numbers from each processor's large bucket, which have already been sorted.
	MPI_Gatherv(largeBucket, totalReceivedInLargeBucket, MPI_FLOAT, data, recvcounts, rdispls, MPI_FLOAT, 0, MPI_COMM_WORLD);  

	t_total_end = MPI_Wtime();
	
	if (myid == 0)
	{
		printf("\nSorting finished. Now, calculate the sum again to check if they are sorted without an error:");
		check(data,nitems);
		printf("\nIs the sum of the numbers equal to the original sum?\n If not, there must be a bug in the program.\n"
			"Note that they need not be exactly the same, due to finite precision of floating-point arithmetic");
	}
  
  MPI_Finalize();

  if (myid == 0)
  {
  	free(data);
	printf("\nTotal time is %f. The amount of time that parallel work takes is %f\n", t_total_end-t_total_start,
		(t_parallel_end1-t_parallel_start1) + (t_parallel_end2-t_parallel_start2));
  }

  free(largeBucket);
  free(smallBuckets);
  free(numOfItemsInEachSmallBucket);
  free(sdispls);
  free(recvcounts);
  free(rdispls);
  
  return 0;
}

