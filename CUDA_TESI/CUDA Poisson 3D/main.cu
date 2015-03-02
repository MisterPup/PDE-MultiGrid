#include "Grid3D.h"
#include "MultiGrid3D.h"
#include "inclusion.h"

//void checkCUDAError(const char *msg);

int main()
{
	/* **********TIMING********** */
	time_t start, stop;
	clock_t ticks;
	time(&start);
	/* **********TIMING********** */	

	float range[] = {0, 1, 0, 1, 0, 1}; //se aumento il range devo aumentare anche i punti 
	//QUI IL RANGE NON DEVE ESSERE MODIFICATO

	int equalSize = 257; //controlla che non sia minore di 5
	int finestGridSizeXYZ[] = {equalSize, equalSize, equalSize};
	int v0 = 2;
	int v1 = 3000;
	int v2 = 3000;

	MultiGrid3D multiGrid3D(finestGridSizeXYZ, range);	


	/* ******************* TEST FULL MULTIGRID CYCLE ***************** */
	
	int finestGridID = 0;

	//multiGrid3D.VCycle(finestGridID, v1, v2);
	multiGrid3D.FullMultiGridVCycle(finestGridID, v0, v1, v2);
	
	//multiGrid3D.PrintAllGrids_v();
	//multiGrid3D.PrintAllGrids_f();
	//multiGrid3D.PrintGrid(0);
	//multiGrid3D.PrintDiff();
	/* ******************* TEST FULL MULTIGRID CYCLE ****************** */


	
	/* **********TIMING********** */
	ticks = clock();
	time(&stop);
	printf("finestGridSize: %d\n", equalSize);
	printf("Used %0.2f seconds of CPU time. \n", (double)ticks/CLOCKS_PER_SEC);
	printf("Finished in about %.0f seconds. \n", difftime(stop, start));
	/* **********TIMING********** */

	return 0;	
}

/*
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }                         
}*/
