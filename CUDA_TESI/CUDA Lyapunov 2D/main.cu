#include "Grid2D.h"
#include "MultiGrid2D.h"
#include "inclusion.h"

int main()
{
	/* **********TIMING********** */
	time_t start, stop;
	clock_t ticks;
	time(&start);
	/* **********TIMING********** */	

	float range[] = {0, 20, 0, 20}; //se aumento il range devo aumentare anche i punti
	int sizeX_A = 2; //matrice 2x2 
	float* A = (float*)malloc(sizeX_A*sizeof(float));
	A[0] = -1.0f;
	A[1] = -2.0f;
	A[2] = 0.0f;
	A[3] = -3.0f;

	int alfa = 2;
	int finestGridSize = 65; //65	
	int v0 = 2;
	int v1 = 500;
	int v2 = 500;	

	MultiGrid2D multiGrid2D(finestGridSize, range, A, sizeX_A, alfa);

	/* ******************* TEST FULL MULTIGRID CYCLE ***************** */

	int finestGridID = 0;	
	//multiGrid2D.VCycle(finestGridID, v1, v2);	
	multiGrid2D.FullMultiGridVCycle(finestGridID, v0, v1, v2);
	multiGrid2D.PrintMeanAbsoluteError();
	//multiGrid2D.PrintDiff();
	//multiGrid2D.PrintGrid(0);

	/* ******************* TEST FULL MULTIGRID CYCLE ***************** */

	
	/* **********TIMING********** */
	ticks = clock();
	time(&stop);
	//printf("range: {%f, %f, %f, %f}\n", range[0], range[1], range[2], range[3]);
	printf("finestGridSize: %d\n", finestGridSize);
	printf("v0: %d\n", v0);
	printf("v1: %d\n", v1);
	printf("v2: %d\n", v2);
	printf("Used %0.2f seconds of CPU time. \n", (double)ticks/CLOCKS_PER_SEC);
	printf("Finished in about %.0f seconds. \n", difftime(stop, start));
	/* **********TIMING********** */

	return 0;	
}
