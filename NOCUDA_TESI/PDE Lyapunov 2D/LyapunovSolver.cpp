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

	float range[] = {0, 1, 0, 1}; //se aumento il range devo aumentare anche i punti

	int size_A = 2;

	float* A = (float*)malloc(size_A*size_A*sizeof(float)); //matrice quadrata

	//int idx = posX + posY * sizeX;
	A[0] = -1.0f; //A[0][0]
	A[1] = -2.0f; //A[0][1]
	A[2] = 0.0f;  //A[1][0]
	A[3] = -3.0f; //A[1][1]

	int equalSize = 1025;
	int finestGridSizeXY[] = {equalSize, equalSize};

	int alfa = 2;
	int v0 = 1;
	int v1 = 500;
	int v2 = 500;	

	//Grid2D grid2D(finestGridSizeXY, range);
	//int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	//grid2D.PrintGrid_v(logfd);

	MultiGrid2D multiGrid2D(finestGridSizeXY, range, A, size_A, alfa);

	/* ******************* TEST FULL MULTIGRID CYCLE ***************** */

	int finestGridID = 0;
	//multiGrid2D.VCycle(finestGridID, v1, v2);
	multiGrid2D.FullMultiGridVCycle(finestGridID, v0, v1, v2);
	multiGrid2D.PrintDiff();
	//multiGrid2D.PrintResidual(0);

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
