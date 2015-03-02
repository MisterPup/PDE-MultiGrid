#include "Grid1D.h"
#include "MultiGrid1D.h"
#include "inclusion.h"

int main()
{
	/* **********TIMING********** */
	time_t start, stop;
	clock_t ticks;
	time(&start);
	/* **********TIMING********** */	

	float range[] = {0, 1}; //range = {xa,xb} xb>xa NON MODIFICARE

	int finestGridSize = 8193; 
	int v0 = 2;
	int v1 = 1000;
	int v2 = 1000;

	MultiGrid1D multiGrid1D(finestGridSize, range);

	/* ******************* TEST FULL MULTIGRID CYCLE ***************** */

	int finestGridID = 0;
	multiGrid1D.FullMultiGridVCycle(finestGridID, v0, v1, v2);
	//multiGrid1D.PrintGrid(0);
	//multiGrid1D.PrintDiff();

	/* ******************* TEST FULL MULTIGRID CYCLE ****************** */
	
	/* **********TIMING********** */
	ticks = clock();
	time(&stop);
	printf("finestGridSize: %d\n", finestGridSize);
	printf("Used %0.2f seconds of CPU time. \n", (double)ticks/CLOCKS_PER_SEC);
	printf("Finished in about %.0f seconds. \n", difftime(stop, start));
	/* **********TIMING********** */

	return 0;	
}
