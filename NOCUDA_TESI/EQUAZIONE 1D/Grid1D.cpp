#include "Grid1D.h"
#include "inclusion.h"

Grid1D::Grid1D(int sizeX, float range[]) //range = {xa,xb} xb>xa
{
	assert(sizeX > 0);
	assert((sizeX-1)%2 == 0); //size=(2^k) - 1
	this->sizeX = sizeX;

	assert(range[1] > range[0]);
	float x_range = range[1] - range[0];
	x_a = range[0];
	x_b = range[1];
	
	h_x = x_range/(float)(sizeX-1);

	h_v = (float*)malloc(sizeX*sizeof(float));
	h_f = (float*)malloc(sizeX*sizeof(float));

	InitV();
	InitF();
}

Grid1D::~Grid1D()
{
	free(h_v);
	free(h_f);
}

void Grid1D::InitV()
{
	h_v[0] = (exp(x_a) + x_a - 3)/(1 + exp(-x_a));
	h_v[sizeX-1] = (exp(x_b) + x_b - 3)/(1 + exp(-x_b));
}

void Grid1D::InitF()
{
	for(int posX = 0; posX < sizeX; posX++)
	{
		float xj = x_a + posX*h_x;
		h_f[posX] = exp(xj);
	}
}

//ha senso solo per griglia più fine (finest)
void Grid1D::PrintDiffApproxReal(int diff_fd)
{
	for(int posX = 0; posX < sizeX; posX++)
	{
		float xj = x_a + posX*h_x;
		float realsol = (exp(xj) + xj - 3)/(1 + exp(-xj));
		float approxsol = h_v[posX];
		float diff = approxsol - realsol;

		char* log = (char*)malloc(100);
		sprintf(log,"xj: %f diff: %f\n", xj, diff);
		write(diff_fd, log, strlen(log));
	}
}

void Grid1D::PrintGrid_v(int logfd)
{
	for(int posX = 0; posX < sizeX; posX++)
	{
		float xj = x_a + posX*h_x;
		float realsol = (exp(xj) + xj - 3)/(1 + exp(-xj));;
		float approxsol = h_v[posX];
		char* log = (char*)malloc(100);
		sprintf(log,"xj: %f approxSol: %f  realSol: %f\n", xj, approxsol, realsol);
		write(logfd, log, strlen(log));
	}
}

void Grid1D::PrintGrid_f(int logfd)
{
	for(int posX = 0; posX < sizeX; posX++)
	{
		float xj = x_a + posX*h_x;

		char* log = (char*)malloc(200);
		sprintf(log, "xj: %f f:%f\n", xj, h_f[posX]);
		write(logfd, log, strlen(log));
	}
}
