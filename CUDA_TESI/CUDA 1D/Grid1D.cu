#include "Grid1D.h"
#include "inclusion.h"

__global__ void CUDASetGrids(float* d_v, float* d_f, int sizeX, float h_x, float x_a, float x_b);

Grid1D::Grid1D(int sizeX, float range[]) //range = {xa,xb,ya,yb} xb>xa yb>ya
{
	assert(sizeX > 0);
	assert((sizeX-1)%2 == 0); //size=(2^k) - 1
	this->sizeX = sizeX;

	assert(range[1] > range[0]);
	float x_range = range[1] - range[0];
	x_a = range[0];
	x_b = range[1];
	
	h_x = x_range/(float)(sizeX-1);

	cudaMalloc(&d_v, sizeX*sizeof(float));
	cudaMalloc(&d_f, sizeX*sizeof(float));

	dim3 threadsPerBlock(NUMTHREAD_X); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)sizeX/threadsPerBlock.x);
	dim3 numBlocks(numBlocks_x);

	CUDASetGrids<<<numBlocks, threadsPerBlock>>>(d_v, d_f, sizeX, h_x, x_a, x_b);
}

Grid1D::~Grid1D()
{
	cudaFree(d_v);
	cudaFree(d_f);
}

//ha senso solo per griglia più fine (finest)
void Grid1D::PrintDiffApproxReal(int diff_fd)
{
	float* h_v = (float*)malloc(sizeX*sizeof(float));
	cudaMemcpy(h_v, d_v, sizeX*sizeof(float), cudaMemcpyDeviceToHost);

	int posX;
	for(posX = 0; posX < sizeX; posX++)
	{
		float xj = x_a + posX*h_x;
		float realsol = (exp(xj) + xj - 3)/(1 + exp(-xj));
		float approxsol = h_v[posX];
		float diff = approxsol - realsol;

		char* log = (char*)malloc(100);
		sprintf(log,"xj: %f diff: %f\n", xj, diff);
		write(diff_fd, log, strlen(log));
	}
	free(h_v);
}

void Grid1D::PrintGrid_v(int logfd)
{
	float* h_v = (float*)malloc(sizeX*sizeof(float));
	cudaMemcpy(h_v, d_v, sizeX*sizeof(float), cudaMemcpyDeviceToHost);

	int posX;
	for(posX = 0; posX < sizeX; posX++)
	{
		float xj = x_a + posX*h_x;
		float realsol = (exp(xj) + xj - 3)/(1 + exp(-xj));;
		float approxsol = h_v[posX];
		char* log = (char*)malloc(100);
		sprintf(log,"xj: %f approxSol: %f  realSol: %f\n", xj, approxsol, realsol);
		write(logfd, log, strlen(log));
	}
	free(h_v);	
}

void Grid1D::PrintGrid_f(int logfd)
{
	float* h_f = (float*)malloc(sizeX*sizeof(float));
	cudaMemcpy(h_f, d_f, sizeX*sizeof(float), cudaMemcpyDeviceToHost);

	int posX;
	for(posX = 0; posX < sizeX; posX++)
	{
		float xj = x_a + posX*h_x;

		char* log = (char*)malloc(200);
		sprintf(log, "xj: %f f:%f\n", xj, h_f[posX]);
		write(logfd, log, strlen(log));
	}
	free(h_f);	
}

/*********************************************************CUDA*********************************************************/

__global__ void CUDASetGrids(float* d_v, float* d_f, int sizeX, float h_x, float x_a, float x_b)
{
	int posX = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	
	if(posX >= sizeX)
		return;
	
	float xj = x_a + posX*h_x;
	d_f[posX] = exp(xj);


	if(posX == 0)
	{
		d_v[posX] = (exp(x_a) + x_a - 3)/(1 + exp(-x_a));
		return;
	}
	if(posX == sizeX - 1)
	{
		d_v[posX] = (exp(x_b) + x_b - 3)/(1 + exp(-x_b));
		return;
	}	
}
