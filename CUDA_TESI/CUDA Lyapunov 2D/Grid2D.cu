#include "Grid2D.h"
#include "inclusion.h"

__global__ void CUDASetBoundaries(float* d_v, float* d_f, size_t pitchInByte, int sizeX, int sizeY, float h_x, float h_y, float x_a, float y_a);
//__global__ void CUDACalculateError(float* d_v, int size, size_t pitch, float h_x, float h_y, float x_a, float y_a, float d_);

//OK
Grid2D::Grid2D(int size, float range[]) //range = {xa,xb,ya,yb} xb>xa yb>ya
{
	assert((size-1)%2 == 0); //size=(2^k) - 1
	this->size = size;

	assert(range[1] > range[0]);
	assert(range[3] > range[2]);
	float x_range = range[1] - range[0];
	float y_range = range[3] - range[2];


	x_a = range[0];
	x_b = range[1];
	y_a = range[2];
	y_b = range[3];

	h_x = x_range/(float)(size-1);
	h_y = y_range/(float)(size-1);	

	cudaMallocPitch(&d_v, &d_pitchByte, size*sizeof(float), size);
	cudaMallocPitch(&d_f, &d_pitchByte, size*sizeof(float), size);

	//Azzero i vettori nel Kernel (in questo modo ci sono meno chiamate)
	//cudaMemset2D(d_v, d_pitchByte, 0, size*sizeof(float), size); //azzera vettore v (poi inizializziamo il bordo)
	//cudaMemset2D(d_f, d_pitchByte, 0, size*sizeof(float), size); //azzera vettore f (deve essere tutto azzerato)		
	//Azzero i vettori nel Kernel (in questo modo ci sono meno chiamate)	

	d_pitch = d_pitchByte/sizeof(float);

	//io ho 16 multiprocessori
	dim3 threadsPerBlock(16, 16); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)size/threadsPerBlock.x);	
	int numBlocks_y = ceil((float)size/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);

	CUDASetBoundaries<<<numBlocks, threadsPerBlock>>>(d_v, d_f, d_pitch, size, size, h_x, h_y, x_a, y_a);
}
//OK
Grid2D::~Grid2D()
{
	cudaFree(d_v);
	cudaFree(d_f);
}

//ha senso solo per griglia più fine (finest)
void Grid2D::PrintDiffApproxReal(int diff_fd)
{
	float* h_v = (float*)malloc(size*size*sizeof(float));
	int h_pitchByte = size*sizeof(float); //il pitch di h_v è pari alla dimensione della riga stessa (stà su HOST)
	cudaMemcpy2D(h_v, h_pitchByte, d_v, d_pitchByte, size*sizeof(float), size, cudaMemcpyDeviceToHost);

	int i, j;
	for(i = 0; i < size; i++)
	{
		float* row = h_v + i*size;
		for(j = 0; j < size; j++)		
		{
			float xj = x_a + j*h_x;
			float yi = y_a + i*h_y;	
			float realSol = 2*xj*xj-4*xj*yi+2*yi*yi;
			float approxSol = row[j];
			float diff = realSol - approxSol;

			char* log = (char*)malloc(200);
			sprintf(log,"yi: %f xj: %f v: %f u: %f e: %f\n", yi, xj, approxSol, realSol, diff);
			write(diff_fd, log, strlen(log));
		}	
	}
	free(h_v);
}

void Grid2D::PrintGrid_v(int logfd)
{
	float* h_v = (float*)malloc(size*size*sizeof(float));
	int h_pitchByte = size*sizeof(float); //il pitch di h_v è pari alla dimensione della riga stessa (stà su HOST)
	cudaMemcpy2D(h_v, h_pitchByte, d_v, d_pitchByte, size*sizeof(float), size, cudaMemcpyDeviceToHost);

	int i, j;
	for(i = 0; i < size; i++)
	{
		float* row = h_v + i*size;
		for(j = 0; j < size; j++)		
		{
			float xj = x_a + j*h_x;
			float yi = y_a + i*h_y;	
			float realsol = 2*xj*xj-4*xj*yi+2*yi*yi;
			float approxSol = row[j];
			char* log = (char*)malloc(200);
			sprintf(log, "i: %d  j: %d  approxSol: %f  realSol: %f\n", i, j, approxSol, realsol);
			write(logfd, log, strlen(log));
		}	
	}
	free(h_v);
}

void Grid2D::PrintGrid_f(int logfd)
{
	float* h_f = (float*)malloc(size*size*sizeof(float));
	int h_pitchByte = size*sizeof(float); //il pitch di h_v è pari alla dimensione della riga stessa (stà su HOST)
	cudaMemcpy2D(h_f, h_pitchByte, d_f, d_pitchByte, size*sizeof(float), size, cudaMemcpyDeviceToHost);

	int i, j;
	for(i = 0; i < size; i++)
	{
		float* row = h_f + i*size;
		for(j = 0; j < size; j++)		
		{
			char* log = (char*)malloc(200);
			sprintf(log, "i: %d j: %d f:%f\n", i, j, row[j]);
			write(logfd, log, strlen(log));
		}
	}
	free(h_f);
}

void Grid2D::PrintMeanAbsoluteError()
{
	float* h_v = (float*)malloc(size*size*sizeof(float));
	int h_pitchByte = size*sizeof(float); //il pitch di h_v è pari alla dimensione della riga stessa (stà su HOST)
	cudaMemcpy2D(h_v, h_pitchByte, d_v, d_pitchByte, size*sizeof(float), size, cudaMemcpyDeviceToHost);

	//senza contare i bordi, su cui vi è la soluzione reale
	int numInternalPoints = 0;
	double totalError = 0;

	int i, j;
	for(i = 1; i < size - 1; i++)
	{
		float* row = h_v + i*size;
		for(j = 1; j < size - 1; j++)		
		{
			numInternalPoints++;

			float xj = x_a + j*h_x;
			float yi = y_a + i*h_y;	
			float realsol = 2*xj*xj-4*xj*yi+2*yi*yi;
			float approxSol = row[j];
			float diff = approxSol - realsol;

			totalError += fabs(diff);
		}	
	}
	free(h_v);

	double meanAbsoluteError = totalError/numInternalPoints;
	printf("MeanAbsoluteError: %f\n", meanAbsoluteError);
}

/*********************************************************CUDA*********************************************************/
//imposta valori sul contorno (PITCH NON VA' IN BYTE)
//OK
__global__ void CUDASetBoundaries(float* d_v, float* d_f, size_t pitch, int sizeX, int sizeY, float h_x, float h_y, float x_a, float y_a)
{
	int iy = blockDim.y * blockIdx.y + threadIdx.y; //indice i
	int ix = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	int idx = iy * pitch + ix;		

	if(idx >= pitch*sizeY) //sizeX qui è sostituito dal pitch (NON IN BYTE)
		return;
	
	//d_v[idx] = 0.0f; //azzero valori interi e esterni (ma poi modifico i bordi)
	d_f[idx] = 0.0f; //azzero valori interi e esterni
	
	if(ix == 0 || ix == sizeX - 1 || iy == 0 || iy == sizeY - 1) //stiamo sui bordi
	{
		float yi = y_a + iy*h_y;
		float xj = x_a + ix*h_x;
		float sol = 2*xj*xj-4*xj*yi+2*yi*yi;
		d_v[idx] = sol;	
		return;
	}
	//else
	d_v[idx] = 0.0f; //azzero valori interni
	
}
/*
__global__ void CUDACalculateError(float* d_v, int size, size_t pitch, float h_x, float h_y, float x_a, float y_a, float sumError)
{
	int iy = blockDim.y * blockIdx.y + threadIdx.y; //indice i
	int ix = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	int idx = iy * pitch + ix;		

	if(idx >= pitch*size) //sizeX qui è sostituito dal pitch (NON IN BYTE)
		return;

	float xj = x_a + ix*h_x;
	float yi = y_a + iy*h_y;	
	float realsol = 2*xj*xj-4*xj*yi+2*yi*yi;
	float diff = d_v[idx] - realsol;
	float percentError = (diff/realsol);

	sumError += percentError;
}*/
