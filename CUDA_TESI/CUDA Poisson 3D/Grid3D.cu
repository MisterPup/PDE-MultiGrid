#include "Grid3D.h"
#include "inclusion.h"

__global__ void CUDAInitGrids(float* d_v, float* d_f, int* d_sizeXYZ, float h_x, float h_y, float h_z, float x_a, float y_a, float z_a);

/*
Ogni riga del Kernel è una matrice, quindi sono le righe della matrice affiancate
Scendendo fra le righe mi sposto lungo le z del vettore 3d
*/


//void checkCUDAError(const char *msg);

Grid3D::Grid3D(int sizeXYZ[], float range[]) //range = {xa,xb,ya,yb,za,zb} xb>xa yb>ya zb>za
{
	int sizeX = sizeXYZ[0];
	int sizeY = sizeXYZ[1];
	int sizeZ = sizeXYZ[2];

	assert(sizeX == sizeY); //per comodità
	assert(sizeX == sizeZ);

	assert((sizeX-1)%2 == 0); //sizeX=(2^k) - 1
	this->sizeX = sizeX;

	assert((sizeY-1)%2 == 0); //sizeY=(2^k) - 1
	this->sizeY = sizeY;

	assert((sizeZ-1)%2 == 0); //sizeZ=(2^k) - 1
	this->sizeZ = sizeZ;
	
	cudaMalloc(&d_sizeXYZ, 3*sizeof(int));

	cudaMemcpy(d_sizeXYZ, sizeXYZ, 3*sizeof(int), cudaMemcpyHostToDevice);
	//checkCUDAError("cudaMemcpy1");
	cudaMemcpy(sizeXYZ, d_sizeXYZ,  3*sizeof(int), cudaMemcpyDeviceToHost);
	//checkCUDAError("cudaMemcpy2");
	assert(range[1] > range[0]);
	assert(range[3] > range[2]);
	assert(range[5] > range[4]);

	float x_range = range[1] - range[0];
	float y_range = range[3] - range[2];
	float z_range = range[5] - range[4];

	x_a = range[0];
	x_b = range[1];
	y_a = range[2];
	y_b = range[3];
	z_a = range[4];
	z_b = range[5];

	h_x = x_range/(float)(sizeX-1);
	h_y = y_range/(float)(sizeY-1);
	h_z = z_range/(float)(sizeZ-1);

	cudaMalloc(&d_v, sizeX*sizeY*sizeZ*sizeof(float));
	cudaMalloc(&d_f, sizeX*sizeY*sizeZ*sizeof(float));	
	//checkCUDAError("cudaMalloc");

	//InitF();

	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(sizeX*sizeY)/threadsPerBlock.x);	
	int numBlocks_y = ceil((float)sizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);
	//printf("numBlocks_x: %d\n", numBlocks_x);
	//printf("numBlocks_y: %d\n", numBlocks_y);	
	
	CUDAInitGrids<<<numBlocks, threadsPerBlock>>>(d_v, d_f, d_sizeXYZ, h_x, h_y, h_z, x_a, y_a, z_a);
	//checkCUDAError("kernel execution");	
}

Grid3D::~Grid3D()
{
	cudaFree(d_v);
	cudaFree(d_f);
	cudaFree(d_sizeXYZ);
}
/*
void Grid3D::InitF()
{
	float* h_f = (float*)malloc(sizeX*sizeY*sizeZ*sizeof(float));
	int posY, posX, posZ;
	for(posY = 0; posY < sizeY; posY++)
	{
		for(posX = 0; posX < sizeX; posX++)
		{
			for(posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				float x = x_a + posX*h_x;
				float y = y_a + posY*h_y;
				float z = z_a + posZ*h_z;
				
				//per far effettuare i calcoli in parallelo alle cuda devi usare sinf (và bene anche sin)
				h_f[idx] = -3*PI_H*PI_H*sin(PI_H*x)*sin(PI_H*y)*sin(PI_H*z);
			}
		}
	}
	cudaMemcpy(d_f, h_f, sizeX*sizeY*sizeZ*sizeof(float), cudaMemcpyHostToDevice);
	free(h_f);
}*/

//ha senso solo per griglia più fine (finest)
void Grid3D::PrintDiffApproxReal(int diff_fd)
{
	float* h_v = (float*)malloc(sizeX*sizeY*sizeZ*sizeof(float));
	cudaMemcpy(h_v, d_v, sizeX*sizeY*sizeZ*sizeof(float), cudaMemcpyDeviceToHost);

	int posY, posX, posZ;
	for(posY = 0; posY < sizeY; posY++)
	{
		for(posX = 0; posX < sizeX; posX++)
		{
			for(posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				float x = x_a + posX*h_x;
				float y = y_a + posY*h_y;
				float z = z_a + posZ*h_z;		

				float realSol = sin(PI_H*x)*sin(PI_H*y)*sin(PI_H*z);
				float approxSol = h_v[idx];
				float diff = realSol - approxSol;
				char* log = (char*)malloc(200);
				sprintf(log,"posY: %d posX: %d posZ: %d diff: %f\n", posY, posX, posZ, diff);
				write(diff_fd, log, strlen(log));
			}
		}
	}
	free(h_v);
}

void Grid3D::PrintGrid_v(int logfd)
{
	float* h_v = (float*)malloc(sizeX*sizeY*sizeZ*sizeof(float));
	cudaMemcpy(h_v, d_v, sizeX*sizeY*sizeZ*sizeof(float), cudaMemcpyDeviceToHost);

	int posY, posX, posZ;
	for(posY = 0; posY < sizeY; posY++)
	{
		for(posX = 0; posX < sizeX; posX++)
		{
			for(posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				float x = x_a + posX*h_x;
				float y = y_a + posY*h_y;
				float z = z_a + posZ*h_z;		

				float realSol = sin(PI_H*x)*sin(PI_H*y)*sin(PI_H*z);
				float approxSol = h_v[idx];

				char* log = (char*)malloc(200);
				sprintf(log,"posY: %d posX: %d posZ: %d value: %f\n", posY, posX, posZ, approxSol);
				//sprintf(log,"posY: %d posX: %d posZ: %d approxSol: %f, realSol: %f\n", posY, posX, posZ, approxSol, realSol);
				write(logfd, log, strlen(log));
				free(log);
			}
		}
	}
	free(h_v);
}

void Grid3D::PrintGrid_f(int logfd)
{
	float* h_f = (float*)malloc(sizeX*sizeY*sizeZ*sizeof(float));
	cudaMemcpy(h_f, d_f, sizeX*sizeY*sizeZ*sizeof(float), cudaMemcpyDeviceToHost);

	int posY, posX, posZ;
	for(posY = 0; posY < sizeY; posY++)
	{
		for(posX = 0; posX < sizeX; posX++)
		{
			for(posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				char* log = (char*)malloc(200);
				sprintf(log,"posY: %d posX: %d posZ: %d value: %f\n", posY, posX, posZ, h_f[idx]);
				write(logfd, log, strlen(log));
			}
		}
	}
	free(h_f);
}

/*********************************************************CUDA*********************************************************/

__global__ void CUDAInitGrids(float* d_v, float* d_f, int* d_sizeXYZ, float h_x, float h_y, float h_z, float x_a, float y_a, float z_a)
{
	int sizeX = d_sizeXYZ[0];
	int sizeY = d_sizeXYZ[1];
	int sizeZ = d_sizeXYZ[2];

	float PI_D = CUDART_PI_F; //PI GRECO SU CUDA (precisione singola)

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int posY = (threadIdx.x + blockIdx.x*blockDim.x)/sizeX; //il numeratore è l'indice assoluto lungo le x
	int posX = (threadIdx.x + blockIdx.x*blockDim.x) - posY*sizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int posZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

	if(posY >= sizeY || posX >= sizeX || posZ >= sizeZ)
		return;

	float x = x_a + posX*h_x;
	float y = y_a + posY*h_y;
	float z = z_a + posZ*h_z;
	
	d_f[idx] = -3*PI_D*PI_D*sin(PI_D*x)*sin(PI_D*y)*sin(PI_D*z);		

	if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1 || posZ == 0 || posZ == sizeZ - 1)
	{
		d_v[idx] = 0.0f; //soluzione nulla lungo il contorno	
		return;
	}
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
