#include "Grid3D.h"
#include "MultiGrid3D.h"
#include "inclusion.h"

__global__ void CUDARestrict(float* fine, int fsizeXYZ[], float* coarse, int csizeXYZ[]);
__global__ void CUDAInterpolate(float* fine, int fsizeXYZ[], float* coarse, int csizeXYZ[]);

__global__ void CUDARelax(float* d_v, float* d_f, float h_x, float h_y, float h_z, int d_sizeXYZ[]); 
__global__ void CUDACalculateResidual(float* d_v, float* d_f, float* d_residual, float h_x, float h_y, float h_z, int d_sizeXYZ[]);
__global__ void CUDAApplyCorrection(float* fine, float* error, int d_sizeXYZ[]);
 

__global__ void CUDASet(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder);
__global__ void CUDASetTESTTEST(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder);





void checkCUDAError(const char *msg);

MultiGrid3D::MultiGrid3D(int finestGridSizeXYZ[], float range[])
{
	InitGrids(finestGridSizeXYZ, range);
}

MultiGrid3D::~MultiGrid3D()
{
	for(int i = 0; i < numGrids; i++)
	{
		grids3D[i]->~Grid3D();
	}
	free(grids3D);
}

void MultiGrid3D::InitGrids(int finestGridSizeXYZ[], float range[])
{
	int finestGridSizeX = finestGridSizeXYZ[0];
	int finestGridSizeY = finestGridSizeXYZ[1];
	int finestGridSizeZ = finestGridSizeXYZ[2];

	int minSize = finestGridSizeX;
	if(finestGridSizeY < minSize)
		minSize = finestGridSizeY;

	if(finestGridSizeZ < minSize)
		minSize = finestGridSizeZ;


	double tempNumGrids = minSize - 1;
	numGrids = log2(tempNumGrids); //finestGrid size = 65 = (2^6)+1 ---> numGrids = 6 (non creo griglie senza punti interni)

	grids3D = (Grid3D**)malloc(numGrids*sizeof(Grid3D*));
	grids3D[0] = new Grid3D(finestGridSizeXYZ, range);
	for(int i = 1; i < numGrids; i++)
	{
		int coarseGridSizeX = ((grids3D[i-1]->sizeX-1)/2)+1; //coarseGridSize = ((2^/finerGridSize)-1)/2)+1
		int coarseGridSizeY = ((grids3D[i-1]->sizeY-1)/2)+1;
		int coarseGridSizeZ = ((grids3D[i-1]->sizeZ-1)/2)+1;

		int coarseGridSizeXYZ[3] = {coarseGridSizeX, coarseGridSizeY, coarseGridSizeZ};
		grids3D[i] = new Grid3D(coarseGridSizeXYZ, range);
	}
}

void MultiGrid3D::Restrict(float* d_fine, int d_fsizeXYZ[], float* d_coarse, int d_csizeXYZ[])
{
	int* h_fsizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_fsizeXYZ, d_fsizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);

	int* h_csizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_csizeXYZ, d_csizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);


	int fsizeX = h_fsizeXYZ[0];
	int fsizeY = h_fsizeXYZ[1];
	int fsizeZ = h_fsizeXYZ[2];

	int csizeX = h_csizeXYZ[0];
	int csizeY = h_csizeXYZ[1];
	int csizeZ = h_csizeXYZ[2];

	free(h_fsizeXYZ);
	free(h_csizeXYZ);

	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni
	assert(csizeY == ((fsizeY-1)/2+1));
	assert(csizeZ == ((fsizeZ-1)/2+1));

	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(csizeX*csizeY)/threadsPerBlock.x); //csize perchè in Restrict riempo una griglia più rada (quindi meno punti)	
	int numBlocks_y = ceil((float)csizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);

	CUDARestrict<<<numBlocks, threadsPerBlock>>>(d_fine, d_fsizeXYZ, d_coarse, d_csizeXYZ);
}



void MultiGrid3D::SetTESTTEST(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder)
{
	int* h_sizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_sizeXYZ, d_sizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);

	int sizeX = h_sizeXYZ[0];
	int sizeY = h_sizeXYZ[1];
	int sizeZ = h_sizeXYZ[2];
	free(h_sizeXYZ);
	
	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(sizeX*sizeY)/threadsPerBlock.x);	
	int numBlocks_y = ceil((float)sizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);

	CUDASetTESTTEST<<<numBlocks, threadsPerBlock>>>(d_v, d_sizeXYZ, value, modifyBorder);
}


void MultiGrid3D::Set(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder)
{
	int* h_sizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_sizeXYZ, d_sizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);

	int sizeX = h_sizeXYZ[0];
	int sizeY = h_sizeXYZ[1];
	int sizeZ = h_sizeXYZ[2];
	free(h_sizeXYZ);
	
	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(sizeX*sizeY)/threadsPerBlock.x);	
	int numBlocks_y = ceil((float)sizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);

	CUDASet<<<numBlocks, threadsPerBlock>>>(d_v, d_sizeXYZ, value, modifyBorder);
}


void MultiGrid3D::Interpolate(float* d_fine, int d_fsizeXYZ[], float* d_coarse, int d_csizeXYZ[])
{
	int* h_fsizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_fsizeXYZ, d_fsizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);

	int* h_csizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_csizeXYZ, d_csizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);


	int fsizeX = h_fsizeXYZ[0];
	int fsizeY = h_fsizeXYZ[1];
	int fsizeZ = h_fsizeXYZ[2];

	int csizeX = h_csizeXYZ[0];
	int csizeY = h_csizeXYZ[1];
	int csizeZ = h_csizeXYZ[2];

	free(h_fsizeXYZ);
	free(h_csizeXYZ);

	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni
	assert(csizeY == ((fsizeY-1)/2+1));
	assert(csizeZ == ((fsizeZ-1)/2+1));

	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(fsizeX*fsizeY)/threadsPerBlock.x); //fsize perchè in Interpolate riempo una griglia più densa (quindi più punti)	
	int numBlocks_y = ceil((float)fsizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);
	
	CUDAInterpolate<<<numBlocks, threadsPerBlock>>>(d_fine, d_fsizeXYZ, d_coarse, d_csizeXYZ);
}

void MultiGrid3D::Relax(Grid3D* curGrid, int ncycles)
{
	float* d_v = curGrid->d_v;
	float* d_f = curGrid->d_f;

	float h_x = curGrid->h_x;
	float h_y = curGrid->h_y;
	float h_z = curGrid->h_z;

	int* h_sizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_sizeXYZ, curGrid->d_sizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);

	int sizeX = h_sizeXYZ[0];
	int sizeY = h_sizeXYZ[1];
	int sizeZ = h_sizeXYZ[2];

	free(h_sizeXYZ);

	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(sizeX*sizeY)/threadsPerBlock.x);	
	int numBlocks_y = ceil((float)sizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);

	for(int k = 0; k < ncycles; k++)                
		CUDARelax<<<numBlocks, threadsPerBlock>>>(d_v, d_f, h_x, h_y, h_z, curGrid->d_sizeXYZ);
}

float* MultiGrid3D::CalculateResidual(Grid3D* curGrid)
{
	float* d_v = curGrid->d_v;
	float* d_f = curGrid->d_f;

	float h_x = curGrid->h_x;
	float h_y = curGrid->h_y;
	float h_z = curGrid->h_z;

	int* h_sizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_sizeXYZ, curGrid->d_sizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);

	int sizeX = h_sizeXYZ[0];
	int sizeY = h_sizeXYZ[1];
	int sizeZ = h_sizeXYZ[2];

	free(h_sizeXYZ);

	float* d_residual; //residuo (DEVICE)

	cudaMalloc(&d_residual, sizeX*sizeY*sizeZ*sizeof(float));

	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(sizeX*sizeY)/threadsPerBlock.x);	
	int numBlocks_y = ceil((float)sizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);
							     
	CUDACalculateResidual<<<numBlocks, threadsPerBlock>>>(d_v, d_f, d_residual, h_x, h_y, h_z, curGrid->d_sizeXYZ);
	checkCUDAError("Residual");
	return d_residual;
}


void MultiGrid3D::ApplyCorrection(float* fine, int d_fsizeXYZ[], float* error, int d_esizeXYZ[])
{
	int* h_fsizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_fsizeXYZ, d_fsizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);

	int* h_esizeXYZ = (int*)malloc(3*sizeof(int));
	cudaMemcpy(h_esizeXYZ, d_esizeXYZ, 3*sizeof(int), cudaMemcpyDeviceToHost);


	int fsizeX = h_fsizeXYZ[0];
	int fsizeY = h_fsizeXYZ[1];
	int fsizeZ = h_fsizeXYZ[2];

	int esizeX = h_esizeXYZ[0];
	int esizeY = h_esizeXYZ[1];
	int esizeZ = h_esizeXYZ[2];

	free(h_fsizeXYZ);
	free(h_esizeXYZ);

	//potresti anche lavorare sui bordi, tanto non vengono modificati perchè l'errore è sicuramente nullo sugli stessi
	assert(fsizeX == esizeX);
	assert(fsizeY == esizeY);
	assert(fsizeZ == esizeZ);

	//4 (MP) x 48 (Cores/MP) = 192 (Cores)		
	dim3 threadsPerBlock(NUMTHREAD_X, NUMTHREAD_Y); //il massimo numero di thread attivi per multiprocessore è 768, ma ogni blocco può avere al max 512 thread 
	int numBlocks_x = ceil((float)(fsizeX*fsizeY)/threadsPerBlock.x);	
	int numBlocks_y = ceil((float)fsizeZ/threadsPerBlock.y);
	dim3 numBlocks(numBlocks_x, numBlocks_y);

	CUDAApplyCorrection<<<numBlocks, threadsPerBlock>>>(fine, error, d_fsizeXYZ);
}

void MultiGrid3D::VCycle(int gridID, int v1, int v2)
{
	Grid3D* fine = grids3D[gridID];
	Relax(fine, v1);
	//printf("Relax\n");
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		float* residual = CalculateResidual(fine);
		int* residualSize = fine->d_sizeXYZ;
		Grid3D* coarse = grids3D[gridID+1];
		Restrict(residual, residualSize, coarse->d_f, coarse->d_sizeXYZ);
		//printf("Restrict\n");	
		Set(coarse->d_v, coarse->d_sizeXYZ, 0.0f, true); //azzero tutto il vettore di errori (sui bordi sicuramente e=0)
		//printf("Set\n");
		VCycle(gridID+1, v1, v2);
		//printf("Salgo\n");
		
		int* fsize = fine->d_sizeXYZ;

		float* d_fine_error;
		cudaMalloc(&d_fine_error, fine->sizeX*fine->sizeY*fine->sizeZ*sizeof(float)); //errore sulla griglia fine
		//coarse->v contiene errore su griglia rada		
		Interpolate(d_fine_error, fsize, coarse->d_v, coarse->d_sizeXYZ); //porti l'errore da lvl kh a lvl (k-1)h, nullo lungo bordi
		//printf("Interpolate\n");
		ApplyCorrection(fine->d_v, fsize, d_fine_error, fsize);
		//printf("ApplyCorrection\n");
	}

	Relax(fine, v2);
	//printf("Relax\n");
}

void MultiGrid3D::FullMultiGridVCycle(int gridID, int v0, int v1, int v2)
{
	Grid3D* fine = grids3D[gridID];
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		Grid3D* coarse = grids3D[gridID + 1];
		Restrict(fine->d_f, fine->d_sizeXYZ, coarse->d_f, coarse->d_sizeXYZ);
		FullMultiGridVCycle(gridID + 1,v0, v1, v2);
		Interpolate(fine->d_v, fine->d_sizeXYZ, coarse->d_v, coarse->d_sizeXYZ);
	}
	else //i contorni non vanno modificati in quanto nel FMG lavoriamo sulla sol approx e non sull'errore	   
		Set(fine->d_v, fine->d_sizeXYZ, 0.0f, false);
	
	for(int i = 0; i < v0; i++)
		VCycle(gridID, v1, v2);
}

void MultiGrid3D::PrintDiff()
{
	int diff_fd = open("log/diff.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids3D[0]->PrintDiffApproxReal(diff_fd);
	close(diff_fd);
}

void MultiGrid3D::PrintGrid(int gridID)
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids3D[gridID]->PrintGrid_v(logfd);
	close(logfd);
}

void MultiGrid3D::PrintAllGrids_v()
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		printf("Grid: %d\n", i);
		grids3D[i]->PrintGrid_v(logfd);
	}

	close(logfd);
}

void MultiGrid3D::PrintAllGrids_f()
{
	int logfd = open("log/log_f.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		printf("Grid: %d\n", i);
		grids3D[i]->PrintGrid_f(logfd);

	}

	close(logfd);
}

/*********************************************************CUDA*********************************************************/

__global__ void CUDARestrict(float* fine, int d_fsizeXYZ[], float* coarse, int d_csizeXYZ[])
{
	int fsizeX = d_fsizeXYZ[0];
	int fsizeY = d_fsizeXYZ[1];
	//int fsizeZ = d_fsizeXYZ[2];

	int csizeX = d_csizeXYZ[0];
	int csizeY = d_csizeXYZ[1];
	int csizeZ = d_csizeXYZ[2];

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int cposY = (threadIdx.x + blockIdx.x*blockDim.x)/csizeX; //il numeratore è l'indice assoluto lungo le x
	int cposX = (threadIdx.x + blockIdx.x*blockDim.x) - cposY*csizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int cposZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY; //indice elemento griglia rada

	int fposX = 2*cposX;
	int fposY = 2*cposY;
	int fposZ = 2*cposZ;
	int fidx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY; //indice elemento griglia fine

	if(cposY >= csizeY || cposX >= csizeX || cposZ >= csizeZ)
		return;

	/* COPIO VALORI SUI BORDI */
	if(cposX == 0 || cposX == csizeX - 1 || cposY == 0 || cposY == csizeY - 1 || cposZ == 0 || cposZ == csizeZ - 1)
	{
		coarse[cidx] = fine[fidx];
		return;
	}
	/* COPIO VALORI SUI BORDI */
	
	//else non stiamo sui bordi
	
	fidx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY;
	float C_C = fine[fidx];
	fidx = fposX + fposY * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float N_C = fine[fidx];
	fidx = fposX + fposY * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float S_C = fine[fidx];
	fidx = (fposX+1) + fposY * fsizeX + fposZ * fsizeX * fsizeY;
	float E_C = fine[fidx];
	fidx = (fposX-1) + fposY * fsizeX + fposZ * fsizeX * fsizeY;
	float O_C = fine[fidx];
	fidx = (fposX+1) + fposY * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float NE_C = fine[fidx];
	fidx = (fposX-1) + fposY * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float NO_C = fine[fidx];
	fidx = (fposX+1) + fposY * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float SE_C = fine[fidx];
	fidx = (fposX-1) + fposY * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float SO_C = fine[fidx];

	fidx = fposX + (fposY-1) * fsizeX + fposZ * fsizeX * fsizeY;
	float C_N = fine[fidx];
	fidx = fposX + (fposY-1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float N_N = fine[fidx];
	fidx = fposX + (fposY-1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float S_N = fine[fidx];
	fidx = (fposX+1) + (fposY-1) * fsizeX + (fposZ) * fsizeX * fsizeY;
	float E_N = fine[fidx];
	fidx = (fposX-1) + (fposY-1) * fsizeX + (fposZ) * fsizeX * fsizeY;
	float O_N = fine[fidx];
	fidx = (fposX+1) + (fposY-1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float NE_N = fine[fidx];
	fidx = (fposX-1) + (fposY-1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float NO_N = fine[fidx];
	fidx = (fposX+1) + (fposY-1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float SE_N = fine[fidx];
	fidx = (fposX-1) + (fposY-1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float SO_N = fine[fidx];

	fidx = fposX + (fposY+1) * fsizeX + fposZ * fsizeX * fsizeY;
	float C_S = fine[fidx];
	fidx = fposX + (fposY+1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float N_S = fine[fidx];
	fidx = fposX + (fposY+1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float S_S = fine[fidx];
	fidx = (fposX+1) + (fposY+1) * fsizeX + (fposZ) * fsizeX * fsizeY;
	float E_S = fine[fidx];
	fidx = (fposX-1) + (fposY+1) * fsizeX + (fposZ) * fsizeX * fsizeY;
	float O_S = fine[fidx];
	fidx = (fposX+1) + (fposY+1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float NE_S = fine[fidx];
	fidx = (fposX-1) + (fposY+1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
	float NO_S = fine[fidx];
	fidx = (fposX+1) + (fposY+1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float SE_S = fine[fidx];
	fidx = (fposX-1) + (fposY+1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
	float SO_S = fine[fidx];

	cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
	coarse[cidx] = (1/8.0f)*(C_C) + (1/16.0f)*((N_C+E_C+S_C+O_C) + (C_N+C_S)) + (1/32.0f)*((NE_C+SE_C+SO_C+NO_C) + (N_N+E_N+S_N+O_N) + (N_S+E_S+S_S+O_S)) + (1/64.0f)*((NE_N+SE_N+SO_N+NO_N) + (NE_S+SE_S+SO_S+NO_S));

}

__global__ void CUDAInterpolate(float* fine, int d_fsizeXYZ[], float* coarse, int d_csizeXYZ[])
{
	int fsizeX = d_fsizeXYZ[0];
	int fsizeY = d_fsizeXYZ[1];
	int fsizeZ = d_fsizeXYZ[2];

	int csizeX = d_csizeXYZ[0];
	int csizeY = d_csizeXYZ[1];
	//int csizeZ = d_csizeXYZ[2];

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int fposY = (threadIdx.x + blockIdx.x*blockDim.x)/fsizeX; //il numeratore è l'indice assoluto lungo le x
	int fposX = (threadIdx.x + blockIdx.x*blockDim.x) - fposY*fsizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int fposZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int fidx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY; //indice elemento griglia rada

	int cposX = fposX/2; //divisione fra interi: 3/2 = 1!
	int cposY = fposY/2;
	int cposZ = fposZ/2;
	int cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY; //indice elemento griglia fine

	if(fposY >= fsizeY || fposX >= fsizeX || fposZ >= fsizeZ)
		return;

	//i bordi non vanno modificati
	if(fposX == 0 || fposX == fsizeX - 1 || fposY == 0 || fposY == fsizeY - 1 || fposZ == 0 || fposZ == fsizeZ - 1)
		return; 

	if(fposY%2 == 0 && fposX%2 == 0 && fposZ%2 == 0) //PPP
	{
		fine[fidx] = coarse[cidx];
		return;
	}

	if(fposY%2 == 0 && fposX%2 != 0 && fposZ%2 == 0) //PDP
	{
		cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
		float O = coarse[cidx];
		cidx = (cposX+1) + cposY * csizeX + cposZ * csizeX * csizeY;
		float E = coarse[cidx];

		fine[fidx] = (1/2.0f)*(O + E);
		return;
	}

	if(fposY%2 != 0 && fposX%2 == 0 && fposZ%2 == 0) //DPP
	{
		cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
		float N = coarse[cidx];
		cidx = cposX + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
		float S = coarse[cidx];

		fine[fidx] = (1/2.0f)*(N + S);
		return;
	}

	if(fposY%2 != 0 && fposX%2 != 0 && fposZ%2 == 0) //DDP
	{
		cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
		float NO = coarse[cidx];
		cidx = (cposX+1) + cposY * csizeX + cposZ * csizeX * csizeY;
		float NE = coarse[cidx];
		cidx = cposX + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
		float SO = coarse[cidx];
		cidx = (cposX+1) + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
		float SE = coarse[cidx];

		fine[fidx] = (1/4.0f)*(NO + NE + SO + SE);
		return;
	}

	/* *******************************************************/

	if(fposY%2 == 0 && fposX%2 == 0 && fposZ%2 != 0) //PPD
	{
		cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
		float S = coarse[cidx];
		cidx = cposX + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
		float N = coarse[cidx];

		fine[fidx] = (1/2.0f)*(S + N);
		return;
	}

	if(fposY%2 == 0 && fposX%2 != 0 && fposZ%2 != 0) //PDD
	{
		cidx = cposX + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
		float NO = coarse[cidx];
		cidx = (cposX+1) + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
		float NE = coarse[cidx];
		cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
		float SO = coarse[cidx];
		cidx = (cposX+1) + cposY * csizeX + cposZ * csizeX * csizeY;
		float SE = coarse[cidx];

		fine[fidx] = (1/4.0f)*(NO + NE + SO + SE);
		return;
	}

	if(fposY%2 != 0 && fposX%2 == 0 && fposZ%2 != 0) //DPD
	{
		cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
		float NO = coarse[cidx];
		cidx = cposX + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
		float NE = coarse[cidx];
		cidx = cposX + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
		float SO = coarse[cidx];
		cidx = cposX + (cposY+1) * csizeX + (cposZ+1) * csizeX * csizeY;
		float SE = coarse[cidx];

		fine[fidx] = (1/4.0f)*(NO + NE + SO + SE);
		return;
	}

	if(fposY%2 != 0 && fposX%2 != 0 && fposZ%2 != 0) //DDD
	{
		cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
		float USO = coarse[cidx]; //UP SUD OVEST

		cidx = cposX + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
		float UNO = coarse[cidx];

		cidx = (cposX+1) + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
		float UNE = coarse[cidx];

		cidx = (cposX+1) + cposY * csizeX + cposZ * csizeX * csizeY;
		float USE = coarse[cidx];


		cidx = cposX + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
		float DSO = coarse[cidx]; //DOWN SUD OVEST

		cidx = cposX + (cposY+1) * csizeX + (cposZ+1) * csizeX * csizeY;
		float DNO = coarse[cidx];

		cidx = (cposX+1) + (cposY+1) * csizeX + (cposZ+1) * csizeX * csizeY;
		float DNE = coarse[cidx];

		cidx = (cposX+1) + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
		float DSE = coarse[cidx];

		fine[fidx] = (1/8.0f)*(USO + UNO + UNE + USE  +  DSO + DNO + DNE + DSE);
		
	}

}

__global__ void CUDARelax(float* d_v, float* d_f, float h_x, float h_y, float h_z, int d_sizeXYZ[])
{
	int sizeX = d_sizeXYZ[0];
	int sizeY = d_sizeXYZ[1];
	int sizeZ = d_sizeXYZ[2];

	float h_x2 = h_x*h_x;
	float h_y2 = h_y*h_y;
	float h_z2 = h_z*h_z;

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int posY = (threadIdx.x + blockIdx.x*blockDim.x)/sizeX; //il numeratore è l'indice assoluto lungo le x
	int posX = (threadIdx.x + blockIdx.x*blockDim.x) - posY*sizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int posZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int idx = posX + posY * sizeX + posZ * sizeX * sizeY; //indice elemento griglia rada


	if(posY >= sizeY || posX >= sizeX || posZ >= sizeZ)
		return;

	//i bordi non vanno modificati
	if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1 || posZ == 0 || posZ == sizeZ - 1)
		return; 

	if((posY+ posX + posZ)%2 == 0) //punti pari = punti rossi
	{
		idx = (posX-1) + posY * sizeX + posZ * sizeX * sizeY;
		float O = d_v[idx];
		idx = (posX+1) + posY * sizeX + posZ * sizeX * sizeY;
		float E = d_v[idx];
		idx = posX + (posY-1) * sizeX + posZ * sizeX * sizeY;
		float N = d_v[idx];
		idx = posX + (posY+1) * sizeX + posZ * sizeX * sizeY;
		float S = d_v[idx];
		idx = posX + posY * sizeX + (posZ-1) * sizeX * sizeY;
		float D = d_v[idx];
		idx = posX + posY * sizeX + (posZ+1) * sizeX * sizeY;
		float U = d_v[idx];

		idx = posX + posY * sizeX + posZ * sizeX * sizeY;
		d_v[idx] = (O*(h_y2*h_z2)+E*(h_y2*h_z2) + N*(h_x2*h_z2)+S*(h_x2*h_z2) + D*(h_x2*h_y2)+U*(h_x2*h_y2) - d_f[idx]*h_x2*h_y2*h_z2)/(2*(h_y2*h_z2 + h_x2*h_z2 + h_x2*h_y2));	
	}

	__syncthreads();//prima di aggiornare i punti neri(che dipendono esclusivamente da quelli rossi), bisogna aspettare che ogni punto rosso sia stato aggiornato
	
	if((posY+ posX + posZ)%2 != 0) //punti dispari = punti neri
	{
		idx = (posX-1) + posY * sizeX + posZ * sizeX * sizeY;
		float O = d_v[idx];
		idx = (posX+1) + posY * sizeX + posZ * sizeX * sizeY;
		float E = d_v[idx];
		idx = posX + (posY-1) * sizeX + posZ * sizeX * sizeY;
		float N = d_v[idx];
		idx = posX + (posY+1) * sizeX + posZ * sizeX * sizeY;
		float S = d_v[idx];
		idx = posX + posY * sizeX + (posZ-1) * sizeX * sizeY;
		float D = d_v[idx];
		idx = posX + posY * sizeX + (posZ+1) * sizeX * sizeY;
		float U = d_v[idx];

		idx = posX + posY * sizeX + posZ * sizeX * sizeY;
		d_v[idx] = (O*(h_y2*h_z2)+E*(h_y2*h_z2) + N*(h_x2*h_z2)+S*(h_x2*h_z2) + D*(h_x2*h_y2)+U*(h_x2*h_y2) - d_f[idx]*h_x2*h_y2*h_z2)/(2*(h_y2*h_z2 + h_x2*h_z2 + h_x2*h_y2));
	}

}

__global__ void CUDASet(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder)
{
	int sizeX = d_sizeXYZ[0];
	int sizeY = d_sizeXYZ[1];
	int sizeZ = d_sizeXYZ[2];

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int posY = (threadIdx.x + blockIdx.x*blockDim.x)/sizeX; //il numeratore è l'indice assoluto lungo le x
	int posX = (threadIdx.x + blockIdx.x*blockDim.x) - posY*sizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int posZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int idx = posX + posY * sizeX + posZ * sizeX * sizeY; //indice elemento griglia rada

	if(posY >= sizeY || posX >= sizeX || posZ >= sizeZ)
		return;

	if(!modifyBorder) //non devo modificare i bordi
		if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1 || posZ == 0 || posZ == sizeZ - 1)	
			return;
	
	//else
	d_v[idx] = value;
}

__global__ void CUDASetTESTTEST(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder)
{
	int sizeX = d_sizeXYZ[0];
	int sizeY = d_sizeXYZ[1];
	int sizeZ = d_sizeXYZ[2];

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int posY = (threadIdx.x + blockIdx.x*blockDim.x)/sizeX; //il numeratore è l'indice assoluto lungo le x
	int posX = (threadIdx.x + blockIdx.x*blockDim.x) - posY*sizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int posZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int idx = posX + posY * sizeX + posZ * sizeX * sizeY; //indice elemento griglia rada

	if(posY >= sizeY || posX >= sizeX || posZ >= sizeZ)
		return;

	d_v[idx] = posX + posY + posZ;
}


__global__ void CUDACalculateResidual(float* d_v, float* d_f, float* d_residual, float h_x, float h_y, float h_z, int d_sizeXYZ[])
{
	int sizeX = d_sizeXYZ[0];
	int sizeY = d_sizeXYZ[1];
	int sizeZ = d_sizeXYZ[2];

	float h_x2 = h_x*h_x;
	float h_y2 = h_y*h_y;
	float h_z2 = h_z*h_z;

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int posY = (threadIdx.x + blockIdx.x*blockDim.x)/sizeX; //il numeratore è l'indice assoluto lungo le x
	int posX = (threadIdx.x + blockIdx.x*blockDim.x) - posY*sizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int posZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int idx = posX + posY * sizeX + posZ * sizeX * sizeY; //indice elemento griglia rada


	if(posY >= sizeY || posX >= sizeX || posZ >= sizeZ)
		return;

	//residuo sui bordi
	if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1 || posZ == 0 || posZ == sizeZ - 1)
	{
		d_residual[idx] = 0.0; //residui sui bordi nullo
		return;
	}

	//else
	idx = (posX-1) + posY * sizeX + posZ * sizeX * sizeY;
	float O = d_v[idx];
	idx = (posX+1) + posY * sizeX + posZ * sizeX * sizeY;
	float E = d_v[idx];
	idx = posX + (posY-1) * sizeX + posZ * sizeX * sizeY;
	float N = d_v[idx];
	idx = posX + (posY+1) * sizeX + posZ * sizeX * sizeY;
	float S = d_v[idx];
	idx = posX + posY * sizeX + (posZ-1) * sizeX * sizeY;
	float D = d_v[idx];
	idx = posX + posY * sizeX + (posZ+1) * sizeX * sizeY;
	float U = d_v[idx];

	idx = posX + posY * sizeX + posZ * sizeX * sizeY;
	d_residual[idx] = d_f[idx] - ((O-2*d_v[idx]+E)/h_x2) - ((N-2*d_v[idx]-S)/h_y2) - ((D-2*d_v[idx]-U)/h_z2);
}
     
__global__ void CUDAApplyCorrection(float* fine, float* error, int d_sizeXYZ[])
{
	int sizeX = d_sizeXYZ[0];
	int sizeY = d_sizeXYZ[1];
	int sizeZ = d_sizeXYZ[2];

	//ogni riga di thread del kernel(riga divisa in più blocchi di thread), costituisce una matrice sizeX*sizeY
	//spostandomi lungo le y del Kernel, mi sposto lungo le z del vettore 3d
	int posY = (threadIdx.x + blockIdx.x*blockDim.x)/sizeX; //il numeratore è l'indice assoluto lungo le x
	int posX = (threadIdx.x + blockIdx.x*blockDim.x) - posY*sizeX; //indice assoluto lungo le x - numero di righe che precedono l'elemento
	int posZ = (blockIdx.y * blockDim.y) + threadIdx.y;

	int idx = posX + posY * sizeX + posZ * sizeX * sizeY; //indice elemento griglia rada


	if(posY >= sizeY || posX >= sizeX || posZ >= sizeZ)
		return;

	//i bordi non vanno modificati! Lì la soluzione è già quella reale
	if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1 || posZ == 0 || posZ == sizeZ - 1)
		return;
	
	fine[idx] = fine[idx] + error[idx]; //u = v + e (u potrebbe essere usato come errore nella griglia un pò più fine)	
}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }                         
}

