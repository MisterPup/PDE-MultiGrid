#include "Grid1D.h"
#include "MultiGrid1D.h"
#include "inclusion.h"

__global__ void CUDARestrict(float* fine, int fsizeX, float* coarse, int csizeX);
__global__ void CUDAInterpolate(float* fine, int fsizeX, float* coarse, int csizeX);
__global__ void CUDARelax(float* d_v, float* d_f, float h_x, int sizeX, float x_a);
__global__ void CUDASet(float* d_v, int sizeX, float value, bool modifyBoundaries);
__global__ void CUDACalculateResidual(float* d_v, float* d_f, float* d_residual, float h_x, int sizeX, float x_a);
__global__ void CUDAApplyCorrection(float* fine, float* error, int sizeX);

MultiGrid1D::MultiGrid1D(int finestGridSize, float range[])
{
	InitGrids(finestGridSize, range);
}

MultiGrid1D::~MultiGrid1D()
{
	for(int i = 0; i < numGrids; i++)
	{
		grids1D[i]->~Grid1D();
	}
	free(grids1D);
}

void MultiGrid1D::InitGrids(int finestGridSize, float range[])
{
	double tempNumGrids = finestGridSize - 1;
	numGrids = log2(tempNumGrids); //finestGrid size = 65 = (2^6)+1 ---> numGrids = 6 (non creo griglie senza punti interni)

	grids1D = (Grid1D**)malloc(numGrids*sizeof(Grid1D*));
	grids1D[0] = new Grid1D(finestGridSize, range);
	for(int i = 1; i < numGrids; i++)
	{
		int coarseGridSize = ((grids1D[i-1]->sizeX-1)/2)+1; //coarseGridSize = ((2^/finerGridSize)-1)/2)+1
		grids1D[i] = new Grid1D(coarseGridSize, range);
	}
}

void MultiGrid1D::Restrict(float* d_fine, int fsizeX, float* d_coarse, int csizeX)
{
	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni

	dim3 threadsPerBlock(NUMTHREAD_X); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)csizeX/threadsPerBlock.x);
	dim3 numBlocks(numBlocks_x);

	CUDARestrict<<<numBlocks, threadsPerBlock>>>(d_fine, fsizeX, d_coarse, csizeX);	
}
//OK
void MultiGrid1D::Interpolate(float* d_fine, int fsizeX, float* d_coarse, int csizeX)
{
	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni

	dim3 threadsPerBlock(NUMTHREAD_X); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)fsizeX/threadsPerBlock.x);
	dim3 numBlocks(numBlocks_x);
	
	CUDAInterpolate<<<numBlocks, threadsPerBlock>>>(d_fine, fsizeX, d_coarse, csizeX);
}

void MultiGrid1D::Relax(Grid1D* curGrid, int ncycles)
{
	float* d_v = curGrid->d_v;
	float* d_f = curGrid->d_f;
	float h_x = curGrid->h_x;
	int sizeX = curGrid->sizeX;
	float x_a = curGrid->x_a;

	dim3 threadsPerBlock(NUMTHREAD_X); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)sizeX/threadsPerBlock.x);
	dim3 numBlocks(numBlocks_x);

	for(int k = 0; k < ncycles; k++)                
		CUDARelax<<<numBlocks, threadsPerBlock>>>(d_v, d_f, h_x, sizeX, x_a);
}

void MultiGrid1D::Set(float* d_v, int sizeX,  float value, bool modifyBoundaries)
{
	dim3 threadsPerBlock(NUMTHREAD_X); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)sizeX/threadsPerBlock.x);
	dim3 numBlocks(numBlocks_x);

	CUDASet<<<numBlocks, threadsPerBlock>>>(d_v, sizeX, value, modifyBoundaries);
}

float* MultiGrid1D::CalculateResidual(Grid1D* curGrid)
{
	int sizeX = curGrid->sizeX;
	float* d_v = curGrid->d_v;
	float* d_f = curGrid->d_f;
	float h_x = curGrid->h_x;
	float x_a = curGrid->x_a;

	float* d_residual; //residuo (DEVICE)
	cudaMalloc(&d_residual, sizeX*sizeof(float));

	dim3 threadsPerBlock(NUMTHREAD_X); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)sizeX/threadsPerBlock.x);
	dim3 numBlocks(numBlocks_x);
							     
	CUDACalculateResidual<<<numBlocks, threadsPerBlock>>>(d_v, d_f, d_residual, h_x, sizeX, x_a);
	
	return d_residual;
}

void MultiGrid1D::ApplyCorrection(float* fine, int fineSizeX, float* error, int errorSizeX)
{
	assert(fineSizeX == errorSizeX); //f_pitch == e_pitch

	dim3 threadsPerBlock(NUMTHREAD_X); //il massimo numero di thread attivi per blocco è 768, ma ogni blocco può avere al max 512 thread 	
	int numBlocks_x = ceil((float)fineSizeX/threadsPerBlock.x);
	dim3 numBlocks(numBlocks_x);

	CUDAApplyCorrection<<<numBlocks, threadsPerBlock>>>(fine, error, fineSizeX);
}

void MultiGrid1D::VCycle(int gridID, int v1, int v2)
{
	Grid1D* fine = grids1D[gridID];
	Relax(fine, v1);
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		float* residual = CalculateResidual(fine);
		int residualSizeX = fine->sizeX;
		Grid1D* coarse = grids1D[gridID+1];
		Restrict(residual, residualSizeX, coarse->d_f, coarse->sizeX);	

		//azzero tutto il vettore di errori (sui bordi sicuramente e=0)
		Set(coarse->d_v, coarse->sizeX, 0.0f, true);

		VCycle(gridID+1, v1, v2);

		int fsizeX = fine->sizeX;
		float* fine_error;  //errore sulla griglia fine
		cudaMalloc(&fine_error, fsizeX*sizeof(float));
		//coarse->v contiene errore su griglia rada
		
		Interpolate(fine_error, fsizeX, coarse->d_v, coarse->sizeX); //porti l'errore da lvl kh a lvl (k-1)h, nullo lungo bordi
		
		ApplyCorrection(fine->d_v, fsizeX, fine_error, fsizeX);
	}

	Relax(fine, v2);

}

void MultiGrid1D::FullMultiGridVCycle(int gridID, int v0, int v1, int v2)
{
	Grid1D* fine = grids1D[gridID];
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		Grid1D* coarse = grids1D[gridID + 1];
		Restrict(fine->d_f, fine->sizeX, coarse->d_f, coarse->sizeX);
		FullMultiGridVCycle(gridID + 1, v0, v1, v2);
		Interpolate(fine->d_v, fine->sizeX, coarse->d_v, coarse->sizeX);
	}
	else //i contorni non vanno modificati in quanto nel FMG lavoriamo sulla sol approx e non sull'errore	   
		Set(fine->d_v, fine->sizeX, 0.0f, false);
	
	for(int i = 0; i < v0; i++)
		VCycle(gridID, v1, v2);
}

void MultiGrid1D::PrintDiff()
{
	int diff_fd = open("log/diff.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids1D[0]->PrintDiffApproxReal(diff_fd);
	close(diff_fd);
}

void MultiGrid1D::PrintGrid(int gridID)
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids1D[gridID]->PrintGrid_v(logfd);
	close(logfd);
}

void MultiGrid1D::PrintAllGrids_v()
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		grids1D[i]->PrintGrid_v(logfd);
	}
	close(logfd);
}

void MultiGrid1D::PrintAllGrids_f()
{
	int logfd = open("log/log_f.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		grids1D[i]->PrintGrid_f(logfd);

	}
	close(logfd);
}

/*********************************************************CUDA*********************************************************/

//f_pitch e c_pitch NON vanno in byte
//OK
__global__ void CUDARestrict(float* fine, int fsizeX, float* coarse, int csizeX)
{
	int cposX = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	
	if(cposX >= csizeX)
		return;

	if(cposX == 0 || cposX == csizeX - 1)
	{
		coarse[cposX] = fine[2*cposX];
		return;
	}

	//else non stiamo sui bordi
	float C = fine[2*cposX];
	float E = fine[2*cposX+1];
	float O = fine[2*cposX-1];

	coarse[cposX] = (1/4.0f)*(O + 2*C + E);

}

//f_pitch e c_pitch NON vanno in byte
//OK
__global__ void CUDAInterpolate(float* fine, int fsizeX, float* coarse, int csizeX)
{
	int fposX = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	
	if(fposX >= fsizeX)
		return;

	int cposX = fposX/2; //divisione fra interi! 3/2=1

	if(fposX == 0 || fposX == fsizeX - 1)
		return; //non modifico punti al contorno

	if(fposX%2 == 0) //punti pari = punti rossi
		fine[fposX] = coarse[cposX];
	else //punto dispari = punti neri
		fine[fposX] = (1/2.0f)*(coarse[cposX] + coarse[cposX+1]);

}

__global__ void CUDARelax(float* d_v, float* d_f, float h_x, int sizeX, float x_a)
{
	int posX = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	
	if(posX >= sizeX)
		return;

	if(posX == 0 || posX == sizeX - 1)
		return; //non modifico punti al contorno
		
	if(posX%2 == 0)  //punti pari = punti rossi
	{
		float xj = x_a + posX*h_x;
		d_v[posX] = (d_v[posX+1]*(exp(xj)+1) - d_f[posX]*h_x*(exp(xj) + 1))/(exp(xj)+1+h_x);
	}

	__syncthreads();//prima di aggiornare i punti neri(che dipendono esclusivamente da quelli rossi), bisogna aspettare che ogni punto rosso sia stato aggiornato 

	if(posX%2 != 0) //punti dispari = punti neri
	{
		float xj = x_a + posX*h_x;
		d_v[posX] = (d_v[posX+1]*(exp(xj)+1) - d_f[posX]*h_x*(exp(xj) + 1))/(exp(xj)+1+h_x);
	}
}

__global__ void CUDASet(float* d_v, int sizeX, float value, bool modifyBoundaries)
{
	int posX = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	
	if(posX >= sizeX)
		return;

	if(!modifyBoundaries)
		if(posX == 0 || posX == sizeX - 1)
			return;

	//else
	d_v[posX] = value;
}

__global__ void CUDACalculateResidual(float* d_v, float* d_f, float* d_residual, float h_x, int sizeX, float x_a)
{
	int posX = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	
	if(posX >= sizeX)
		return;

	//residuo sul contorno è nullo, in quanto lì la soluzione è esatta
	if(posX == 0 || posX == sizeX - 1)
	{
		d_residual[posX] = 0;
		return;
	}
	
	//else
	float xj = x_a + posX*h_x;
	d_residual[posX]= d_f[posX] - (d_v[posX+1] - d_v[posX])/h_x - d_v[posX]/(exp(xj)+1);
}
                                   
__global__ void CUDAApplyCorrection(float* fine, float* error, int sizeX)
{
	int posX = blockDim.x * blockIdx.x + threadIdx.x; //indice j
	
	if(posX >= sizeX)
		return;

	if(posX == 0 || posX == sizeX - 1)
		return;

	fine[posX] = fine[posX] + error[posX]; //u = v + e (u potrebbe essere usato come errore nella griglia un pò più fine)
}

