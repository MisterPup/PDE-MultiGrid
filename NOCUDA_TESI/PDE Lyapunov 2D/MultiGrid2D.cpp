#include "MultiGrid2D.h"
#include "inclusion.h"

//Corretto
MultiGrid2D::MultiGrid2D(int finestGridSizeXY[], float range[], float* _A, int A_size, int alfa)
{
		InitGrids(finestGridSizeXY, range);
		InitA(_A, A_size, alfa);
}
//Corretto
MultiGrid2D::~MultiGrid2D()
{
	for(int i = 0; i < numGrids; i++)
	{
		grids2D[i]->~Grid2D();
	}
	free(grids2D);
}
//Corretto
void MultiGrid2D::InitGrids(int finestGridSizeXY[], float range[])
{
	int finestGridSizeX = finestGridSizeXY[0];
	int finestGridSizeY = finestGridSizeXY[1];

	int minSize = finestGridSizeX;
	if(finestGridSizeY < minSize)
		minSize = finestGridSizeY;

	int tempNumGrids = minSize - 1;
	numGrids = log2(tempNumGrids); //finestGrid size = 65 = (2^6)+1 ---> numGrids = 6 (non creo griglie senza punti interni)

	grids2D = (Grid2D**)malloc(numGrids*sizeof(Grid2D*));
	grids2D[0] = new Grid2D(finestGridSizeXY, range);
	for(int i = 1; i < numGrids; i++)
	{
		int coarseGridSizeX = ((grids2D[i-1]->sizeX-1)/2)+1; //coarseGridSize = ((2^/finerGridSize)-1)/2)+1
		int coarseGridSizeY = ((grids2D[i-1]->sizeY-1)/2)+1;

		int coarseGridSizeXY[2] = {coarseGridSizeX, coarseGridSizeY};
		grids2D[i] = new Grid2D(coarseGridSizeXY, range);
	}
}

//Corretto
void MultiGrid2D::InitA(float* _A, int A_size, int alfa)
{
	this->alfa = alfa;
	this->sizeA = A_size;

	matrixA = (float*)malloc(A_size*sizeof(float));

	for(int posY = 0; posY < sizeA; posY++)
	{
		for(int posX = 0; posX < sizeA; posX++)
		{
			int idx = posX + posY * sizeA;
			matrixA[idx] = _A[idx];
		}
	}
}

//fine->coarse
void MultiGrid2D::Restrict(float* fine, int fsizeXY[], float* coarse, int csizeXY[])
{
	int fsizeX = fsizeXY[0];
	int fsizeY = fsizeXY[1];

	int csizeX = csizeXY[0];
	int csizeY = csizeXY[1];

	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni
	assert(csizeY == ((fsizeY-1)/2+1));

	float C = 0;
	float N = 0;
	float S = 0;
	float E = 0;
	float O = 0;
	float NE = 0;
	float NO = 0;
	float SE = 0;
	float SO = 0;

	//valorizza punti interni ed esterni alla grid
	for(int cposY = 0; cposY < csizeY; cposY++)
	{
		for(int cposX = 0; cposX < csizeX; cposX++)
		{
			int fposX = 2*cposX;
			int fposY = 2*cposY;
			int cidx = 0;
			int fidx = 0;

			/* COPIO VALORI SUI BORDI */
			if(cposX == 0 || cposX == csizeX - 1 || cposY == 0 || cposY == csizeY - 1)
			{
				cidx = cposX + cposY * csizeX;
				fidx = fposX + fposY * fsizeX;
				coarse[cidx] = fine[fidx];
				continue;
			}

			fidx = fposX + fposY * fsizeX;
			C = fine[fidx];
			fidx = fposX + (fposY-1) * fsizeX;
			N = fine[fidx];
			fidx = fposX + (fposY+1) * fsizeX;
			S = fine[fidx];
			fidx = (fposX+1) + fposY * fsizeX;
			E = fine[fidx];
			fidx = (fposX-1) + fposY * fsizeX;
			O = fine[fidx];
			fidx = (fposX+1) + (fposY-1) * fsizeX;
			NE = fine[fidx];
			fidx = (fposX-1) + (fposY-1) * fsizeX;
			NO = fine[fidx];
			fidx = (fposX+1) + (fposY+1) * fsizeX;
			SE = fine[fidx];
			fidx = (fposX-1) + (fposY+1) * fsizeX;
			SO = fine[fidx];

			cidx = cposX + cposY * csizeX;
			coarse[cidx] = (1/16.0f)*(NO+NE+SO+SE + 2*(O+E+N+S) + 4*C);
		}
	}
}

void MultiGrid2D::Interpolate(float* fine, int fsizeXY[], float* coarse, int csizeXY[])
{
	int fsizeX = fsizeXY[0];
	int fsizeY = fsizeXY[1];

	int csizeX = csizeXY[0];
	int csizeY = csizeXY[1];

	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni
	assert(csizeY == ((fsizeY-1)/2+1));

	//valorizza punti interni alla grid, i valori sul contorno rimangono invariati
	for(int fposY = 0; fposY < fsizeY; fposY++)
	{
		for(int fposX = 0; fposX < fsizeX; fposX++)
		{
			int cposX = fposX/2; //divisione fra interi! 3/2=1
			int cposY = fposY/2;

			int fidx = fposX + fposY * fsizeX;
			int cidx = cposX + cposY * csizeX;

			if(fposX == 0 || fposX == fsizeX - 1 || fposY == 0 || fposY == fsizeY - 1)
				continue; //non modifico punti al contorno

			if(fposY%2 == 0 && fposX%2 == 0)
			{
				fine[fidx] = coarse[cidx];
				continue;
			}
			if(fposY%2 != 0 && fposX%2 == 0)
			{
				cidx = cposX + cposY * csizeX;
				float N = coarse[cidx];
				cidx = cposX + (cposY+1) * csizeX;
				float S = coarse[cidx];

				fine[fidx] = (1/2.0f)*(N + S);
				continue;
			}

			if(fposY%2 == 0 && fposX%2 != 0)
			{
				cidx = cposX + cposY * csizeX;
				float O = coarse[cidx];
				cidx = (cposX+1) + cposY * csizeX;
				float E = coarse[cidx];

				fine[fidx] = (1/2.0f)*(O + E);
				continue;
			}

			if(fposY%2 != 0 && fposX%2 != 0)
			{
				cidx = cposX + cposY * csizeX;
				float NO = coarse[cidx];
				cidx = (cposX+1) + cposY * csizeX;
				float NE = coarse[cidx];
				cidx = cposX + (cposY+1) * csizeX;
				float SO = coarse[cidx];
				cidx = (cposX+1) + (cposY+1) * csizeX;
				float SE = coarse[cidx];

				fine[fidx] = (1/4.0f)*(NO + NE + SO + SE);
				continue;
			}
		}
	}
}


void MultiGrid2D::Relax(Grid2D* curGrid, int ncycles)
{
	float* h_v = curGrid->h_v;
	float* f = curGrid->h_f;

	float h_x = curGrid->h_x;
	float h_y = curGrid->h_y;

	int sizeX = curGrid->sizeX;
	int sizeY = curGrid->sizeY;

	float x_a = curGrid->x_a;
	float y_a = curGrid->y_a;

	float K1 = 0.0f;
	float K2 = 0.0f;

	for(int k = 0; k < ncycles; k++)
	{
		//i valori al limite non devono essere modificati
		for(int posY = 0; posY < sizeY; posY++)
		{
			for(int posX = 0; posX < sizeX; posX++)
			{
				if((posY+ posX)%2 == 0) //punti pari
				{
					int idx = posX + posY * sizeX;

					if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1)
						continue; //non modifico punti al contorno

					float xj = x_a + posX*h_x;
					float yi = y_a + posY*h_y;

					K1 = matrixA[0]*xj + matrixA[1]*yi;
					K2 = matrixA[2]*xj + matrixA[3]*yi;

					float den = K1*h_y+K2*h_x-alfa*h_x*h_y;

					int idxVarX = (posX+1) + posY * sizeX;
					int idxVarY = posX + (posY+1) * sizeX;

					h_v[idx] = (h_y*K1*h_v[idxVarX] + h_x*K2*h_v[idxVarY] - f[idx]*h_x*h_y)/(den);
				}
			}
		}
		//i valori al limite non devono essere modificati
		for(int posY = 0; posY < sizeY; posY++)
		{
			for(int posX = 0; posX < sizeX; posX++)
			{
				if((posY+ posX)%2 != 0) //punti dispari
				{
					int idx = posX + posY * sizeX;

					if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1)
						continue; //non modifico punti al contorno

					float xj = x_a + posX*h_x;
					float yi = y_a + posY*h_y;

					K1 = matrixA[0]*xj + matrixA[1]*yi;
					K2 = matrixA[2]*xj + matrixA[3]*yi;

					float den = K1*h_y+K2*h_x-alfa*h_x*h_y;

					int idxVarX = (posX+1) + posY * sizeX;
					int idxVarY = posX + (posY+1) * sizeX;

					h_v[idx] = (h_y*K1*h_v[idxVarX] + h_x*K2*h_v[idxVarY] - f[idx]*h_x*h_y)/(den);
				}
			}
		}
	}
}

void MultiGrid2D::setToValue(float* grid, int sizeXY[], float value, bool modifyBoundaries)
{
	int sizeX = sizeXY[0];
	int sizeY = sizeXY[1];

	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			if(!modifyBoundaries)
				if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1)
					continue;

			int idx = posX + posY * sizeX;
			grid[idx] = value;
		}
	}
}



void MultiGrid2D::FullMultiGridVCycle(int gridID, int v0, int v1, int v2)
{
	Grid2D* fine = grids2D[gridID];
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		Grid2D* coarse = grids2D[gridID + 1];
		Restrict(fine->h_f, fine->sizeXY, coarse->h_f,coarse->sizeXY);
		FullMultiGridVCycle(gridID + 1, v0, v1, v2);
		Interpolate(fine->h_v, fine->sizeXY, coarse->h_v,coarse->sizeXY);
	}
	else
	{   //i contorni non vanno modificati in quanto nel FMG lavoriamo sulla sol approx e non sull'errore
		setToValue(fine->h_v, fine->sizeXY, 0.0f, false);
	}
	for(int i = 0; i < v0; i++)
		VCycle(gridID, v1, v2);
}

void MultiGrid2D::VCycle(int gridID, int v1, int v2)
{
	Grid2D* fine = grids2D[gridID];
	Relax(fine, v1);
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		float* residual = CalculateResidual(fine);
		int* residualSize = fine->sizeXY;
		Grid2D* coarse = grids2D[gridID+1];
		Restrict(residual, residualSize, coarse->h_f, coarse->sizeXY);

		//azzero tutto il vettore di errori (sui bordi sicuramente e=0)
		setToValue(coarse->h_v, coarse->sizeXY, 0.0f, true);

		VCycle(gridID+1, v1, v2);

		int* fsizeXY = fine->sizeXY;
		float* fine_error = (float*)malloc(fsizeXY[0]*fsizeXY[1]*sizeof(float)); //errore sulla griglia fine
		//coarse->v contiene errore su griglia rada
		Interpolate(fine_error, fsizeXY, coarse->h_v, coarse->sizeXY); //porti l'errore da lvl kh a lvl (k-1)h, nullo lungo bordi

		ApplyCorrection(fine->h_v, fsizeXY, fine_error, fsizeXY);
	}

	Relax(fine, v2);

}


void MultiGrid2D::ApplyCorrection(float* fine,  int fsizeXY[], float* error,  int esizeXY[])
{
	int fsizeX = fsizeXY[0];
	int fsizeY = fsizeXY[1];

	int esizeX = esizeXY[0];
	int esizeY = esizeXY[1];

	assert(fsizeX == esizeX);
	assert(fsizeY == esizeY);

	//non bisogna modificare le soluzioni al contorno perchè lì la soluzione è esatta
	for(int fposY = 0; fposY < fsizeY; fposY++)
	{
		for(int fposX = 0; fposX < fsizeX; fposX++)
		{
			if(fposX == 0 || fposX == fsizeX - 1 || fposY == 0 || fposY == fsizeY - 1)
				continue;

			int idx = fposX + fposY * fsizeX;
			fine[idx] = fine[idx] + error[idx]; //u = v + e (u potrebbe essere usato come errore nella griglia un pò più fine)
		}
	}
}
float* MultiGrid2D::CalculateResidual(Grid2D* fine)
{
	float* h_v = fine->h_v;
	float* h_f = fine->h_f;

	float h_x = fine->h_x;
	float h_y = fine->h_y;

	int sizeX = fine->sizeX;
	int sizeY = fine->sizeY;

	float y_a = fine->y_a;
	float x_a = fine->x_a;

	float* residual = (float*)malloc(sizeX*sizeY*sizeof(float));

	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			int idx = posX + posY * sizeX;

			if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1)
			{
				residual[idx] = 0.0f; //residuo nullo sui bordi
			}
			else
			{
				float xj = x_a + posX*h_x;
				float yi = y_a + posY*h_y;
				float K1 = matrixA[0]*xj + matrixA[1]*yi;
				float K2 = matrixA[2]*xj + matrixA[3]*yi;

				int idxVarX = (posX+1) + posY * sizeX;
				int idxVarY = posX + (posY+1) * sizeX;

				residual[idx] = h_f[idx] - (h_y*K1*h_v[idxVarX] + h_x*K2*h_v[idxVarY] - h_v[idx]*(h_y*K1+h_x*K2-alfa*h_x*h_y))/(h_x*h_y);
			}
		}
	}
	return residual;
}

void MultiGrid2D::PrintDiff()
{
	int diff_fd = open("log/diff.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids2D[0]->PrintDiffApproxReal(diff_fd);
}

void MultiGrid2D::PrintGrid(int gridID)
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids2D[gridID]->PrintGrid_v(logfd);
}

void MultiGrid2D::PrintAllGrids_v()
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		printf("Grid: %d\n", i);
		grids2D[i]->PrintGrid_v(logfd);

	}
}

void MultiGrid2D::PrintAllGrids_f()
{
	int logfd = open("log/log_f.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		printf("Grid: %d\n", i);
		grids2D[i]->PrintGrid_f(logfd);

	}
}

void MultiGrid2D::PrintResidual(int gridID)
{
	int logfd = open("log/residual.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids2D[gridID]->PrintResidual(logfd, matrixA, alfa);
}

