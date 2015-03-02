#include "MultiGrid1D.h"
#include "inclusion.h"


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
	int tempNumGrids = finestGridSize - 1;
	numGrids = log2(tempNumGrids); //finestGrid size = 65 = (2^6)+1 ---> numGrids = 6 (non creo griglia senza punti interni)

	grids1D = (Grid1D**)malloc(numGrids*sizeof(Grid1D*));
	grids1D[0] = new Grid1D(finestGridSize, range);
	for(int i = 1; i < numGrids; i++)
	{
		int coarseGridSize = ((grids1D[i-1]->sizeX-1)/2)+1; //coarseGridSize = ((2^/finerGridSize)-1)/2)+1
		grids1D[i] = new Grid1D(coarseGridSize, range);
	}
}

//fine->coarse X
void MultiGrid1D::Restrict(float* fine, int fsizeX, float* coarse, int csizeX)
{
	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni

	float C = 0;
	float E = 0;
	float O = 0;

	for(int cposX = 0; cposX < csizeX; cposX++)
	{
		/* COPIO VALORI SUI BORDI */
		if(cposX == 0 || cposX == csizeX - 1)
		{
			coarse[cposX] = fine[2*cposX];
			continue;
		}
		/* COPIO VALORI SUI BORDI */

		C = fine[2*cposX];
		E = fine[2*cposX+1];
		O = fine[2*cposX-1];

		coarse[cposX] = (1/4.0f)*(O + 2*C + E);
	}
}

void MultiGrid1D::Interpolate(float* fine, int fsizeX, float* coarse, int csizeX)
{
	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni

	for(int fposX = 0; fposX < fsizeX; fposX++)
	{
		int cposX = fposX/2; //divisione fra interi! 3/2=1

		if(fposX == 0 || fposX == fsizeX - 1)
			continue; //non modifico punti al contorno

		if(fposX%2 == 0) //punti pari = punti rossi
			fine[fposX] = coarse[cposX];
		else //punto dispari = punti neri
			fine[fposX] = (1/2.0f)*(coarse[cposX] + coarse[cposX+1]);
	}

}

void MultiGrid1D::Relax(Grid1D* curGrid, int ncycles)
{
	float* h_v = curGrid->h_v;
	float* h_f = curGrid->h_f;
	float h_x = curGrid->h_x;
	float h_x2 = h_x*h_x;
	int sizeX = curGrid->sizeX;

	float x_a = curGrid->x_a;

	for(int k = 0; k < ncycles; k++)
	{
		//i valori al limite non devono essere modificati
		for(int posX = 0; posX < sizeX; posX++)
		{
			if(posX%2 == 0) //punti pari = punti rossi
			{
				if(posX == 0 || posX == sizeX - 1)
					continue; //non modifico punti al contorno

				float xj = x_a + posX*h_x;

				h_v[posX] = (h_v[posX+1]*(exp(xj)+1) - h_f[posX]*h_x*(exp(xj) + 1))/(exp(xj)+1+h_x);
			}
		}

		//i valori al limite non devono essere modificati
		for(int posX = 0; posX < sizeX; posX++)
		{
			if(posX%2 != 0) //punti dispari = punti neri
			{
				if(posX == 0 || posX == sizeX - 1)
					continue; //non modifico punti al contorno

				float xj = x_a + posX*h_x;
				h_v[posX] = (h_v[posX+1]*(exp(xj)+1) - h_f[posX]*h_x*(exp(xj) + 1))/(exp(xj)+1+h_x);
			}
		}
	}
}

void MultiGrid1D::setToValue(float* grid, int sizeX, float value, bool modifyBoundaries)
{
	for(int posX = 0; posX < sizeX; posX++)
	{
			if(!modifyBoundaries)
				if(posX == 0 || posX == sizeX - 1)
					continue;

			grid[posX] = value;
	}
}

void MultiGrid1D::FullMultiGridVCycle(int gridID, int v0, int v1, int v2)
{
	Grid1D* fine = grids1D[gridID];
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		Grid1D* coarse = grids1D[gridID + 1];
		Restrict(fine->h_f, fine->sizeX, coarse->h_f,coarse->sizeX);
		FullMultiGridVCycle(gridID + 1,v0, v1, v2);
		Interpolate(fine->h_v, fine->sizeX, coarse->h_v, coarse->sizeX);
	}
	else
	{   //i contorni non vanno modificati in quanto nel FMG lavoriamo sulla sol approx e non sull'errore
		setToValue(fine->h_v, fine->sizeX, 0.0f, false);
	}
	for(int i = 0; i < v0; i++)
		VCycle(gridID, v1, v2);
}

void MultiGrid1D::VCycle(int gridID, int v1, int v2)
{
	Grid1D* fine = grids1D[gridID];
	Relax(fine, v1);
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		float* residual = CalculateResidual(fine);
		int residualSize = fine->sizeX;
		Grid1D* coarse = grids1D[gridID+1];
		Restrict(residual, residualSize, coarse->h_f, coarse->sizeX);

		//azzero tutto il vettore di errori (sui bordi sicuramente e=0)
		setToValue(coarse->h_v, coarse->sizeX, 0.0f, true);

		VCycle(gridID+1, v1, v2);

		int fsizeX = fine->sizeX;
		float* fine_error = (float*)malloc(fsizeX*sizeof(float)); //errore sulla griglia fine

		//coarse->v contiene errore su griglia rada
		Interpolate(fine_error, fsizeX, coarse->h_v, coarse->sizeX); //porti l'errore da lvl kh a lvl (k-1)h, nullo lungo bordi

		ApplyCorrection(fine->h_v, fsizeX, fine_error, fsizeX);
	}
	Relax(fine, v2);
}

void MultiGrid1D::ApplyCorrection(float* fine, int fineSizeX, float* fine_error, int errorSizeX)
{
	assert(fineSizeX == errorSizeX);
	//non bisogna modificare le soluzioni al contorno perchè lì la soluzione è esatta
	for(int posX = 0; posX < fineSizeX; posX++)
	{
		if(posX == 0 || posX == fineSizeX - 1)
			continue;

		fine[posX] = fine[posX] + fine_error[posX]; //u = v + e (u potrebbe essere usato come errore nella griglia un pò più fine)
	}
}

float* MultiGrid1D::CalculateResidual(Grid1D* fine)
{
	int sizeX = fine->sizeX;
	float* h_v = fine->h_v;
	float* h_f = fine->h_f;
	float h_x = fine->h_x;
	float x_a = fine->x_a;


	float* residual = (float*)malloc(sizeX*sizeof(float));

	for(int posX = 0; posX < sizeX; posX++)
	{
		//residuo sul contorno è nullo, in quanto lì la soluzione è esatta
		if(posX == 0 || posX == sizeX - 1)
		{
			residual[posX] = 0;
			continue;
		}
		float xj = x_a + posX*h_x;
		residual[posX]= h_f[posX] - (h_v[posX+1] - h_v[posX])/h_x - h_v[posX]/(exp(xj)+1);
	}

	return residual;
}

void MultiGrid1D::PrintDiff()
{
	int diff_fd = open("log/diff.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids1D[0]->PrintDiffApproxReal(diff_fd);
}

void MultiGrid1D::PrintGrid(int gridID)
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids1D[gridID]->PrintGrid_v(logfd);
}

void MultiGrid1D::PrintAllGrids_v()
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		grids1D[i]->PrintGrid_v(logfd);

	}
}

void MultiGrid1D::PrintAllGrids_f()
{
	int logfd = open("log/log_f.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		grids1D[i]->PrintGrid_f(logfd);
	}
}
