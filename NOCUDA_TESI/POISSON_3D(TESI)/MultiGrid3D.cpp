#include "MultiGrid3D.h"
#include "inclusion.h"

//Corretto
MultiGrid3D::MultiGrid3D(int finestGridSizeXYZ[], float range[])
{
		InitGrids(finestGridSizeXYZ, range);
}
//Corretto
MultiGrid3D::~MultiGrid3D()
{
	for(int i = 0; i < numGrids; i++)
	{
		grids3D[i]->~Grid3D();
	}
	free(grids3D);
}
//Corretto
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


	int tempNumGrids = minSize - 1;
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

//fine->coarse Corretto
void MultiGrid3D::Restrict(float* fine, int fsizeXYZ[], float* coarse, int csizeXYZ[])
{
	int fsizeX = fsizeXYZ[0];
	int fsizeY = fsizeXYZ[1];
	int fsizeZ = fsizeXYZ[2];

	int csizeX = csizeXYZ[0];
	int csizeY = csizeXYZ[1];
	int csizeZ = csizeXYZ[2];

	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni
	assert(csizeY == ((fsizeY-1)/2+1));
	assert(csizeZ == ((fsizeZ-1)/2+1));

	/* CENTRO */
	float C_C = 0;
	float N_C = 0;
	float S_C = 0;
	float E_C = 0;
	float O_C = 0;
	float NE_C = 0;
	float NO_C = 0;
	float SE_C = 0;
	float SO_C = 0;
	/* CENTRO */

	/* NORD */
	float C_N = 0;
	float N_N = 0;
	float S_N = 0;
	float E_N = 0;
	float O_N = 0;
	float NE_N = 0;
	float NO_N = 0;
	float SE_N = 0;
	float SO_N = 0;
	/* NORD */

	/* SUD */
	float C_S = 0;
	float N_S = 0;
	float S_S = 0;
	float E_S = 0;
	float O_S = 0;
	float NE_S = 0;
	float NO_S = 0;
	float SE_S = 0;
	float SO_S = 0;
	/* SUD */

	//valorizza punti interni ed esterni alla grid
	for(int cposY = 0; cposY < csizeY; cposY++)
	{
		for(int cposX = 0; cposX < csizeX; cposX++)
		{
			for(int cposZ = 0; cposZ < csizeZ; cposZ++)
			{
				int fposX = 2*cposX;
				int fposY = 2*cposY;
				int fposZ = 2*cposZ;
				int idx = 0;

				/* COPIO VALORI SUI BORDI */
				if(cposX == 0 || cposX == csizeX - 1 || cposY == 0 || cposY == csizeY - 1 || cposZ == 0 || cposZ == csizeZ - 1)
				{
					int cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
					int fidx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY;
					coarse[cidx] = fine[fidx];
					continue;
				}
				/* COPIO VALORI SUI BORDI */

				idx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY;
				C_C = fine[idx];
				idx = fposX + fposY * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				N_C = fine[idx];
				idx = fposX + fposY * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				S_C = fine[idx];
				idx = (fposX+1) + fposY * fsizeX + fposZ * fsizeX * fsizeY;
				E_C = fine[idx];
				idx = (fposX-1) + fposY * fsizeX + fposZ * fsizeX * fsizeY;
				O_C = fine[idx];
				idx = (fposX+1) + fposY * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				NE_C = fine[idx];
				idx = (fposX-1) + fposY * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				NO_C = fine[idx];
				idx = (fposX+1) + fposY * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				SE_C = fine[idx];
				idx = (fposX-1) + fposY * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				SO_C = fine[idx];

				idx = fposX + (fposY-1) * fsizeX + fposZ * fsizeX * fsizeY;
				C_N = fine[idx];
				idx = fposX + (fposY-1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				N_N = fine[idx];
				idx = fposX + (fposY-1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				S_N = fine[idx];
				idx = (fposX+1) + (fposY-1) * fsizeX + (fposZ) * fsizeX * fsizeY;
				E_N = fine[idx];
				idx = (fposX-1) + (fposY-1) * fsizeX + (fposZ) * fsizeX * fsizeY;
				O_N = fine[idx];
				idx = (fposX+1) + (fposY-1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				NE_N = fine[idx];
				idx = (fposX-1) + (fposY-1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				NO_N = fine[idx];
				idx = (fposX+1) + (fposY-1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				SE_N = fine[idx];
				idx = (fposX-1) + (fposY-1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				SO_N = fine[idx];

				idx = fposX + (fposY+1) * fsizeX + fposZ * fsizeX * fsizeY;
				C_S = fine[idx];
				idx = fposX + (fposY+1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				N_S = fine[idx];
				idx = fposX + (fposY+1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				S_S = fine[idx];
				idx = (fposX+1) + (fposY+1) * fsizeX + (fposZ) * fsizeX * fsizeY;
				E_S = fine[idx];
				idx = (fposX-1) + (fposY+1) * fsizeX + (fposZ) * fsizeX * fsizeY;
				O_S = fine[idx];
				idx = (fposX+1) + (fposY+1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				NE_S = fine[idx];
				idx = (fposX-1) + (fposY+1) * fsizeX + (fposZ+1) * fsizeX * fsizeY;
				NO_S = fine[idx];
				idx = (fposX+1) + (fposY+1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				SE_S = fine[idx];
				idx = (fposX-1) + (fposY+1) * fsizeX + (fposZ-1) * fsizeX * fsizeY;
				SO_S = fine[idx];

				idx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
				coarse[idx] = (1/8.0f)*(C_C) + (1/16.0f)*((N_C+E_C+S_C+O_C) + (C_N+C_S)) + (1/32.0f)*((NE_C+SE_C+SO_C+NO_C) + (N_N+E_N+S_N+O_N) + (N_S+E_S+S_S+O_S)) + (1/64.0f)*((NE_N+SE_N+SO_N+NO_N) + (NE_S+SE_S+SO_S+NO_S));
			}
		}
	}
}
//MESSA IO ORA!!!!! NON MODIFICA I BORDI!!!! RISULTATO UGUALE A QUELLA SOTTO DI INTERPOLAZIONE
void MultiGrid3D::Interpolate(float* fine, int fsizeXYZ[], float* coarse, int csizeXYZ[])
{
	int fsizeX = fsizeXYZ[0];
	int fsizeY = fsizeXYZ[1];
	int fsizeZ = fsizeXYZ[2];

	int csizeX = csizeXYZ[0];
	int csizeY = csizeXYZ[1];
	int csizeZ = csizeXYZ[2];

	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni
	assert(csizeY == ((fsizeY-1)/2+1));
	assert(csizeZ == ((fsizeZ-1)/2+1));


	//valorizza punti interni alla grid, i valori sul contorno rimangono invariati
	for(int fposY = 1; fposY < fsizeY - 1; fposY++)
	{
		for(int fposX = 1; fposX < fsizeX - 1; fposX++)
		{
			for(int fposZ = 1; fposZ < fsizeZ - 1; fposZ++)
			{
				int cposX = fposX/2; //divisione fra interi! 3/2=1
				int cposY = fposY/2;
				int cposZ = fposZ/2;

				int fidx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY;
				int cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;


				if(fposY%2 == 0 && fposX%2 == 0 && fposZ%2 == 0) //PPP
				{
					fine[fidx] = coarse[cidx];
					continue;
				}

				if(fposY%2 == 0 && fposX%2 != 0 && fposZ%2 == 0) //PDP
				{
					cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
					float O = coarse[cidx];
					cidx = (cposX+1) + cposY * csizeX + cposZ * csizeX * csizeY;
					float E = coarse[cidx];

					fine[fidx] = (1/2.0f)*(O + E);
					continue;
				}

				if(fposY%2 != 0 && fposX%2 == 0 && fposZ%2 == 0) //DPP
				{
					cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
					float N = coarse[cidx];
					cidx = cposX + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
					float S = coarse[cidx];

					fine[fidx] = (1/2.0f)*(N + S);
					continue;
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
					continue;
				}

				/* *******************************************************/

				if(fposY%2 == 0 && fposX%2 == 0 && fposZ%2 != 0) //PPD
				{
					cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
					float S = coarse[cidx];
					cidx = cposX + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
					float N = coarse[cidx];

					fine[fidx] = (1/2.0f)*(S + N);
					continue;
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
					continue;
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
					continue;
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
					continue;
				}
			}
		}
	}
}
/*
void MultiGrid3D::Interpolate(float* fine, int fsizeXYZ[], float* coarse, int csizeXYZ[])
{
	int fsizeX = fsizeXYZ[0];
	int fsizeY = fsizeXYZ[1];
	int fsizeZ = fsizeXYZ[2];

	int csizeX = csizeXYZ[0];
	int csizeY = csizeXYZ[1];
	int csizeZ = csizeXYZ[2];

	assert(csizeX == ((fsizeX-1)/2+1)); //controllo su dimensioni
	assert(csizeY == ((fsizeY-1)/2+1));
	assert(csizeZ == ((fsizeZ-1)/2+1));

	//NON DOVREI NON MODIFICARE I BORDI???????

	//valorizza punti interni alla grid, i valori sul contorno rimangono invariati
	for(int fposY = 0; fposY < fsizeY; fposY++)
	{
		for(int fposX = 0; fposX < fsizeX; fposX++)
		{
			for(int fposZ = 0; fposZ < fsizeZ; fposZ++)
			{
				int cposX = fposX/2; //divisione fra interi! 3/2=1
				int cposY = fposY/2;
				int cposZ = fposZ/2;

				int fidx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY;
				int cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;


				if(fposY%2 == 0 && fposX%2 == 0 && fposZ%2 == 0) //PPP
				{
					fine[fidx] = coarse[cidx];
					continue;
				}

				if(fposY%2 == 0 && fposX%2 != 0 && fposZ%2 == 0) //PDP
				{
					cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
					float O = coarse[cidx];
					cidx = (cposX+1) + cposY * csizeX + cposZ * csizeX * csizeY;
					float E = coarse[cidx];

					fine[fidx] = (1/2.0f)*(O + E);
					continue;
				}

				if(fposY%2 != 0 && fposX%2 == 0 && fposZ%2 == 0) //DPP
				{
					cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
					float N = coarse[cidx];
					cidx = cposX + (cposY+1) * csizeX + cposZ * csizeX * csizeY;
					float S = coarse[cidx];

					fine[fidx] = (1/2.0f)*(N + S);
					continue;
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
					continue;
				}
*/
				/* *******************************************************/
/*
				if(fposY%2 == 0 && fposX%2 == 0 && fposZ%2 != 0) //PPD
				{
					cidx = cposX + cposY * csizeX + cposZ * csizeX * csizeY;
					float S = coarse[cidx];
					cidx = cposX + cposY * csizeX + (cposZ+1) * csizeX * csizeY;
					float N = coarse[cidx];

					fine[fidx] = (1/2.0f)*(S + N);
					continue;
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
					continue;
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
					continue;
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
					continue;
				}
			}
		}
	}
}*/

void MultiGrid3D::Relax(Grid3D* curGrid, int ncycles)
{
	float* v = curGrid->h_v;
	float* f = curGrid->h_f;

	float h_x = curGrid->h_x;
	float h_y = curGrid->h_y;
	float h_z = curGrid->h_z;

	float h_x2 = h_x*h_x;
	float h_y2 = h_y*h_y;
	float h_z2 = h_z*h_z;

	int sizeX = curGrid->sizeX;
	int sizeY = curGrid->sizeY;
	int sizeZ = curGrid->sizeZ;

	for(int k = 0; k < ncycles; k++)
	{
		//i valori al limite non devono essere modificati
		for(int posY = 1; posY < sizeY - 1; posY++)
		{
			for(int posX = 1; posX < sizeX - 1; posX++)
			{
				for(int posZ = 1; posZ < sizeZ - 1; posZ++)
				{
					if((posY+ posX + posZ)%2 == 0) //punti pari
					{
						int idx = 0;
						idx = (posX-1) + posY * sizeX + posZ * sizeX * sizeY;
						float O = v[idx];
						idx = (posX+1) + posY * sizeX + posZ * sizeX * sizeY;
						float E = v[idx];
						idx = posX + (posY-1) * sizeX + posZ * sizeX * sizeY;
						float N = v[idx];
						idx = posX + (posY+1) * sizeX + posZ * sizeX * sizeY;
						float S = v[idx];
						idx = posX + posY * sizeX + (posZ-1) * sizeX * sizeY;
						float D = v[idx];
						idx = posX + posY * sizeX + (posZ+1) * sizeX * sizeY;
						float U = v[idx];

						idx = posX + posY * sizeX + posZ * sizeX * sizeY;
						v[idx] = (O*(h_y2*h_z2)+E*(h_y2*h_z2) + N*(h_x2*h_z2)+S*(h_x2*h_z2) + D*(h_x2*h_y2)+U*(h_x2*h_y2) - f[idx]*h_x2*h_y2*h_z2)/(2*(h_y2*h_z2 + h_x2*h_z2 + h_x2*h_y2));
					}
				}
			}
		}

		for(int posY = 1; posY < sizeY - 1; posY++)
		{
			for(int posX = 1; posX < sizeX - 1; posX++)
			{
				for(int posZ = 1; posZ < sizeZ - 1; posZ++)
				{
					if((posY+ posX + posZ)%2 != 0) //punti dispari
					{
						int idx = 0;
						idx = (posX-1) + posY * sizeX + posZ * sizeX * sizeY;
						float O = v[idx];
						idx = (posX+1) + posY * sizeX + posZ * sizeX * sizeY;
						float E = v[idx];
						idx = posX + (posY-1) * sizeX + posZ * sizeX * sizeY;
						float N = v[idx];
						idx = posX + (posY+1) * sizeX + posZ * sizeX * sizeY;
						float S = v[idx];
						idx = posX + posY * sizeX + (posZ-1) * sizeX * sizeY;
						float D = v[idx];
						idx = posX + posY * sizeX + (posZ+1) * sizeX * sizeY;
						float U = v[idx];

						idx = posX + posY * sizeX + posZ * sizeX * sizeY;
						v[idx] = (O*(h_y2*h_z2)+E*(h_y2*h_z2) + N*(h_x2*h_z2)+S*(h_x2*h_z2) + D*(h_x2*h_y2)+U*(h_x2*h_y2) - f[idx]*h_x2*h_y2*h_z2)/(2*(h_y2*h_z2 + h_x2*h_z2 + h_x2*h_y2));
					}
				}
			}
		}
	}
}

void MultiGrid3D::FullMultiGridVCycle(int gridID, int v0, int v1, int v2)
{
	Grid3D* fine = grids3D[gridID];
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		Grid3D* coarse = grids3D[gridID + 1];
		Restrict(fine->h_f, fine->sizeXYZ, coarse->h_f,coarse->sizeXYZ);
		FullMultiGridVCycle(gridID + 1,v0, v1, v2);
		Interpolate(fine->h_v, fine->sizeXYZ, coarse->h_v,coarse->sizeXYZ);
	}
	else
	{   //i contorni non vanno modificati in quanto nel FMG lavoriamo sulla sol approx e non sull'errore
		setToValue(fine->h_v, fine->sizeXYZ, 0.0f, false);
	}
	for(int i = 0; i < v0; i++)
		VCycle(gridID, v1, v2);
}

void MultiGrid3D::setToValue(float* grid, int sizeXYZ[], float value, bool modifyBoundaries)
{
	int sizeX = sizeXYZ[0];
	int sizeY = sizeXYZ[1];
	int sizeZ = sizeXYZ[2];

	if(modifyBoundaries)
	{
		for(int posY = 0; posY < sizeY; posY++)
		{
			for(int posX = 0; posX < sizeX; posX++)
			{
				for(int posZ = 0; posZ < sizeZ; posZ++)
				{
					int idx = posX + posY * sizeX + posZ * sizeX * sizeY;
					grid[idx] = value;
				}
			}
		}
	}
	if(!modifyBoundaries)
	{
		for(int posY = 1; posY < sizeY - 1; posY++)
		{
			for(int posX = 1; posX < sizeX - 1; posX++)
			{
				for(int posZ = 1; posZ < sizeZ - 1; posZ++)
				{
					int idx = posX + posY * sizeX + posZ * sizeX * sizeY;
					grid[idx] = value;
				}
			}
		}
	}
}

void MultiGrid3D::VCycle(int gridID, int v1, int v2)
{
	Grid3D* fine = grids3D[gridID];
	Relax(fine, v1);
	if(gridID != numGrids - 1) //la griglia corrente non è la più rada
	{
		float* residual = CalculateResidual(fine);
		int* residualSize = fine->sizeXYZ;
		Grid3D* coarse = grids3D[gridID+1];
		Restrict(residual, residualSize, coarse->h_f, coarse->sizeXYZ);
		//azzero tutto il vettore di errori (sui bordi sicuramente e=0)
		setToValue(coarse->h_v, coarse->sizeXYZ, 0.0f, true); //forse sui bordi non è proprio così per via delle condizioni inziali strane
		VCycle(gridID+1, v1, v2);

		int* fsize = fine->sizeXYZ;
		float* fine_error = (float*)malloc(fsize[0]*fsize[1]*fsize[2]*sizeof(float)); //errore sulla griglia fine
		//coarse->v contiene errore su griglia rada
		Interpolate(fine_error, fsize, coarse->h_v, coarse->sizeXYZ); //porti l'errore da lvl kh a lvl (k-1)h, nullo lungo bordi

		ApplyCorrection(fine->h_v, fsize, fine_error, fsize);
	}

	Relax(fine, v2);

}

void MultiGrid3D::ApplyCorrection(float* fine,  int fsizeXYZ[], float* error,  int esizeXYZ[])
{
	int fsizeX = fsizeXYZ[0];
	int fsizeY = fsizeXYZ[1];
	int fsizeZ = fsizeXYZ[2];

	int esizeX = esizeXYZ[0];
	int esizeY = esizeXYZ[1];
	int esizeZ = esizeXYZ[2];

	//potresti anche lavorare sui bordi, tanto non vengono modificati perchè l'errore è sicuramente nullo sugli stessi
	assert(fsizeX == esizeX);
	assert(fsizeY == esizeY);
	assert(fsizeZ == esizeZ);

	//non bisogna toccare le soluzioni al contorno perchè lì la soluzione è esatta
	for(int fposY = 1; fposY < fsizeY - 1; fposY++)
	{
		for(int fposX = 1; fposX < fsizeX - 1; fposX++)
		{
			for(int fposZ = 1; fposZ < fsizeZ - 1; fposZ++)
			{
				int idx = fposX + fposY * fsizeX + fposZ * fsizeX * fsizeY;
				fine[idx] = fine[idx] + error[idx]; //u = v + e (u potrebbe essere usato come errore nella griglia un pò più fine)
			}
		}
	}
}

float* MultiGrid3D::CalculateResidual(Grid3D* fine)
{
	float* v = fine->h_v;
	float* f = fine->h_f;

	float h_x = fine->h_x;
	float h_y = fine->h_y;
	float h_z = fine->h_z;

	float h_x2 = h_x*h_x;
	float h_y2 = h_y*h_y;
	float h_z2 = h_z*h_z;

	int sizeX = fine->sizeX;
	int sizeY = fine->sizeY;
	int sizeZ = fine->sizeZ;

	float* residual = (float*)malloc(sizeX*sizeY*sizeZ*sizeof(float));

	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			for(int posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;
				if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1 || posZ == 0 || posZ == sizeZ - 1)
					residual[idx] = 0.0f; //mi trovo sui bordi
				else
				{
					idx = (posX-1) + posY * sizeX + posZ * sizeX * sizeY;
					float O = v[idx];
					idx = (posX+1) + posY * sizeX + posZ * sizeX * sizeY;
					float E = v[idx];
					idx = posX + (posY-1) * sizeX + posZ * sizeX * sizeY;
					float N = v[idx];
					idx = posX + (posY+1) * sizeX + posZ * sizeX * sizeY;
					float S = v[idx];
					idx = posX + posY * sizeX + (posZ-1) * sizeX * sizeY;
					float D = v[idx];
					idx = posX + posY * sizeX + (posZ+1) * sizeX * sizeY;
					float U = v[idx];

					//int idx = posX + posY * sizeX + posZ * sizeX * sizeY;
					idx = posX + posY * sizeX + posZ * sizeX * sizeY;
					residual[idx] = f[idx] - ((O-2*v[idx]+E)/h_x2) - ((N-2*v[idx]-S)/h_y2) - ((D-2*v[idx]-U)/h_z2);
				}
			}
		}
	}

	return residual;
}

void MultiGrid3D::PrintGrid(int gridID)
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids3D[gridID]->PrintGrid_v(logfd);
}

void MultiGrid3D::PrintAllGrids_v()
{
	int logfd = open("log/log_v.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		//printf("Grid: %d\n", i);
		grids3D[i]->PrintGrid_v(logfd);
	}
}

void MultiGrid3D::PrintAllGrids_f()
{
	int logfd = open("log/log_f.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);

	for(int i = 0; i < numGrids; i++)
	{
		//printf("Grid: %d\n", i);
		grids3D[i]->PrintGrid_f(logfd);
	}
}

void MultiGrid3D::PrintDiff()
{
	int logfd = open("log/diff.txt",O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
	grids3D[0]->PrintDiff(logfd);
}
