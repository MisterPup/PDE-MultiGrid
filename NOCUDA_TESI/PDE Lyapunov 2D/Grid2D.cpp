#include "Grid2D.h"
#include "inclusion.h"

Grid2D::Grid2D(int sizeXY[], float range[]) //range = {xa,xb,ya,yb} xb>xa yb>ya
{
	int sizeX = sizeXY[0];
	int sizeY = sizeXY[1];

	assert(sizeX == sizeY);

	assert((sizeX-1)%2 == 0); //size=(2^k) - 1
	this->sizeX = sizeX;

	assert((sizeY-1)%2 == 0); //size=(2^k) - 1
	this->sizeY = sizeY;

	this->sizeXY = (int*)malloc(2*sizeof(int));
	this->sizeXY[0] = sizeX;
	this->sizeXY[1] = sizeY;

	assert(range[1] > range[0]);
	assert(range[3] > range[2]);

	float x_range = range[1] - range[0];
	float y_range = range[3] - range[2];


	x_a = range[0];
	x_b = range[1];
	y_a = range[2];
	y_b = range[3];


	h_x = x_range/(float)(sizeX-1);
	h_y = y_range/(float)(sizeY-1);
	
	h_v = (float*)malloc(sizeX*sizeY*sizeof(float));
	h_f = (float*)malloc(sizeX*sizeY*sizeof(float));

	InitV();
	InitF();
}

Grid2D::~Grid2D()
{
	free(h_v);
	free(h_f);
}

void Grid2D::InitV()
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
				int idx = posX + posY * sizeX;
				if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1)
				{
					float yi = y_a + posY*h_y;
					float xj = x_a + posX*h_x;
					float sol = 2*xj*xj-4*xj*yi+2*yi*yi;
					h_v[idx] = sol;
				}
				else
					h_v[idx] = 0.0f;
		}
	}
}

void Grid2D::InitF()
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			int idx = posX + posY * sizeX;
			h_f[idx] = 0.0f;
		}
	}
}

//ha senso solo per griglia più fine (finest)
void Grid2D::PrintDiffApproxReal(int diff_fd)
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
				int idx = posX + posY * sizeX;

				float xj = x_a + posX*h_x;
				float yi = y_a + posY*h_y;

				float realsol = 2*xj*xj-4*xj*yi+2*yi*yi;
				float approxSol = h_v[idx];
				float diff = approxSol - realsol;

				char* log = (char*)malloc(200);
				sprintf(log,"yi: %f xj: %f diff: %f\n", yi, xj, diff);
				//sprintf(log,"posY: %d posX: %d diff: %f\n", posY, posX, diff);
				write(diff_fd, log, strlen(log));
		}
	}
}

void Grid2D::PrintGrid_v(int logfd)
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
				int idx = posX + posY * sizeX;

				float xj = x_a + posX*h_x;
				float yi = y_a + posY*h_y;
				float realsol = 2*xj*xj-4*xj*yi+2*yi*yi;
				float approxSol = h_v[idx];

				char* log = (char*)malloc(200);
				sprintf(log, "i: %d  j: %d  approxSol: %f  realSol: %f\n", posY, posX, approxSol, realsol);
				write(logfd, log, strlen(log));
		}
	}
}

void Grid2D::PrintGrid_f(int logfd)
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
				int idx = posX + posY * sizeX;

				float xj = x_a + posX*h_x;
				float yi = y_a + posY*h_y;

				char* log = (char*)malloc(200);
				sprintf(log,"yi: %f xj: %f value: %f\n", yi, xj, h_f[idx]);
				//sprintf(log,"posY: %d posX: %d value: %f\n", posY, posX, h_f[idx]);
				write(logfd, log, strlen(log));

		}
	}
}

void Grid2D::PrintResidual(int logfd, float* matrixA, int alfa)
{
	float* residual = (float*)malloc(sizeX*sizeY*sizeof(float));
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			int idx = posX + posY * sizeX;
			//residuo sul contorno è nullo, in quanto lì la soluzione è esatta
			if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1)
			{
				residual[idx] = 0.0f;
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

	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			int idx = posX + posY * sizeX;
			char* log = (char*)malloc(200);
			sprintf(log, "posY: %d posX: %d residual:%f\n", posY, posX, residual[idx]);
			write(logfd, log, strlen(log));
		}
	}
}
