#include "Grid3D.h"
#include "inclusion.h"

Grid3D::Grid3D(int sizeXYZ[], float range[]) //range = {xa,xb,ya,yb,za,zb} xb>xa yb>ya zb>za
{
	int sizeX = sizeXYZ[0];
	int sizeY = sizeXYZ[1];
	int sizeZ = sizeXYZ[2];

	assert(sizeX == sizeY);
	assert(sizeX == sizeZ);

	assert((sizeX-1)%2 == 0); //sizeX=(2^k) - 1
	this->sizeX = sizeX;

	assert((sizeY-1)%2 == 0); //sizeY=(2^k) - 1
	this->sizeY = sizeY;

	assert((sizeZ-1)%2 == 0); //sizeZ=(2^k) - 1
	this->sizeZ = sizeZ;

	this->sizeXYZ = (int*)malloc(3*sizeof(int));
	this->sizeXYZ[0] = sizeX;
	this->sizeXYZ[1] = sizeY;
	this->sizeXYZ[2] = sizeZ;

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
	
	h_v = (float*)malloc(sizeX*sizeY*sizeZ*sizeof(float));
	h_f = (float*)malloc(sizeX*sizeY*sizeZ*sizeof(float));

	InitV();
	InitF();

}

Grid3D::~Grid3D()
{
	free(h_v);
	free(h_f);
}

void Grid3D::InitV()
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			for(int posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				if(posX == 0 || posX == sizeX - 1 || posY == 0 || posY == sizeY - 1 || posZ == 0 || posZ == sizeZ - 1)
					h_v[idx] = 0.0f;
			}
		}
	}
}

void Grid3D::InitF()
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			for(int posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				float x = x_a + posX*h_x;
				float y = y_a + posY*h_y;
				float z = z_a + posZ*h_z;

				h_f[idx] = -3*PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z);
			}
		}
	}
}

void Grid3D::PrintGrid_v(int logfd)
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			for(int posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				char* log = (char*)malloc(200);
				printf("posY: %d posX: %d posZ: %d value: %f\n", posY, posX, posZ, h_v[idx]);
				sprintf(log,"posY: %d posX: %d posZ: %d value: %f\n", posY, posX, posZ, h_v[idx]);
				write(logfd, log, strlen(log));
			}
		}
	}
}

void Grid3D::PrintGrid_f(int logfd)
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			for(int posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				char* log = (char*)malloc(200);
				//printf("posY: %d posX: %d posZ: %d value: %f\n", posY, posX, posZ, h_f[idx]);
				sprintf(log,"posY: %d posX: %d posZ: %d value: %f\n", posY, posX, posZ, h_f[idx]);
				write(logfd, log, strlen(log));
			}
		}
	}
}

void Grid3D::PrintDiff(int logfd)
{
	for(int posY = 0; posY < sizeY; posY++)
	{
		for(int posX = 0; posX < sizeX; posX++)
		{
			for(int posZ = 0; posZ < sizeZ; posZ++)
			{
				int idx = posX + posY * sizeX + posZ * sizeX * sizeY;

				float x = x_a + posX*h_x;
				float y = y_a + posY*h_y;
				float z = z_a + posZ*h_z;

				float realSol = sin(PI*x)*sin(PI*y)*sin(PI*z);
				float approxSol = h_v[idx];
				float diff = realSol - approxSol;
				char* log = (char*)malloc(200);
				sprintf(log,"posY: %d posX: %d posZ: %d diff: %f\n", posY, posX, posZ, diff);
				write(logfd, log, strlen(log));
			}
		}
	}
}
