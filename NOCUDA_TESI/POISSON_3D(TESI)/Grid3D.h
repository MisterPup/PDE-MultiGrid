#ifndef GRID3D_H
#define GRID3D_H

class Grid3D
{
	public:
		float* h_v; //griglia contenente le soluzioni approssimate
		float* h_f; //griglia contenente il valore esatto di f

		//int size; //numero di punti lungo asse X: sizeX=sizeY=sizeZ = (2^k) + 1
		int sizeX; //numero di punti lungo asse X: sizeX=sizeY=sizeZ = (2^k) + 1
		int sizeY; //numero di punti lungo asse Y: sizeX=sizeY=sizeZ = (2^k) + 1
		int sizeZ; //numero di punti lungo asse Z: sizeX=sizeY=sizeZ = (2^k) + 1
		int* sizeXYZ; //{sizeX, sizeY, sizeZ}

		float h_x; //distanza fra due punti lungo l'asse x
		float h_y; //distanza fra due punti lungo l'asse y
		float h_z; //distanza fra due punti lungo l'asse z

		//estremi intervallo:
		float x_a;
		float x_b;
		float y_a;
		float y_b;
		float z_a;
		float z_b;

		Grid3D(int sizeXYZ[], float range[]);
		~Grid3D();

		void InitV();
		void InitF();

		void PrintGrid_v(int logfd);
		void PrintGrid_f(int logfd);
		void PrintDiff(int logfd);
		void PrintResidual(int logfd);
};
#endif
