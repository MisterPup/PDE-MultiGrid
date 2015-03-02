#ifndef GRID2D_H
#define GRID2D_H

class Grid2D
{
	public:
		float* h_v; //griglia contenente le soluzioni approssimate
		float* h_f; //griglia contenente il valore esatto di f

		int sizeX; //numero di punti lungo asse X (gli altri assi hanno lo stesso numero) = (2^k) + 1
		int sizeY; //numero di punti lungo asse X (gli altri assi hanno lo stesso numero) = (2^k) + 1
		int* sizeXY; //{sizeX, sizeY}

		float h_x; //distanza fra due punti lungo l'asse x
		float h_y; //distanza fra due punti lungo l'asse y

		//estremi intervallo:
		float x_a;
		float x_b;
		float y_a;
		float y_b;

		Grid2D(int sizeXY[], float range[]);
		~Grid2D();

		void InitV();
		void InitF();

		void PrintDiffApproxReal(int diff_fd); //stampa differenza valore approssimato, valore reale, ha senso solo per griglia finest
		void PrintGrid_v(int logfd);
		void PrintGrid_f(int logfd);
		void PrintResidual(int logfd, float* matrixA, int alfa);
};
#endif
