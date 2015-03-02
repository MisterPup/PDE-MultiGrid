#ifndef GRID1D_H
#define GRID1D_H

class Grid1D
{
	public:
		int sizeX; //numero di punti lungo un asse (gli altri assi hanno lo stesso numero) = (2^k) + 1

		float h_x; //distanza fra due punti lungo l'asse x

		//estremi intervallo:
		float x_a;
		float x_b;

		/* **********CUDA********** */
		float* d_v; //griglia contenente le soluzioni approssimate (DEVICE)
		float* d_f; //griglia contenente il valore esatto di f (DEVICE)
		/* **********CUDA********** */

		Grid1D(int sizeX, float range[]) ;
		~Grid1D();
		void PrintDiffApproxReal(int diff_fd); //stampa differenza valore approssimato, valore reale, ha senso solo per griglia finest
		void PrintGrid_v(int logfd);
		void PrintGrid_f(int logfd);
};
#endif
