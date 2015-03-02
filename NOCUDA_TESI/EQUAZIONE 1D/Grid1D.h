#ifndef GRID1D_H
#define GRID1D_H

class Grid1D
{
	public:
		float* h_v; //griglia contenente le soluzioni approssimate
		float* h_f; //griglia contenente il valore esatto di f

		int sizeX; //numero di punti lungo asse x = (2^k) + 1
		float h_x; //distanza fra due punti adiacenti

		//estremi intervallo:
		float x_a;
		float x_b;

		Grid1D(int sizeX, float range[]) ;
		~Grid1D();

		void InitV();
		void InitF();

		void PrintDiffApproxReal(int diff_fd); //stampa differenza valore approssimato, valore reale, ha senso solo per griglia finest
		void PrintGrid_v(int logfd);
		void PrintGrid_f(int logfd);
};
#endif
