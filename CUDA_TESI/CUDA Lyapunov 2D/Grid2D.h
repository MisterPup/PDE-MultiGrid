#ifndef GRID2D_H
#define GRID2D_H

class Grid2D
{
	public:
		int size; //numero di punti lungo un asse (gli altri assi hanno lo stesso numero) = (2^k) + 1

		float h_x; //distanza fra due punti lungo l'asse x
		float h_y; //distanza fra due punti lungo l'asse y

		//estremi intervallo:
		float x_a;
		float x_b;
		float y_a;
		float y_b;

		/* **********CUDA********** */
		float* d_v; //griglia contenente le soluzioni approssimate (DEVICE)
		float* d_f; //griglia contenente il valore esatto di f (DEVICE)
		size_t d_pitchByte; //dimensione in byte che una riga deve avere per mantenere la coalescenza
		int d_pitch; //pitch in numero di elementi (d_pitch = d_pitchByte/sizeof(tipoDiDato)
		/* **********CUDA********** */

		Grid2D(int size, float range[]) ;
		~Grid2D();
		void PrintDiffApproxReal(int diff_fd); //stampa differenza valore approssimato, valore reale, ha senso solo per griglia finest
		void PrintGrid_v(int logfd);
		void PrintGrid_f(int logfd);
		void PrintMeanAbsoluteError();
};
#endif
