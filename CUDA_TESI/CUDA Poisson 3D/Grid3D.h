#ifndef GRID3D_H
#define GRID3D_H

class Grid3D
{
	public:
		int sizeX; //numero di punti lungo asseX (gli altri assi hanno lo stesso numero) = (2^k) + 1
		int sizeY;
		int sizeZ;
		/* **********CUDA********** */
		int* d_sizeXYZ;
		/* **********CUDA********** */

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

		/* **********CUDA********** */
		float* d_v; //griglia contenente le soluzioni approssimate (DEVICE)
		float* d_f; //griglia contenente il valore esatto di f (DEVICE)
		//size_t d_pitchByte; //dimensione in byte che una riga deve avere per mantenere la coalescenza
		//int d_pitch; //pitch in numero di elementi (d_pitch = d_pitchByte/sizeof(tipoDiDato)
		/* **********CUDA********** */

		Grid3D(int sizeXYZ[], float range[]) ;
		~Grid3D();
		//void InitF();
		void PrintDiffApproxReal(int diff_fd); //stampa differenza valore approssimato, valore reale, ha senso solo per griglia finest
		void PrintGrid_v(int logfd);
		void PrintGrid_f(int logfd);
};
#endif
