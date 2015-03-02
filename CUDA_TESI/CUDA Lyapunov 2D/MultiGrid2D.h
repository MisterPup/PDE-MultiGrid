#ifndef MULTIGRID2D_H
#define MULTIGRID2D_H

#include "Grid2D.h"

class MultiGrid2D
{
	public:
		Grid2D** grids2D;
		int numGrids;

		float* d_matrixA; //matrice A (DEVICE)
		int sizeX_A;
		int alfa;

		MultiGrid2D(int finestGridSize, float range[], float* _A, int A_sizeX, int alfa);
		~MultiGrid2D();
		void InitGrids(int finestGridSize, float range[]);
		void InitA(float* _A, int A_size, int alfa);
		void Restrict(float* fine, int fsize, int f_pitch, float* coarse, int csize, int c_pitch); //fine-->coarse
		void Interpolate(float* fine, int fsize, int f_pitch, float* coarse, int csize, int c_pitch); //coarse-->fine
		void Relax(Grid2D* curGrid, int ncycles); //applica Red-Black Gauss Seidel		
		float* CalculateResidual(Grid2D* fine);
		void ApplyCorrection(float* fine, int fineSize, int f_pitch, float* error, int errorSize, int e_pitch); //Correzione soluzione/errore su griglia più fine
		void Set(float* v, int size, int pitch, float value, bool modifyBorder); //se true modifica anche i bordi		
		void VCycle(int gridID, int v1, int v2);
		void FullMultiGridVCycle(int gridID, int v0, int v1, int v2);
		void PrintDiff();
		void PrintGrid(int gridID);
		void PrintAllGrids_v();
		void PrintAllGrids_f();
		void PrintMeanAbsoluteError();
};
#endif
