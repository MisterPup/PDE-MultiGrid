#ifndef MULTIGRID2D_H
#define MULTIGRID2D_H

#include "Grid2D.h"

class MultiGrid2D
{
	public:
		Grid2D** grids2D;
		int numGrids;

		float* matrixA;
		int sizeA;
		int alfa;

		MultiGrid2D(int finestGridSizeXY[], float range[], float* _A, int A_size, int alfa);
		~MultiGrid2D();
		void InitGrids(int finestGridSizeXY[], float range[]);
		void InitA(float* _A, int A_size, int alfa);

		void Restrict(float* fine, int fsizeXY[], float* coarse, int csizeXY[]); //fine-->coarse
		void Interpolate(float* fine, int fsizeXY[], float* coarse, int csizeXY[]); //coarse-->fine

		void Relax(Grid2D* curGrid, int ncycles); //applica Red-Black Gauss Seidel
		float* CalculateResidual(Grid2D* fine);
		void ApplyCorrection(float* fine, int fsizeXY[], float* error, int esizeXY[]); //Correzione soluzione/errore su griglia più fine
		void setToValue(float* grid, int sizeXY[], float value, bool modifyBoundaries);

		void VCycle(int gridID, int v1, int v2);
		void FullMultiGridVCycle(int gridID, int v0, int v1, int v2);

		void PrintDiff();
		void PrintGrid(int gridID);
		void PrintAllGrids_v();
		void PrintAllGrids_f();
		void PrintResidual(int gridID);
};
#endif
