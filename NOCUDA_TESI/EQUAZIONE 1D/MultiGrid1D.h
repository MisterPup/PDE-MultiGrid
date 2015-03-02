#ifndef MULTIGRID1D_H
#define MULTIGRID1D_H

#include "Grid1D.h"

class MultiGrid1D
{
	public:
		Grid1D** grids1D;
		int numGrids;

		MultiGrid1D(int finestGridSize, float range[]);
		~MultiGrid1D();
		void InitGrids(int finestGridSize, float range[]);

		void Restrict(float* fine, int fsize, float* coarse, int csize); //fine-->coarse
		void Interpolate(float* fine, int fsize, float* coarse, int csize); //coarse-->fine

		void Relax(Grid1D* curGrid, int ncycles); //applica Red-Black Gauss Seidel
		float* CalculateResidual(Grid1D* fine);
		void ApplyCorrection(float* fine, int fineSize, float* error, int errorSize); //Correzione soluzione/errore su griglia più fine
		void setToValue(float* grid, int sizeX, float value, bool modifyBoundaries);

		void VCycle(int gridID, int v1, int v2);
		void FullMultiGridVCycle(int gridID, int v0, int v1, int v2);

		void PrintDiff();
		void PrintGrid(int gridID);
		void PrintAllGrids_v();
		void PrintAllGrids_f();
};
#endif
