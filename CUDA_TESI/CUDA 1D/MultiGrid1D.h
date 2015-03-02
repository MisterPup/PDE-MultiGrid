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

		void Restrict(float* fine, int fsizeX, float* coarse, int csizeX); //fine-->coarse
		void Interpolate(float* fine, int fsizeX, float* coarse, int csizeX); //coarse-->fine
		void Relax(Grid1D* curGrid, int ncycles); //applica Red-Black Gauss Seidel
		
		float* CalculateResidual(Grid1D* fine);
		void ApplyCorrection(float* fine, int fineSizeX, float* error, int errorSizeX); //Correzione soluzione/errore su griglia più fine
		void Set(float* d_v, int sizeX,  float value, bool modifyBoundaries); //se true modifica anche i bordi
		
		void VCycle(int gridID, int v1, int v2);
		void FullMultiGridVCycle(int gridID, int v0, int v1, int v2);

		void PrintDiff();
		void PrintGrid(int gridID);
		void PrintAllGrids_v();
		void PrintAllGrids_f();
};
#endif
