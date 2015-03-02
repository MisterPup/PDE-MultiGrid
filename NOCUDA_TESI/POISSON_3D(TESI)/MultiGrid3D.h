#ifndef MULTIGRID3D_H
#define MULTIGRID3D_H

#include "Grid3D.h"

class MultiGrid3D
{
	public:
		Grid3D** grids3D;
		int numGrids;

		MultiGrid3D(int finestGridSizeXYZ[], float range[]); //finestGridSizeXYZ array di finestGridSize relativi ad ogni asse
		~MultiGrid3D();
		void InitGrids(int finestGridSizeXYZ[], float range[]);

		void Restrict(float* fine, int fsizeXYZ[], float* coarse, int csizeXYZ[]); //fine-->coarse
		void Interpolate(float* fine, int fsizeXYZ[], float* coarse, int csizeXYZ[]); //coarse-->fine
		void Relax(Grid3D* curGrid, int ncycles); //applica Red-Black Gauss Seidel

		void setToValue(float* grid, int sizeXYZ[], float value, bool modifyBoundaries);
		float* CalculateResidual(Grid3D* fine);
		void ApplyCorrection(float* fine,  int fsizeXYZ[], float* error,  int esizeXYZ[]); //Correzione soluzione/errore su griglia più fine

		void VCycle(int gridID, int v1, int v2);
		void FullMultiGridVCycle(int gridID, int v0, int v1, int v2);

		void PrintGrid(int gridID);
		void PrintAllGrids_v();
		void PrintAllGrids_f();
		void PrintResidual(int gridID);
		void PrintDiff();
		//void PrintAllGrids_r();
};
#endif
