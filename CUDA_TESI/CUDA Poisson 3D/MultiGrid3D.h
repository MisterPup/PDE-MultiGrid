#ifndef MULTIGRID3D_H
#define MULTIGRID3D_H

#include "Grid3D.h"

class MultiGrid3D
{
	public:
		Grid3D** grids3D;
		int numGrids;

		MultiGrid3D(int finestGridSizeXYZ[], float range[]);
		~MultiGrid3D();
		void InitGrids(int finestGridSizeXYZ[], float range[]);

		void Restrict(float* d_fine, int d_fsizeXYZ[], float* d_coarse, int d_csizeXYZ[]); //fine-->coarse
		void Interpolate(float* d_fine, int d_fsizeXYZ[], float* d_coarse, int d_csizeXYZ[]); //coarse-->fine
		
		void Relax(Grid3D* curGrid, int ncycles); //applica Red-Black Gauss Seidel

		void Set(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder); //se true modifica anche i bordi
		void SetTESTTEST(float* d_v, int d_sizeXYZ[], float value, bool modifyBorder);
		
		float* CalculateResidual(Grid3D* curGrid);				
		void ApplyCorrection(float* fine, int d_fsizeXYZ[], float* error, int d_esizeXYZ[]); //Correzione soluzione/errore su griglia più fine

		void VCycle(int gridID, int v1, int v2);
				
		
		
		void FullMultiGridVCycle(int gridID, int v0, int v1, int v2);
		
		void PrintDiff();
		void PrintGrid(int gridID);
		void PrintAllGrids_v();
		void PrintAllGrids_f();
};
#endif
