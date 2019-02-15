#ifndef _DUGKS_DECLARATION_H_
#define _DUGKS_DECLARATION_H_

#include "Mesh_2D.h"

//--------------------------------------------------
extern int Faces,Nodes,Cells;

extern Node_2D *NodeArray;
extern Face_2D *FaceArray;
extern Cell_2D *CellArray;
//--------------------------------------------------
//
//--------------------------------------------------
extern int InteriorFaceNum,PeriodicFaceNum,WallFaceNum,P_InletFaceNum,
		   P_OutletFaceNum,	SymmetryFaceNum, P_FarfieldFaceNum, V_InletFaceNum,BoundFaceNum;

extern Face_2D **InteriorFaceA,**WallFaceA,**PeriodicFaceA,**P_InletFaceA,
	           **P_OutletFaceA,**SymmetryFaceA,**P_FarfieldFaceA,**V_InletFaceA,**BoundFaceA;

extern Cell_2D *PeriodicShadowCA, *WallShadowCA,*P_InletShadowCA,*P_OutletShadowCA,
			   *SymmetryShadowCA, *P_FarfieldShadowCA,*V_InletShadowCA;

extern Cell_2D *PeriodicShadowC_NE, *PeriodicShadowC_NW,
			   *PeriodicShadowC_SE, *PeriodicShadowC_SW;

extern Cell_2D ***CarCellArray;

extern double * const xi_u, * const xi_v;
//---------------------------------------------------
extern bool LeftRight,TopBottom;

extern unsigned left,right,top,bottom;

extern int step;

extern double SumRho, SumT, SumEk;
//------------------------function declaration---------------------

//----------DmVn.cpp----------------------------
extern void Update_DVDF_Eq(Cell_2D &cell);

extern void Update_DVDF_Eqh(Face_2D &face);

extern void Update_MacroVar(Cell_2D& cell);

extern void Update_MacroVar_h(Face_2D& face);

extern void Update_DVDF_Eqh(Face_2D &face,int k);

extern void Update_DVDF_h(Face_2D& face,int k);

extern void Update_DVDF_Flux_h(Face_2D &face,int k);

//-----------------------output-------------------------
extern void printSplitLine(char c = ' ');

extern void printErrorLine(char c = ' ');

extern void printErrorMessage(int line,const char *file,const char *func);

//
template<typename T>
inline int MeshIndex(const T &End,const T &Beg)
{
	return (nullptr == End ? 0 : End - Beg + 1);
}
#endif
