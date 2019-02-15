#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include "DUGKSDeclaration.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::ostringstream;
//--------------------------------------------------------
int Faces = 0,Nodes = 0,Cells = 0;

Node_2D *NodeArray = nullptr;
Face_2D *FaceArray = nullptr;
Cell_2D *CellArray = nullptr;
//-------------------------------------------------------
//
//------------------------------------------------------
int InteriorFaceNum = 0, PeriodicFaceNum = 0, 
	WallFaceNum = 0, P_InletFaceNum = 0, P_OutletFaceNum = 0,
	SymmetryFaceNum = 0,P_FarfieldFaceNum = 0,V_InletFaceNum = 0, BoundFaceNum = 0;

Face_2D **InteriorFaceA = nullptr, **WallFaceA = nullptr,**PeriodicFaceA = nullptr,
		**P_InletFaceA = nullptr, **P_OutletFaceA = nullptr,**SymmetryFaceA = nullptr,
		**P_FarfieldFaceA = nullptr,**V_InletFaceA = nullptr,**BoundFaceA = nullptr;

Cell_2D *PeriodicShadowCA = nullptr, *WallShadowCA = nullptr,*P_InletShadowCA = nullptr,
		*P_OutletShadowCA = nullptr,*SymmetryShadowCA = nullptr,*P_FarfieldShadowCA = nullptr,
		*V_InletShadowCA = nullptr;

//
Cell_2D cellCorners[4];

Cell_2D *PeriodicShadowC_NE = cellCorners, *PeriodicShadowC_NW = cellCorners + 1,
	    *PeriodicShadowC_SW = cellCorners + 2, *PeriodicShadowC_SE = cellCorners + 3;

Cell_2D ***CarCellArray = nullptr;

bool LeftRight = false,TopBottom = false;

unsigned left = 0,right = 0,top = 0,bottom = 0;

int step = 0;
double SumRho = 0.0, SumT = 0.0, SumEk = 0.0;
//--------------MeshConstruct.cpp----------------
extern int MeshConstruct(const string& s);
extern int MeshArea();
extern void FacesClassify();
extern void NeighbourCellConstruct();
extern void  Grad_LSMatrix();
//
//--------------MeshCheck.cpp--------------------
extern int MeshCheck();
extern int MeshOutput(const string& s);
//
//---------------GhostCells.cpp------------------
extern int ShadowCellConstruct();
extern void ShadowCellCornerConstruct();
//
//!----------------MeshCartesian.cpp--------------
extern void DiagonalCellConstruct();
extern void CarfaceCellsConstruct();
extern void CarfaceFacesConstruct();
extern void AllocateCarCellArray();
extern void DeallocateCarCellArray();
extern void SetFace_dxdy();
//
//--------------Qmodel.cpp-------------------
extern void DiscreteVelocityAssign();
extern void setXiDotdS();
//--------------main.cpp---------------------
extern void AllocateResource();
extern void DeallocateResource();
extern void AllocateShadowCorners();
//-------------Preprocess.cpp-------------------
extern void TG_Initialization();
extern void TaylorCouetteInitialization();
extern void ShockStructure();
extern void UniformFlow();
extern void Riemann2D();
extern void LidDrivenSquare();
extern void ForceDrivenTG();
extern void unsteadyTaylorGreen();
extern void SCMP_Drop();
extern void SCMP_FlatInterface();
extern void AC_Drop();
extern void AC_RisingBubble();
extern void AC_LayeredPoiseuille();
extern void AC_InterfaceElongation();
extern void AC_Smoothed();
extern void AC_ZalesakDisk();
extern void AC_DiagTrans();
extern void AC_SpinodalDecomposition();
extern void AC_RayleighTaylor();
//----------------Inc_2DSolver.cpp----------------
extern void DUGKS2DSolver();
//--------------Output.cpp------------------------
extern void OutputCase();
extern void Output_Convergence();
extern void SelfCheck();
#ifdef _ZERO_NDEBUG_FLIP
extern void printCorners(Cell_2D *cellptr);
#endif
int main()
{
	string MeshName(_MESHFILE_NAME_ARK);
	SelfCheck();
	#include "CreatMesh.h"
//- - - - - - - - -Mesh- - - - - - - - - - - - -
	// MeshConstruct(MeshName);
	// MeshArea();
	// FacesClassify();
	// ShadowCellConstruct();
	// NeighbourCellConstruct();
	// OutputCase();
	// #ifdef _CARTESIAN_MESH_FLIP
	// SetFace_dxdy();
	// ShadowCellCornerConstruct();
	// DiagonalCellConstruct();
	// CarfaceCellsConstruct();
	// #endif
//-----------------Preprocess--------------------
	AllocateResource();
	Grad_LSMatrix();
	MeshCheck();

//-------------------Initialization-------------------
	#include "ZeroCondition.h"
	OutputCase();
//------------------Solve-------------------
	#ifndef _ZERO_NDEBUG_FLIP
	DUGKS2DSolver();
	#endif
//------------------Afterprocess----------------
	#ifdef _ZERO_NDEBUG_FLIP
	MeshOutput(MeshName);
	getchar();
	#endif
	DeallocateResource();
}
//
// void AllocateInCells(const int& CellNum,Cell_2D* &ptrCell)
// {
// 	if(CellNum == 0) return;
// 	for(int k = 0;k < CellNum;++k)
// 		ptrCell[k].AllocateInCell();
// }
void AllocateResource()
{
	DiscreteVelocityAssign();
//
	// cout <<"Allocating Resources for Faces..."<<endl;
	// for(int n = 0;n < Faces;++n)
	// 	FaceArray[n].AllocateInFace();
	// cout <<"Allocating Resources for Faces Done"<<endl;
//
	printSplitLine();
	cout <<"Calculating xi*n*dS..."<<endl;
	setXiDotdS();
	cout <<"Calculating xi*n*dS Done"<<endl;
	printSplitLine(nl);
//
	// cout <<"Allocating Resources for Cells..."<<endl;
	// AllocateInCells(Cells,CellArray);
	// AllocateInCells(WallFaceNum,WallShadowCA);
	// AllocateInCells(PeriodicFaceNum,PeriodicShadowCA);
	// AllocateInCells(P_InletFaceNum,P_InletShadowCA);
	// AllocateInCells(P_OutletFaceNum,P_OutletShadowCA);
	// AllocateInCells(SymmetryFaceNum,SymmetryShadowCA);
	// AllocateInCells(P_FarfieldFaceNum,P_FarfieldShadowCA);
	// AllocateInCells(V_InletFaceNum,V_InletShadowCA);
	// for(int n = 0;n < PeriodicFaceNum;++n)
	// {
	// 	PeriodicShadowCA[n] = *PeriodicShadowCA[n].ShadowC;
	// }
	printSplitLine();
	cout <<"Resources Allocating..."<<endl;
	#ifdef _CARTESIAN_MESH_FLIP
	//AllocateShadowCorners();
	AllocateCarCellArray();
	#endif
	cout <<"Resources Allocated"<<endl;
	printSplitLine(nl);
}
void AllocateShadowCorners()
{
	*PeriodicShadowC_SW = *(PeriodicShadowC_SW->ShadowC);
	*PeriodicShadowC_SE = *(PeriodicShadowC_SE->ShadowC);
	*PeriodicShadowC_NW = *(PeriodicShadowC_NW->ShadowC);
	*PeriodicShadowC_NE = *(PeriodicShadowC_NE->ShadowC);
}
void DeallocateShadowCorners()
{
	delete PeriodicShadowC_SW;
	delete PeriodicShadowC_SE;
	delete PeriodicShadowC_NW;
	delete PeriodicShadowC_NE;
}
void DeallocateFaces(int BoundFaceNum, Face_2D** &ptrBoundFaceA)
{
	if (0 == BoundFaceNum) return;
	delete[] ptrBoundFaceA;
	ptrBoundFaceA = nullptr;
}
void DeallocateCells(int BoundFaceNum, Cell_2D* &ptrBoundShadowCA)
{
	if (0 == BoundFaceNum) return;
	delete[] ptrBoundShadowCA;
	ptrBoundShadowCA = nullptr;
}
void DeallocateResource()
{
	delete[] xi_u;
	delete[] xi_v;
	delete[] NodeArray; NodeArray = nullptr;
	delete[] FaceArray; FaceArray = nullptr;
	delete[] CellArray; CellArray = nullptr;
	DeallocateFaces(InteriorFaceNum,InteriorFaceA);
	DeallocateFaces(PeriodicFaceNum,PeriodicFaceA);
	DeallocateFaces(WallFaceNum,WallFaceA);
	DeallocateFaces(P_InletFaceNum,P_InletFaceA);
	DeallocateFaces(P_OutletFaceNum,P_OutletFaceA);
	DeallocateFaces(SymmetryFaceNum,SymmetryFaceA);
	DeallocateFaces(P_FarfieldFaceNum,P_FarfieldFaceA);
	DeallocateFaces(V_InletFaceNum,V_InletFaceA);
//
	DeallocateCells(WallFaceNum,WallShadowCA);
	DeallocateCells(P_InletFaceNum,P_InletShadowCA);
	DeallocateCells(P_OutletFaceNum,P_OutletShadowCA);
	DeallocateCells(SymmetryFaceNum,SymmetryShadowCA);
	DeallocateCells(P_FarfieldFaceNum,P_FarfieldShadowCA);
	DeallocateCells(V_InletFaceNum,V_InletShadowCA);
	DeallocateCells(PeriodicFaceNum,PeriodicShadowCA);
	#ifdef _CARTESIAN_MESH_FLIP
	DeallocateCarCellArray();
	#endif
	printSplitLine();
	cout <<"Resources Deallocated."<<endl;
	cout <<"Programme END."<<endl;	 
	printSplitLine(nl);
}