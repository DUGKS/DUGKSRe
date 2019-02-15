#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include "DUGKSDeclaration.h"

using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;


//index of diagonal cell--index of neighbour cell-------index of the neighbour faces
/*
---------------				----------------				 ------
| D1 |   | D0 |				|    | C1 |    |			     | F1 |     
---------------				----------------			----------------
|    | C |    |				| C2 |  C | C0 |			  F2 |  C | F0  
---------------				----------------			----------------
| D2 |   | D3 |				|    | C3 |    |			     | F3 |     
---------------				----------------				 ------
*/

void sortFacesInThisCell(Cell_2D &cell)
{
	Face_2D *Facetmp[4] = {nullptr,nullptr,nullptr,nullptr};
	for(int i = 0;i < cell.celltype;++i)
	{
		if(fabs(cell.Face_C[i]->xf - cell.xc) > infinitesimal)
		{
			if(cell.Face_C[i]->xf > cell.xc)
			{
				Facetmp[0] = cell.Face_C[i];
			}
			else
			{
				Facetmp[2] = cell.Face_C[i];
			}
		}
		else
		{
			if(cell.Face_C[i]->yf > cell.yc)
			{
				Facetmp[1] = cell.Face_C[i];
			}
			else
			{
				Facetmp[3] = cell.Face_C[i];
			}
		}
	}
	for(int i = 0;i < cell.celltype;++i)
	{
		cell.Face_C[i] = Facetmp[i];
	}
}
void DiagonalCellConstruct()
{
	int ne = 0,nw = 0,se = 0,sw = 0;
	for(int n = 0;n < Cells;++n)
	{
	//---------------------------------Diagonal 0 3------------------------
		if(nullptr == CellArray[n].Cell_C[0]->ShadowC)
		{
			CellArray[n].Cell_Diag[0] = CellArray[n].Cell_C[0]->Cell_C[1];
			CellArray[n].Cell_Diag[3] = CellArray[n].Cell_C[0]->Cell_C[3];
		}
		else
		{
			if(nullptr == CellArray[n].Cell_C[1]->ShadowC)
			{
				CellArray[n].Cell_Diag[0] = CellArray[n].Cell_C[1]->Cell_C[0];
			}
			else
			{
				CellArray[n].Cell_Diag[0] = PeriodicShadowC_NE;
				++ne;
			}
			if(nullptr == CellArray[n].Cell_C[3]->ShadowC)
			{
				CellArray[n].Cell_Diag[3] = CellArray[n].Cell_C[3]->Cell_C[0];
			}
			else
			{
				CellArray[n].Cell_Diag[3] = PeriodicShadowC_SE;
				++se;
			}
		}
	//---------------------------------Diagonal 1 2------------------------
		if(nullptr == CellArray[n].Cell_C[2]->ShadowC)
		{
			CellArray[n].Cell_Diag[1] = CellArray[n].Cell_C[2]->Cell_C[1];
			CellArray[n].Cell_Diag[2] = CellArray[n].Cell_C[2]->Cell_C[3];
		}
		else
		{
			if(nullptr == CellArray[n].Cell_C[1]->ShadowC)
			{
				CellArray[n].Cell_Diag[1] = CellArray[n].Cell_C[1]->Cell_C[2];
			}
			else
			{
				CellArray[n].Cell_Diag[1] = PeriodicShadowC_NW;
				++nw;
			}
			if(nullptr == CellArray[n].Cell_C[3]->ShadowC)
			{
				CellArray[n].Cell_Diag[2] = CellArray[n].Cell_C[3]->Cell_C[2];
			}
			else
			{
				CellArray[n].Cell_Diag[2] = PeriodicShadowC_SW;
				++sw;
			}
		}
	}
}
/*
****               ****
*0 *               * 2*
****               ****
----------face--------->
****       |       ****
*1 *     normal    * 3*
****       |       ****
*/
void CarfaceCellsConstruct()
{
	LoopPS(Faces)
	{
		Face_2D &face = FaceArray[n];
		if(EqualZero(FaceArray[n].Vx))
		{
			if(0 != face._dx)
			{
				cout <<"Collision : Vx == 0 && face dx != 0 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
			face._dx = 2*(face.faceNodes[1]->xN - face.faceNodes[0]->xN);
			if(doubleEqual(FaceArray[n].Vy,1))
			{
				if(face.neigh != face.owner->Cell_C[1])
				{
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[0];
				face.faceCells[1] = face.owner->Cell_Diag[0];
				face.faceCells[2] = face.owner->Cell_C[2];
				face.faceCells[3] = face.owner->Cell_Diag[1];
			}
			else if(doubleEqual(FaceArray[n].Vy,-1))
			{
				if(face.neigh != face.owner->Cell_C[3])
				{
					_PRINT_SPLITLINE_ARK
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					cout <<"face xf : "<<face.xf<<"    "<<"face.yf : "<<face.yf<<endl;
					cout <<"face Vx : "<<face.Vx<<"    "<<"face.Vy : "<<face.Vy<<endl;
					cout <<face.neigh->xc - face.owner->xc <<"    "
						 <<face.neigh->yc - face.owner->yc<<endl;
					cout <<"face neigh : "<<face.neigh<<"    "<<"face.owner : "<<face.owner<<endl;
					for(int k = 0;k < 4;++k)
					{
						cout <<"Cell_C["<<k<<"] : "<<face.owner->Cell_C[k]<<"    ";
					}
					cout << endl;
					_PRINT_ERROR_MSG_FLIP
					_PRINT_SPLITLINE_ARK
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[2];
				face.faceCells[1] = face.owner->Cell_Diag[2];
				face.faceCells[2] = face.owner->Cell_C[0];
				face.faceCells[3] = face.owner->Cell_Diag[3];
			}
			else
			{
				cout <<"Vy != 1 && Vy != -1 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
		}
		else if(EqualZero(FaceArray[n].Vy))
		{
			if(0 != face._dy)
			{
				cout <<"Collision : Vx == 0 && face dx != 0 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
			face._dy = 2*(face.faceNodes[1]->yN - face.faceNodes[0]->yN);
			if(doubleEqual(FaceArray[n].Vx,1))
			{
				if(face.neigh != face.owner->Cell_C[0])
				{
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[3];
				face.faceCells[1] = face.owner->Cell_Diag[3];
				face.faceCells[2] = face.owner->Cell_C[1];
				face.faceCells[3] = face.owner->Cell_Diag[0];
			}
			else if(doubleEqual(FaceArray[n].Vx,-1))
			{
				if(face.neigh != face.owner->Cell_C[2])
				{
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[1];
				face.faceCells[1] = face.owner->Cell_Diag[1];
				face.faceCells[2] = face.owner->Cell_C[3];
				face.faceCells[3] = face.owner->Cell_Diag[2];
			}
			else
			{
				cout <<"Vx != 1 && Vx != -1 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
		}
		else
		{
			cout <<"Attempting to construct faceCells for non Cartesian Mesh "<<endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(-1);
		}
	}
}
//!----------------------construct neighbour Faces in each face---------------
int locateFaceIndex(Face_2D const *const faceptr)
{
	for(int nF = 0;nF < ARK::cellF;++nF)
	{
		if(faceptr == faceptr->owner->Face_C[nF])
			return nF;
	}
	return -1;
}
Cell_2D const* targetCell(Cell_2D const *const cellptr)
{
	if(nullptr == cellptr->ShadowC)
		return cellptr;
	else
		return targetCell(cellptr->ShadowC);
}
void checkIf_NeighMatchOwner
(
	const Face_2D *const faceInowner,
	const Face_2D *const faceInneigh,
	int iFace
)
{
	if(12 == faceInowner->bc_type || 8 ==faceInowner->bc_type)
		return;
	if(faceInowner != faceInneigh)
	{
		printErrorLine();
		_PRINT_ERROR_MSG_FLIP
		cout <<"owner and neigh doesn't share the same face. iFace = "<<iFace<<endl;
		printErrorLine('\n');
		exit(-1);
	}
}
void DiagonalFacesCheck(Face_2D &face,int const k)
{
	if(12 != face.faceFaces[k]->bc_type && 8 != face.faceFaces[k]->bc_type)
	return;
	Face_2D *faceptrA = face.faceFaces[k], *faceptrB = face.faceFaces[k]->shadowF;
	double a = fabs(face.xf - faceptrA->xf) + fabs(face.yf - faceptrA->yf);
	double b = fabs(face.xf - faceptrB->xf) + fabs(face.yf - faceptrB->yf);
	if(a > b)
	{
		face.faceFaces[k] = face.faceFaces[k]->shadowF;
	}
}
void DiagonalfaceFacesConstruct(Face_2D &face,int const iFace)
{
	face.faceFaces[2] = targetCell(face.owner->Cell_Diag[0])->Face_C[iFace];
	face.faceFaces[4] = targetCell(face.owner->Cell_Diag[1])->Face_C[iFace];
	face.faceFaces[6] = targetCell(face.owner->Cell_Diag[2])->Face_C[iFace];
	face.faceFaces[8] = targetCell(face.owner->Cell_Diag[3])->Face_C[iFace];
	DiagonalFacesCheck(face,2);
	DiagonalFacesCheck(face,4);
	DiagonalFacesCheck(face,6);
	DiagonalFacesCheck(face,8);
}
void CarfaceFacesConstruct()
{
	LoopPS(Faces)
	{
		Face_2D &face = FaceArray[n];
		face.faceFaces[0] = &face;
		//!
		int iFace =  locateFaceIndex(&face);
		int iPair = (iFace > 1 ? iFace - 2 : iFace + 2);
		//!
		if(-1 == iFace)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"couldn't find face in its owner cell"<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		//!
		checkIf_NeighMatchOwner(&face,face.neigh->Face_C[iPair],iFace);
		//
		if(0 == iFace)
		{
			face.faceFaces[1] = targetCell(face.neigh)->Face_C[iFace];
			face.faceFaces[3] = targetCell(face.owner->Cell_C[1])->Face_C[iFace];
			face.faceFaces[5] = face.owner->Face_C[iPair];
			face.faceFaces[7] = targetCell(face.owner->Cell_C[3])->Face_C[iFace];	
		}
		else if(1 == iFace)
		{
			face.faceFaces[1] = targetCell(face.owner->Cell_C[0])->Face_C[iFace];
			face.faceFaces[3] = targetCell(face.neigh)->Face_C[iFace];
			face.faceFaces[5] = targetCell(face.owner->Cell_C[2])->Face_C[iFace];
			face.faceFaces[7] = face.owner->Face_C[iPair];
		}
		else if(2 == iFace)
		{
			face.faceFaces[1] = face.owner->Face_C[iPair];
			face.faceFaces[3] = targetCell(face.owner->Cell_C[1])->Face_C[iFace];
			face.faceFaces[5] = targetCell(face.neigh)->Face_C[iFace];
			face.faceFaces[7] = targetCell(face.owner->Cell_C[3])->Face_C[iFace];
		}
		else if(3 == iFace)
		{
			face.faceFaces[1] = targetCell(face.owner->Cell_C[0])->Face_C[iFace];
			face.faceFaces[3] = face.owner->Face_C[iPair];
			face.faceFaces[5] = targetCell(face.owner->Cell_C[2])->Face_C[iFace];
			face.faceFaces[7] = targetCell(face.neigh)->Face_C[iFace];
		}
		//!diagonal faceFaces construct
		DiagonalfaceFacesConstruct(face,iFace);
		//! face type checking	
		if(2 != face.bc_type && 12 != face.bc_type && 8 != face.bc_type)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"wrong face type while constructing faceFaces"<<endl;
			printErrorLine('\n');
			exit(-1);
		}
	}
}
void AllocateCarCellArray()
{
	const double MinDx = dx, MinDy = dy;
//
	CarCellArray = new Cell_2D** [Nxp2];
	for(int i = 0;i < Nxp2;++i)
		CarCellArray[i] = new Cell_2D* [Nyp2];
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		CarCellArray[i][j]=0;
	}
//
	for(int n = 0;n < Cells;++n)
	{
		int i = static_cast<int>(std::round((CellArray[n].xc - X_Beg + 0.5*MinDx)/MinDx));
		int j = static_cast<int>(std::round((CellArray[n].yc - Y_Beg + 0.5*MinDy)/MinDy));
		CarCellArray[i][j] = &CellArray[n];
	}
	#ifdef _PERIODIC_12_8_BCs_FLIP
	for(int k = 0;k < PeriodicFaceNum;++k)
	{
		int i = static_cast<int>(std::round((PeriodicShadowCA[k].xc - X_Beg + 0.5*MinDx)/MinDx));
		int j = static_cast<int>(std::round((PeriodicShadowCA[k].yc - Y_Beg + 0.5*MinDy)/MinDy));
		CarCellArray[i][j] = &PeriodicShadowCA[k];
	}
	#endif
	#ifdef _Wall_3_BCs_FLIP
	for(int k = 0;k < WallFaceNum;++k)
	{
		int i = static_cast<int>(std::round((WallShadowCA[k].xc - X_Beg + 0.5*MinDx)/MinDx));
		int j = static_cast<int>(std::round((WallShadowCA[k].yc - Y_Beg + 0.5*MinDy)/MinDy));
		CarCellArray[i][j] = &WallShadowCA[k];
	}
	#endif
	CarCellArray[0][0] = PeriodicShadowC_SW;
	CarCellArray[0][Nyp1] = PeriodicShadowC_NW;
	CarCellArray[Nxp1][0] = PeriodicShadowC_SE;
	CarCellArray[Nxp1][Nyp1] = PeriodicShadowC_NE;
}
void DeallocateCarCellArray()
{
	for(int i = 0;i < Nxp2;++i)
	{
		delete[] CarCellArray[i];
		CarCellArray[i] = nullptr;
	}
	delete[] CarCellArray;
	CarCellArray = nullptr;
}
void SetFace_dxdy()
{
	for(int n = 0;n < Faces;++n)
	{
		FaceArray[n]._dx = FaceArray[n].owner->xc - FaceArray[n].neigh->xc;
		FaceArray[n]._dy = FaceArray[n].owner->yc - FaceArray[n].neigh->yc; 
		SetZero(FaceArray[n]._dx);
		SetZero(FaceArray[n]._dy);
	}
}