#include <iostream>
#include <iomanip>
#include "DUGKSDeclaration.h"

using std::cout;
using std::endl;
using std::setiosflags;
using std::setprecision;

void checkIfShadowCell(Cell_2D *cellptr)
{
	if(nullptr == cellptr->ShadowC)
	{
		cout <<"fatal error : non shadow cell exits when constructing shadow corner cells";
		_PRINT_ERROR_MSG_FLIP
		getchar();
		exit(0);
	}
}
void ShadowCellAssign(Cell_2D *&lhs,Cell_2D *&rhs)
{
	lhs->ShadowC = rhs;
	lhs->xc = rhs->xc;
	lhs->yc = rhs->yc;
	lhs->volume = rhs->volume;
	lhs->celltype = rhs->celltype;
	*lhs = *rhs;
}
bool locateFace(Face_2D const *faceptr,double const Y,double const X)
{
	return
	(
		fabs(faceptr->yf - Y) < infinitesimal
		&&
		fabs(faceptr->xf - X) < infinitesimal
	);
}
void LeftRightPeriodic()
{
	double Halfdx = ::dx/2.0;
	for(int n = 0;n < BoundFaceNum;++n)
	{
		Cell_2D *target = BoundFaceA[n]->neigh;
		checkIfShadowCell(target);
		if(locateFace(BoundFaceA[n],Y_Beg,X_Beg + Halfdx))
		{
			ShadowCellAssign(PeriodicShadowC_SE,target);
			PeriodicShadowC_SE->xc += Lx;
		}
		else if(locateFace(BoundFaceA[n],Y_Beg,X_End - Halfdx))
		{
			ShadowCellAssign(PeriodicShadowC_SW,target);
			PeriodicShadowC_SW->xc -= Lx;
		}
		else if(locateFace(BoundFaceA[n],Y_End,X_Beg + Halfdx))
		{
			ShadowCellAssign(PeriodicShadowC_NE,target);
			PeriodicShadowC_NE->xc += Lx;
		}
		else if(locateFace(BoundFaceA[n],Y_End,X_End - Halfdx))
		{
			ShadowCellAssign(PeriodicShadowC_NW,target);
			PeriodicShadowC_NW->xc -= Lx;
		}
	}
}
void TopBottomPeriodic()
{
	double Halfdy = dy/2.0;
	for(int n = 0;n < BoundFaceNum;++n)
	{
		Cell_2D *target = BoundFaceA[n]->neigh;
		checkIfShadowCell(target);
		if(locateFace(BoundFaceA[n],Y_Beg+Halfdy,X_Beg))
		{
			ShadowCellAssign(PeriodicShadowC_NW,target);
			PeriodicShadowC_NW->yc += Ly;
		}
		else if(locateFace(BoundFaceA[n],Y_Beg+Halfdy,X_End))
		{
			ShadowCellAssign(PeriodicShadowC_NE,target);
			PeriodicShadowC_NE->yc += Ly;
		}
		else if(locateFace(BoundFaceA[n],Y_End-Halfdy,X_Beg))
		{
			ShadowCellAssign(PeriodicShadowC_SW,target);
			PeriodicShadowC_SW->yc -= Ly;
		}
		else if(locateFace(BoundFaceA[n],Y_End-Halfdy,X_End))
		{
			ShadowCellAssign(PeriodicShadowC_SE,target);
			PeriodicShadowC_SE->yc -= Ly;
		}
	}
}
void ShadowCellCornerConstruct()
{
	if(LeftRight)
	{
		LeftRightPeriodic();
	}
	else if(TopBottom)
	{
		TopBottomPeriodic();
	}
	else
	{
		cout <<"warning : No Periodic Boundary founded"<<endl;
		exit(-1);
	}
}
void lrtbConstruct(Face_2D const* faceptr)
{
	if(X_Beg == faceptr->xf)
	{
		left = faceptr->zone;
	}
	else if(X_End == faceptr->xf)
	{
		right = faceptr->zone;
	}
	else if(Y_Beg == faceptr->yf)
	{
		bottom = faceptr->zone;
	}
	else if(Y_End == faceptr->yf)
	{
		top = faceptr->zone;
	}
}
void ShadowCPeriodicConstruct(int PeriodicFaceNum)
{
	if(0 == PeriodicFaceNum) return;
	PeriodicShadowCA = new Cell_2D[PeriodicFaceNum];
	int k = 0;
	for(int i = 0;i != Faces;++i)
	{
		if(12 == FaceArray[i].bc_type || 8 == FaceArray[i].bc_type)
		{
			FaceArray[i].neigh = PeriodicShadowCA + k;
			PeriodicShadowCA[k].zone = FaceArray[i].zone;
			PeriodicShadowCA[k].Face_C[0] =  FaceArray + i;
			PeriodicShadowCA[k].Cell_C[0] =  FaceArray[i].owner;
			PeriodicShadowCA[k].BoundF = FaceArray + i;
			FaceArray[i].owner->BoundF = FaceArray + i;
			ShadowCellAssign(FaceArray[i].neigh,FaceArray[i].shadowF->owner);
			if(X_Beg == FaceArray[i].xf)
			{
				PeriodicShadowCA[k].xc -= Lx;
				LeftRight = true;
			}
			else if(X_End == FaceArray[i].xf)
			{
				PeriodicShadowCA[k].xc += Lx;
			}
			else if(Y_Beg == FaceArray[i].yf)
			{
				PeriodicShadowCA[k].yc -= Ly;
				TopBottom = true;
			}
			else if(Y_End == FaceArray[i].yf)
			{
				PeriodicShadowCA[k].yc += Ly;
			}
			else
			{
				cout <<"Construct Shadow Cell Failed"<<endl;
				cout <<"xf : "<<FaceArray[i].xf<<"    yf : "<<FaceArray[i].yf<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
			}
			lrtbConstruct(FaceArray + i);
			//PeriodicShadowCA[k] = *PeriodicShadowCA[k].ShadowC;
			++k;
		}
	}
}
void ShadowCBoundConstruct(int BoundFaceNum,Face_2D** &ptrFaceA,Cell_2D* &ptrShadowCA)
{
	if(0 == BoundFaceNum) return;
	ptrShadowCA = new Cell_2D[BoundFaceNum]();
	for(int i = 0;i != BoundFaceNum;++i)
	{
		ptrFaceA[i]->neigh = ptrShadowCA + i;
		ptrShadowCA[i].zone = ptrFaceA[i]->zone;
		ptrShadowCA[i].ShadowC = ptrFaceA[i]->owner;
		ptrShadowCA[i].celltype = ptrShadowCA[i].ShadowC->celltype;
		ptrShadowCA[i].volume = ptrShadowCA[i].ShadowC->volume;
		ptrShadowCA[i].Face_C[0] = ptrFaceA[i];
		//PushCell(*ptrFaceA[i]->owner,ptrShadowCA + i);
		ptrShadowCA[i].Cell_C[0] = ptrFaceA[i]->owner;
		ptrShadowCA[i].BoundF = ptrFaceA[i];
		ptrFaceA[i]->owner->BoundF = ptrFaceA[i];
		if(3 == ptrFaceA[i]->owner->celltype)
		{
			ptrShadowCA[i].xc = ptrFaceA[i]->xf;
			ptrShadowCA[i].yc = ptrFaceA[i]->yf;
		}
		else if(4 == ptrFaceA[i]->owner->celltype)
		{
			ptrShadowCA[i].xc = 2.0*ptrFaceA[i]->xf - ptrShadowCA[i].Cell_C[0]->xc;
			ptrShadowCA[i].yc = 2.0*ptrFaceA[i]->yf - ptrShadowCA[i].Cell_C[0]->yc;	
		}
		else
		{
			cout <<"wrong element type : "<<ptrFaceA[i]->owner->celltype<<endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
		lrtbConstruct(ptrFaceA[i]);
	}
}
int ShadowCellConstruct()
{
	//ShadowCellMacroCheck();
	ShadowCPeriodicConstruct(PeriodicFaceNum);
	ShadowCBoundConstruct(WallFaceNum,WallFaceA,WallShadowCA);
	ShadowCBoundConstruct(P_InletFaceNum,P_InletFaceA,P_InletShadowCA);
	ShadowCBoundConstruct(P_OutletFaceNum,P_OutletFaceA,P_OutletShadowCA);
	ShadowCBoundConstruct(SymmetryFaceNum,SymmetryFaceA,SymmetryShadowCA);
	ShadowCBoundConstruct(P_FarfieldFaceNum,P_FarfieldFaceA,P_FarfieldShadowCA);
	ShadowCBoundConstruct(V_InletFaceNum,V_InletFaceA,V_InletShadowCA);
	return 0;
}