#include <iostream>
#include <iomanip>
#include <omp.h>
#include "DUGKSDeclaration.h"
using std::ios;
using std::setiosflags;
using std::setprecision;
using std::cout;
using std::endl;

extern void SetwRBF(double *wRBF,Cell_2D *cellptr,Cell_2D::DVDF Cell_2D::*dvdf,int const k);

extern void SetwRBF(double *wRBF,Cell_2D *cellptr,double MacroQuantity::*var);

extern double valueRBF(double const *wRBF,Cell_2D* cellptr,double const x,double const y);

extern Cell_2D* targetCell(Cell_2D* cellptr);

inline
double VenkatakrishnanExpression(double a,double b)
{
	return (a*a + 2*a*b + Kh3)/(a*a + 2*b*b + a*b + Kh3);
} 
void VenkatakrishnanFluxLimiter(Cell_2D &cell,int const k)
{
	#ifdef _ARK_MOMENTUM_FLIP
	double GradfBPDotDelta[4] = {0},LfBP[4] = {0};
	double MaxfBP = cell.f.BarP[k],MinfBP = cell.f.BarP[k];
	double MinfBPLimiter = 1;
	#endif
//
	// double MaxfBP = -1E+5,MinfBP = -MaxfBP, MaxgBP = -1E+5,MingBP = -MaxgBP;
	#ifdef _ARK_THERMAL_FLIP
	double GradgBPDotDelta[4] = {0},LgBP[4] = {0};
	double MaxgBP = cell.g.BarP[k],MingBP = cell.g.BarP[k];
	double MingBPLimiter = 1;
	#endif
//	
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		if(cell.Cell_C[iFace]->f.BarP[k] > MaxfBP)
			MaxfBP = cell.Cell_C[iFace]->f.BarP[k];
		if(cell.Cell_C[iFace]->f.BarP[k] < MinfBP)
			MinfBP = cell.Cell_C[iFace]->f.BarP[k];
		#endif
		//isothermal
		#ifdef _ARK_THERMAL_FLIP
		if(cell.Cell_C[iFace]->g.BarP[k] > MaxgBP)
			MaxgBP = cell.Cell_C[iFace]->g.BarP[k];
		if(cell.Cell_C[iFace]->g.BarP[k] < MingBP)
			MingBP = cell.Cell_C[iFace]->g.BarP[k];
		#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		GradfBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.f.BarP_x[k]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.f.BarP_y[k];					        
		if(GradfBPDotDelta[iFace] > 0)
		{
			LfBP[iFace]
			=
			VenkatakrishnanExpression(MaxfBP-cell.f.BarP[k],GradfBPDotDelta[iFace]);
		}
		else if(GradfBPDotDelta[iFace] < 0)
		{
			LfBP[iFace]
			=
			VenkatakrishnanExpression(MinfBP-cell.f.BarP[k],GradfBPDotDelta[iFace]);
		}
		else
			LfBP[iFace] = 1;
		#endif
		//isothermal
		#ifdef _ARK_THERMAL_FLIP
		GradgBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.g.BarP_x[k]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.g.BarP_y[k];
		if(GradgBPDotDelta[iFace] > 0)
		{
			LgBP[iFace]
			=  
			VenkatakrishnanExpression(MaxgBP-cell.g.BarP[k],GradgBPDotDelta[iFace]);
		}
		else if(GradgBPDotDelta[iFace] < 0)
		{
			LgBP[iFace]
			=
			VenkatakrishnanExpression(MingBP-cell.g.BarP[k],GradgBPDotDelta[iFace]);
		}
		else
			LgBP[iFace] = 1;
		#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		if(LfBP[iFace] < MinfBPLimiter) MinfBPLimiter = LfBP[iFace];
		#endif
		//isothermal
		#ifdef _ARK_THERMAL_FLIP
		if(LgBP[iFace] < MingBPLimiter) MingBPLimiter = LgBP[iFace];
		#endif
	}
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	cell.fBPLimiter = MinfBPLimiter;
	#endif
	//isothermal
	#ifdef _ARK_THERMAL_FLIP
	cell.gBPLimiter = MingBPLimiter;
	#endif
}
void UW_Interior_DVDF_Bh_Limiter(Face_2D& face,Cell_2D* ptr_C,int const k)
{
	double dx = face.xf - hDt*xi_u[k] - ptr_C->xc;
	double dy = face.yf - hDt*xi_v[k] - ptr_C->yc;
	VenkatakrishnanFluxLimiter(*ptr_C,k);
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = ptr_C->h.BarP[k]
				   + (dx*ptr_C->h.BarP_x[k] + dy*ptr_C->h.BarP_y[k]);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.f.BhDt[k] = ptr_C->f.BarP[k]
	               + ptr_C->fBPLimiter*(dx*ptr_C->f.BarP_x[k] + dy*ptr_C->f.BarP_y[k]);
	#endif
	//isothermal flip
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = ptr_C->g.BarP[k]
				   + ptr_C->gBPLimiter*(dx*ptr_C->g.BarP_x[k] + dy*ptr_C->g.BarP_y[k]);
	#endif
}
void UW_Interior_DVDF_Bh(Face_2D& face,Cell_2D const* ptr_C,int const k)
{
	double dx = face.xf - hDt*xi_u[k] - ptr_C->xc;
	double dy = face.yf - hDt*xi_v[k] - ptr_C->yc;
	//
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = ptr_C->h.BarP[k]
					  + (dx*ptr_C->h.BarP_x[k] + dy*ptr_C->h.BarP_y[k]);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP				  
	face.f.BhDt[k] = ptr_C->f.BarP[k]
					  + (dx*ptr_C->f.BarP_x[k] + dy*ptr_C->f.BarP_y[k]);
	#endif
	//isothermal flip
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = ptr_C->g.BarP[k]
					  + (dx*ptr_C->g.BarP_x[k] + dy*ptr_C->g.BarP_y[k]);
	#endif
}
void UW_Interior_DVDF_Bh(Face_2D& face,int const k)
{
	Cell_2D const *ow = face.owner, *ne = face.neigh; 
	double dxow = face.xf - hDt*xi_u[k] - ow->xc;
	double dyow = face.yf - hDt*xi_v[k] - ow->yc;

	double dxne = face.xf - hDt*xi_u[k] - ne->xc;
	double dyne = face.yf - hDt*xi_v[k] - ne->yc;
	#ifdef _ARK_ALLENCAHN_FLIP
	// face.h.BhDt[k] = 0.5*(face.owner->h.BarP[k]+face.neigh->h.BarP[k]);
	face.h.BhDt[k] = (ow->h.BarP[k] + (dxow*ow->h.BarP_x[k] + dyow*ow->h.BarP_y[k])
	+ ne->h.BarP[k] + (dxne*ne->h.BarP_x[k] + dyne*ne->h.BarP_y[k]))/2;
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP			  
	face.f.BhDt[k] = (ow->f.BarP[k] + (dxow*ow->f.BarP_x[k] + dyow*ow->f.BarP_y[k])
	+ ne->f.BarP[k] + (dxne*ne->f.BarP_x[k] + dyne*ne->f.BarP_y[k]))/2;
	// face.f.BhDt[k] = 0.5*(face.owner->f.BarP[k]+face.neigh->f.BarP[k]);
	#endif
	//isothermal flip
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = 0.5*(face.owner->g.BarP[k] + face.neigh->g.BarP[k]);
	#endif
}

extern double update_DVDF_x(Cell_2D const *cellptr,Cell_2D::DVDF Cell_2D::*dvdf, int const k);
extern double update_DVDF_y(Cell_2D const *cellptr,Cell_2D::DVDF Cell_2D::*dvdf, int const k);
extern double update_DVDF_xx(Cell_2D const *cellptr,Cell_2D::DVDF Cell_2D::*dvdf, int const k);
extern double update_DVDF_yy(Cell_2D const *cellptr,Cell_2D::DVDF Cell_2D::*dvdf, int const k);
extern double update_DVDF_xy(Cell_2D const *cellptr,Cell_2D::DVDF Cell_2D::*dvdf, int const k);

void UW3rd_Interior_DVDF_Bh(Face_2D &face, Cell_2D* cellptr, int const k)
{
	double dx = face.xf - hDt*xi_u[k] - cellptr->xc;
	double dy = face.yf - hDt*xi_v[k] - cellptr->yc;
	Cell_2D const *cellPointer = targetCell(cellptr);

	#ifdef _ARK_ALLENCAHN_FLIP
	double hBP_x = 0.0, hBP_y = 0.0, hBP_xx = 0.0, hBP_yy = 0.0, hBP_xy = 0.0;

	hBP_x  = update_DVDF_x(cellPointer, &Cell_2D::h, k);
	hBP_y  = update_DVDF_y(cellPointer, &Cell_2D::h, k);
	hBP_xx = update_DVDF_xx(cellPointer, &Cell_2D::h, k);
	hBP_yy = update_DVDF_yy(cellPointer, &Cell_2D::h, k);
	hBP_xy = update_DVDF_xy(cellPointer, &Cell_2D::h, k);

	face.h.BhDt[k] = cellPointer->h.BarP[k] + dx*hBP_x + dy*hBP_y
				   + dx*dx*hBP_xx/2 + dy*dy*hBP_yy/2 + dx*dy*hBP_xy;
	#endif

	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	double fBP_x = 0.0, fBP_y = 0.0, fBP_xx = 0.0, fBP_yy = 0.0, fBP_xy = 0.0;

	fBP_x  = update_DVDF_x(cellPointer, &Cell_2D::f, k);
	fBP_y  = update_DVDF_y(cellPointer, &Cell_2D::f, k);
	fBP_xx = update_DVDF_xx(cellPointer, &Cell_2D::f, k);
	fBP_yy = update_DVDF_yy(cellPointer, &Cell_2D::f, k);
	fBP_xy = update_DVDF_xy(cellPointer, &Cell_2D::f, k);

	face.f.BhDt[k] = cellPointer->f.BarP[k] + dx*fBP_x + dy*fBP_y
				   + dx*dx*fBP_xx/2 + dy*dy*fBP_yy/2 + dx*dy*fBP_xy;
	#endif
	//!thermal
	#ifdef _ARK_THERMAL_FLIP
	double gBP_x = 0.0, gBP_y = 0.0, gBP_xx = 0.0, gBP_yy = 0.0, gBP_xy = 0.0;
	#endif


}
void UW3rd_Interior_DVDF_Bh(Face_2D &face, int const k)
{
	double dxOwner = face.xf - hDt*xi_u[k] - face.owner->xc;
	double dyOwner = face.yf - hDt*xi_v[k] - face.owner->yc;

	double dxNeigh = face.xf - hDt*xi_u[k] - face.neigh->xc;
	double dyNeigh = face.yf - hDt*xi_v[k] - face.neigh->yc;

	Cell_2D const *cellptrOwner = face.owner;
	Cell_2D const *cellptrNeigh = targetCell(face.neigh);

	#ifdef _ARK_ALLENCAHN_FLIP
	double 
	hBP_xOwner = 0.0, hBP_yOwner = 0.0, 
	hBP_xxOwner = 0.0, hBP_yyOwner = 0.0, 
	hBP_xyOwner = 0.0;

	double 
	hBP_xNeigh = 0.0, hBP_yNeigh = 0.0, 
	hBP_xxNeigh = 0.0, hBP_yyNeigh = 0.0, 
	hBP_xyNeigh = 0.0;

	hBP_xOwner  = update_DVDF_x(cellptrOwner, &Cell_2D::h, k);
	hBP_yOwner  = update_DVDF_y(cellptrOwner, &Cell_2D::h, k);
	hBP_xxOwner = update_DVDF_xx(cellptrOwner, &Cell_2D::h, k);
	hBP_yyOwner = update_DVDF_yy(cellptrOwner, &Cell_2D::h, k);
	hBP_xyOwner = update_DVDF_xy(cellptrOwner, &Cell_2D::h, k);

	hBP_xNeigh  = update_DVDF_x(cellptrNeigh, &Cell_2D::h, k);
	hBP_yNeigh  = update_DVDF_y(cellptrNeigh, &Cell_2D::h, k);
	hBP_xxNeigh = update_DVDF_xx(cellptrNeigh, &Cell_2D::h, k);
	hBP_yyNeigh = update_DVDF_yy(cellptrNeigh, &Cell_2D::h, k);
	hBP_xyNeigh = update_DVDF_xy(cellptrNeigh, &Cell_2D::h, k);

	face.h.BhDt[k] = cellptrOwner->h.BarP[k]
				   + dxOwner*hBP_xOwner + dyOwner*hBP_yOwner
				   + dxOwner*dxOwner*hBP_xxOwner/2 + dyOwner*dyOwner*hBP_yyOwner/2
				   + dxOwner*dyOwner*hBP_xyOwner;

	face.h.BhDt[k] += cellptrNeigh->h.BarP[k]
				   + dxNeigh*hBP_xNeigh + dyNeigh*hBP_yNeigh
				   + dxNeigh*dxNeigh*hBP_xxNeigh/2 + dyNeigh*dyNeigh*hBP_yyNeigh/2
				   + dxNeigh*dyNeigh*hBP_xyNeigh;

	face.h.BhDt[k] /= 2;
	#endif

	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP

	double 
	fBP_xOwner = 0.0, fBP_yOwner = 0.0, 
	fBP_xxOwner = 0.0, fBP_yyOwner = 0.0, 
	fBP_xyOwner = 0.0;

	double 
	fBP_xNeigh = 0.0, fBP_yNeigh = 0.0, 
	fBP_xxNeigh = 0.0, fBP_yyNeigh = 0.0, 
	fBP_xyNeigh = 0.0;

	fBP_xOwner  = update_DVDF_x(cellptrOwner, &Cell_2D::f, k);
	fBP_yOwner  = update_DVDF_y(cellptrOwner, &Cell_2D::f, k);
	fBP_xxOwner = update_DVDF_xx(cellptrOwner, &Cell_2D::f, k);
	fBP_yyOwner = update_DVDF_yy(cellptrOwner, &Cell_2D::f, k);
	fBP_xyOwner = update_DVDF_xy(cellptrOwner, &Cell_2D::f, k);

	fBP_xNeigh  = update_DVDF_x(cellptrNeigh, &Cell_2D::f, k);
	fBP_yNeigh  = update_DVDF_y(cellptrNeigh, &Cell_2D::f, k);
	fBP_xxNeigh = update_DVDF_xx(cellptrNeigh, &Cell_2D::f, k);
	fBP_yyNeigh = update_DVDF_yy(cellptrNeigh, &Cell_2D::f, k);
	fBP_xyNeigh = update_DVDF_xy(cellptrNeigh, &Cell_2D::f, k);

	face.f.BhDt[k] = cellptrOwner->f.BarP[k]
				   + dxOwner*fBP_xOwner + dyOwner*fBP_yOwner
				   + dxOwner*dxOwner*fBP_xxOwner/2 + dyOwner*dyOwner*fBP_yyOwner/2
				   + dxOwner*dyOwner*fBP_xyOwner;

	face.f.BhDt[k] += cellptrNeigh->f.BarP[k]
				   + dxNeigh*fBP_xNeigh + dyNeigh*fBP_yNeigh
				   + dxNeigh*dxNeigh*fBP_xxNeigh/2 + dyNeigh*dyNeigh*fBP_yyNeigh/2
				   + dxNeigh*dyNeigh*fBP_xyNeigh;

	face.f.BhDt[k] /= 2;
	#endif
	//!thermal
	#ifdef _ARK_THERMAL_FLIP
	double gBP_x = 0.0, gBP_y = 0.0, gBP_xx = 0.0, gBP_yy = 0.0, gBP_xy = 0.0;
	#endif
}

void RBF_Interior_DVDF_Bh(Face_2D& face,Cell_2D* cellptr,int const k)
{
	double dx = face.xf - hDt*xi_u[k] - cellptr->xc;
	double dy = face.yf - hDt*xi_v[k] - cellptr->yc;
	double wRBF[RBF::DIM] = {0.0};

	#ifdef _ARK_ALLENCAHN_FLIP
	SetwRBF(wRBF,cellptr,&Cell_2D::h,k);
	face.h.BhDt[k] = valueRBF(wRBF,cellptr,dx,dy);
	#endif

	#ifdef _ARK_MOMENTUM_FLIP
	SetwRBF(wRBF,cellptr,&Cell_2D::f,k);
	face.f.BhDt[k] = valueRBF(wRBF,cellptr,dx,dy);
	#endif
}
void RBF_Interior_DVDF_Bh(Face_2D& face,int const k)
{
	double dxOwner = face.xf - hDt*xi_u[k] - face.owner->xc;
	double dyOwner = face.yf - hDt*xi_v[k] - face.owner->yc;

	double dxNeigh = face.xf - hDt*xi_u[k] - face.neigh->xc;
	double dyNeigh = face.yf - hDt*xi_v[k] - face.neigh->yc;

	double wRBFOwner[RBF::DIM] = {0.0};
	double wRBFNeigh[RBF::DIM] = {0.0};

	#ifdef _ARK_ALLENCAHN_FLIP
	SetwRBF(wRBFOwner,face.owner,&Cell_2D::h,k);
	SetwRBF(wRBFNeigh,face.neigh,&Cell_2D::h,k);
	face.h.BhDt[k] = (valueRBF(wRBFOwner,face.owner,dxOwner,dyOwner)
				   + valueRBF(wRBFNeigh,face.neigh,dxNeigh,dyNeigh))/2.0;
	#endif

	#ifdef _ARK_MOMENTUM_FLIP
	SetwRBF(wRBFOwner,face.owner,&Cell_2D::f,k);
	SetwRBF(wRBFNeigh,face.neigh,&Cell_2D::f,k);
	face.f.BhDt[k] = (valueRBF(wRBFOwner,face.owner,dxOwner,dyOwner)
				   + valueRBF(wRBFNeigh,face.neigh,dxNeigh,dyNeigh))/2.0;
	#endif
}
void CD_Interior_DVDF_Bh(Face_2D &face,int const k)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	double _hBP_xF = 1000,_hBP_yF = 1000;
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	double _fBP_xF = 1000,_fBP_yF = 1000;
	#endif
	//!isothermal
	#ifdef _ARK_THERMAL_FLIP
	double _gBP_xF = 1000,_gBP_yF = 1000;
	#endif
	//!central difference
	if(0.0 == face.Vx)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		//_hBP_xF = 0.5*(face.owner->h.BarP_x[k] + face.neigh->h.BarP_x[k]);
		_hBP_xF = (0.5*(face.faceCells[3]->h.BarP[k] - face.faceCells[1]->h.BarP[k])
				+ 0.5*(face.faceCells[2]->h.BarP[k] - face.faceCells[0]->h.BarP[k]))
				/face._dx;
		_hBP_yF = (face.owner->h.BarP[k] - face.neigh->h.BarP[k])/face._dy;
		#endif
		//
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		_fBP_xF = (0.5*(face.faceCells[3]->f.BarP[k] - face.faceCells[1]->f.BarP[k])
				+ 0.5*(face.faceCells[2]->f.BarP[k] - face.faceCells[0]->f.BarP[k]))
				/face._dx;
		_fBP_yF = (face.owner->f.BarP[k] - face.neigh->f.BarP[k])/face._dy;
		#endif
		//!isothemal
		#ifdef _ARK_THERMAL_FLIP
		//_gBP_xF =  0.5*(face.owner->g.BarP_x[k] + face.neigh->g.BarP_x[k]);
		_gBP_xF = (0.5*(face.faceCells[3]->g.BarP[k] - face.faceCells[1]->g.BarP[k])
				+ 0.5*(face.faceCells[2]->g.BarP[k] - face.faceCells[0]->g.BarP[k]))
				/face._dx;
		_gBP_yF = (face.owner->g.BarP[k] - face.neigh->g.BarP[k])/face._dy;
		#endif
	}
	else if(0.0 == face.Vy)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		//_hBP_yF = 0.5*(face.owner->h.BarP_y[k] + face.neigh->h.BarP_y[k]);
		_hBP_yF = (0.5*(face.faceCells[3]->h.BarP[k] - face.faceCells[1]->h.BarP[k])
				+ 0.5*(face.faceCells[2]->h.BarP[k] - face.faceCells[0]->h.BarP[k]))
				/face._dy;
		_hBP_xF = (face.owner->h.BarP[k] - face.neigh->h.BarP[k])/face._dx;
		#endif
		//
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		_fBP_yF = (0.5*(face.faceCells[3]->f.BarP[k] - face.faceCells[1]->f.BarP[k])
				+ 0.5*(face.faceCells[2]->f.BarP[k] - face.faceCells[0]->f.BarP[k]))
				/face._dy;
		_fBP_xF = (face.owner->f.BarP[k] - face.neigh->f.BarP[k])/face._dx;
		#endif
		//
		#ifdef _ARK_THERMAL_FLIP
		_gBP_yF = (0.5*(face.faceCells[3]->g.BarP[k] - face.faceCells[1]->g.BarP[k])
				+ 0.5*(face.faceCells[2]->g.BarP[k] - face.faceCells[0]->g.BarP[k]))
				/face._dy;
		_gBP_xF = (face.owner->g.BarP[k] - face.neigh->g.BarP[k])/face._dx;
		#endif
	}
	else
	{
		cout <<"wrong interface"<<endl;
		getchar();
	}
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = 0.5*(face.owner->h.BarP[k] + face.neigh->h.BarP[k])
			- hDt*(_hBP_xF*xi_u[k] + _hBP_yF*xi_v[k]);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.f.BhDt[k] = 0.5*(face.owner->f.BarP[k] + face.neigh->f.BarP[k])
			- hDt*(_fBP_xF*xi_u[k] + _fBP_yF*xi_v[k]);
	#endif
	//!isothemal
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = 0.5*(face.owner->g.BarP[k] + face.neigh->g.BarP[k])
			- hDt*(_gBP_xF*xi_u[k] + _gBP_yF*xi_v[k]);
	#endif	
}