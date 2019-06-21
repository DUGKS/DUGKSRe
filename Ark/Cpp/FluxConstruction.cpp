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

double value_DVDF_Bh_3rd
(
	Cell_2D const* cellptr, Cell_2D::DVDF const Cell_2D::*dvdf,int const k,double const dx,double const dy
)
{
	double xBh = 0.0;
	xBh = (cellptr->*dvdf).BarP[k]
		+ dx*(cellptr->*dvdf).BarP_x[k] + dy*(cellptr->*dvdf).BarP_y[k]
		+ dx*dx*(cellptr->*dvdf).BarP_xx[k]/2.0
		+ dy*dy*(cellptr->*dvdf).BarP_yy[k]/2.0
		+ dx*dy*(cellptr->*dvdf).BarP_xy[k];
	return xBh;
}
void UW3rd_Interior_DVDF_Bh
(
	Face_2D &face, Cell_2D const* cellptr, int const k
)
{
	double dx = face.xf - ::hDt*xi_u[k] - cellptr->xc;
	double dy = face.yf - ::hDt*xi_v[k] - cellptr->yc;

	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = value_DVDF_Bh_3rd(cellptr,&Cell_2D::h,k,dx,dy);
	#endif

	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.f.BhDt[k] = value_DVDF_Bh_3rd(cellptr,&Cell_2D::f,k,dx,dy);
	#endif
	//!thermal
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = value_DVDF_Bh_3rd(cellptr,&Cell_2D::g,k,dx,dy);
	#endif


}
void UW3rd_Interior_DVDF_Bh(Face_2D &face, int const k)
{
	double dxOwner = face.xf - ::hDt*xi_u[k] - face.owner->xc;
	double dyOwner = face.yf - ::hDt*xi_v[k] - face.owner->yc;

	double dxNeigh = face.xf - ::hDt*xi_u[k] - face.neigh->xc;
	double dyNeigh = face.yf - ::hDt*xi_v[k] - face.neigh->yc;

	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = 
	value_DVDF_Bh_3rd(face.owner,&Cell_2D::h,k,dxOwner,dyOwner);

	face.h.BhDt[k] += 
	value_DVDF_Bh_3rd(face.neigh,&Cell_2D::h,k,dxNeigh,dyNeigh);

	face.h.BhDt[k] /= 2.0;
	#endif

	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.f.BhDt[k] = 
	value_DVDF_Bh_3rd(face.owner,&Cell_2D::f,k,dxOwner,dyOwner);

	face.f.BhDt[k] += 
	value_DVDF_Bh_3rd(face.neigh,&Cell_2D::f,k,dxNeigh,dyNeigh);

	face.f.BhDt[k] /= 2.0;
	#endif
	//!thermal
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = 
	value_DVDF_Bh_3rd(face.owner,&Cell_2D::g,k,dxOwner,dyOwner);

	face.g.BhDt[k] += 
	value_DVDF_Bh_3rd(face.neigh,&Cell_2D::g,k,dxNeigh,dyNeigh);

	face.g.BhDt[k] /= 2.0;
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