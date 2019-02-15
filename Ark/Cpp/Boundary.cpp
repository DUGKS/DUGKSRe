#include "DUGKSDeclaration.h"
#include <iostream>
using std::cout;
using std::endl;

void P_Inlet_4_Boundary()
{
	LoopPSB(P_InletFaceNum)
	LoopVS(Q)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		P_InletShadowCA[nB].f.BarP[k] = P_InletShadowCA[nB].Cell_C[0]->f.BarP[k];
		#endif
		//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		P_InletShadowCA[nB].g.BarP[k] = P_InletShadowCA[nB].Cell_C[0]->g.BarP[k];
		#endif
	}
}
void P_Outlet_5_Boundary()
{
	LoopPSB(P_OutletFaceNum)
	LoopVS(Q)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		P_OutletShadowCA[nB].f.BarP[k] = P_OutletShadowCA[nB].Cell_C[0]->f.BarP[k];
		#endif
	//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		P_OutletShadowCA[nB].g.BarP[k] = P_OutletShadowCA[nB].Cell_C[0]->g.BarP[k];
		#endif
	}
}
void ExtrapolationfBP(Cell_2D *shadowCell,Cell_2D const *neighb,Cell_2D const *nextNeighb)
{
	// shadowCell->MsQ().Phi = 2*neighb->MsQ().Phi - nextNeighb->MsQ().Phi;
	// shadowCell->MsQ().Rho = 2*neighb->MsQ().Rho - nextNeighb->MsQ().Rho;
	//shadowCell->MsQ() = 2.0*neighb->MsQ() - nextNeighb->MsQ();
	shadowCell->MsQ() = neighb->MsQ();
	shadowCell->MsQ().U = -neighb->MsQ().U;
	shadowCell->MsQ().V = -neighb->MsQ().V;
	shadowCell->MsQ().Lambda = neighb->MsQ().Lambda;
	shadowCell->MsQ().Mu = neighb->MsQ().Mu;
	//!linear extrapolation
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		shadowCell->h.BarP[k] = 2*neighb->h.BarP[k] - nextNeighb->h.BarP[k];
		#endif
		//
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		shadowCell->f.BarP[k] = 2*neighb->f.BarP[k] - nextNeighb->f.BarP[k];
		#endif
		//
		#ifndef _ARK_ISOTHERMAL_FLIP
		shadowCell->g.BarP[k] = 2*neighb->g.BarP[k] - nextNeighb->g.BarP[k];
		#endif
	}
	//!non-equilibrium extrapolation
	LoopVS(Q)
	{

	}
}
void WallShadowC_fBP(Cell_2D &shadowCell)
{
// used for GradfBP
	unsigned const TOP = top,BOTTOM = bottom,LEFT = left,RIGHT = right;
	#ifdef _CARTESIAN_MESH_FLIP
	if(TOP == shadowCell.zone)
	{
		ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[3]);
	}
	else if(BOTTOM == shadowCell.zone)
	{
		ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[1]);
	}
	else if(RIGHT == shadowCell.zone)
	{
		ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[2]);
	}
	else if(LEFT == shadowCell.zone)
	{
		ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[0]);
	}
	#endif
		//
	#ifndef _CARTESIAN_MESH_FLIP
	for(int k = 0;k < Q;++k)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		shadowCell.h.BarP[k] = shadowCell.h.Eq[k]
		+ cell->h.aBP*(cell->h.Tilde[k] - cell->h.Eq[k] + ::hDt*cell.h.So[k]);
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		shadowCell.f.BarP[k] = shadowCell.f.Eq[k]
		+ cell->f.aBP*(cell->f.Tilde[k] - cell->f.Eq[k] + ::hDt*cell.f.So[k]);
		#endif
		//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP	
		shadowCell.g.BarP[k] = shadowCell.g.Eq[k]
		+ cell->g.aBP*(cell->g.Tilde[k] - cell->g.Eq[k] + ::hDt*cell.g.So[k]);
		#endif
	}
	#endif
}
extern void Update_phiFlux_h(Face_2D& face);
//
void Wall_3_BB(Face_2D &face)
{
	using D2Q9::_BB;
	LoopVS(Q)
	{
		if(face.xi_n_dS[k] < 0)
		{
			int k_BB = _BB[k];
			#ifdef _ARK_ALLENCAHN_FLIP
			face.h.hDt[k] = face.h.hDt[k_BB];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			face.f.hDt[k] = face.f.hDt[k_BB];
			#endif
			//!isothemal
			#ifndef _ARK_ISOTHERMAL_FLIP
			face.g.hDt[k] = face.g.hDt[k_BB];
			#endif
		}
	}
}
void Wall_3_NEE(Face_2D &face)
{
	Cell_2D &cell = *face.owner;
	face.MsQh().Phi = cell.MsQ().Phi;
					// + cell.MsQ().Phi_x*(face.xf - cell.xc)
					// + cell.MsQ().Phi_y*(face.yf - cell.yc);
	face.MsQh().Rho = cell.MsQ().Rho;
					// + cell.MsQ().Rho_x*(face.xf - cell.xc)
					// + cell.MsQ().Rho_y*(face.yf - cell.yc);
	face.MsQh().p = cell.MsQ().p;
	// face.MsQh().U = 0.0;
	// face.MsQh().V = 0.0;
	Update_DVDF_Eqh(face);
	//
	#ifdef _ARK_ALLENCAHN_FLIP
	double hNEq = 2.0*cell.h.tau/(2.0*cell.h.tau + ::dt);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	double fNEq = 2.0*cell.f.tau/(2.0*cell.f.tau + ::dt);
	#endif
	//
	#ifndef _ARK_ISOTHERMAL_FLIP	
	double gNEq = 2.0*cell.g.tau/(2.0*cell.g.tau + ::dt);
	#endif
	//
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.hDt[k] = face.h.EqhDt[k]
		+ hNEq*(cell.h.Tilde[k] - cell.h.Eq[k] + hDt*cell.h.So[k]);
		#endif
	//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		face.f.hDt[k] = face.f.EqhDt[k]
		+ fNEq*(cell.f.Tilde[k] - cell.f.Eq[k] + hDt*cell.f.So[k]);
		#endif
	//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP	
		face.g.hDt[k] = face.g.EqhDt[k]
		+ gNEq*(cell.g.Tilde[k] - cell.gEq[k] + hDt*cell.g.So[k]);
		#endif
	}
	Wall_3_BB(face);
}
void fluxCheck(Face_2D const* faceptr)
{
	Face_2D const &face = *faceptr;
	#ifdef _ARK_ALLENCAHN_FLIP
	double hFluxSum = 0.0;
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	double fFluxSum = 0.0;
	#endif
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		hFluxSum += face.h.hDt[k];
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		fFluxSum += face.f.hDt[k];
		#endif
	}
	#ifdef _ARK_ALLENCAHN_FLIP
	if(!EqualZero(hFluxSum))
	{
		cout <<"hfluxSum of Wall Boundary is nonzero"<<endl;
		cout <<"xf : "<<face.xf<<fs<<"yf : "<<face.yf<<fs<<"hFluxSum : "<<hFluxSum<<endl;
		getchar();
	}
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	if(!EqualZero(fFluxSum))
	{
		
		cout <<"ffluxSum of Wall Boundary is nonzero"<<endl;
		cout <<"xf : "<<face.xf<<fs<<"yf : "<<face.yf<<fs<<"fFluxSum : "<<fFluxSum<<endl;
		getchar();
	}
	#endif
}