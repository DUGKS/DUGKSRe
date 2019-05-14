#include "DUGKSDeclaration.h"

using PhaseFieldAC::PhiL;
using PhaseFieldAC::PhiV;
using PhaseFieldAC::RhoL;
using PhaseFieldAC::RhoV;
using PhaseFieldAC::wI;
using PhaseFieldAC::Gx;
using PhaseFieldAC::Gy;

double const Kp = RT*9/5.0, omega0 = 4.0/9.0;

extern void Grad_Phi_6points(Cell_2D *cellptr);

extern void Grad_Phi_CD(Cell_2D *center);

inline double SourcePhi(double phi)
{
	return -4*(phi-PhiL)*(phi-PhiV)/(wI*(PhiL-PhiV));
}
inline double ChemicalPotential(double Phi,double laplacianPhi)
{
	return
	(
	2*PhaseFieldAC::Beta*(Phi-PhiL)*(Phi-PhiV)*(2*Phi-PhiL-PhiV)
	-PhaseFieldAC::Kappa*laplacianPhi
	);
}
void SourceMomentum(Cell_2D *cellptr)
{
	cellptr->msq->calcRho_xRho_y();
	//
	double F = ChemicalPotential(cellptr->MsQ().Phi,cellptr->MsQ().laplacianPhi);
	cellptr->msq->calcFxFy(F);
	cellptr->msq->Fy += Gy*cellptr->msq->Rho;//Gy*(cellptr->msq->Rho-RhoL);//
	cellptr->msq->Fx += Gx;
	//
	cellptr->msq->U += hDt*cellptr->msq->Fx;
	cellptr->msq->V += hDt*cellptr->msq->Fy;
	cellptr->msq->U /= cellptr->msq->Rho;
	cellptr->msq->V /= cellptr->msq->Rho;
}
void MacroSource(Cell_2D *cellptr)
{
//
	Grad_Phi_6points(cellptr);
//	Grad_Phi_CD(cellptr)
//
	#ifdef _ARK_MOMENTUM_FLIP
	#ifdef _ARK_FORCE_FLIP
	SourceMomentum(cellptr);
	#endif

	cellptr->MsQ().p += (hDt*cellptr->MsQ().RhoXUYV()
			  		- cellptr->MsQ().Rho*omega0*cellptr->MsQ().SqUV()*Lambda0);
	cellptr->MsQ().p *= Kp;

	// cellptr->MsQ().p += ::hDt*cellptr->MsQ().RhoXUYV()*RT;
	#endif

	double L = sqrt(cellptr->MsQ().SqPhixPhiy());
	if(L!= 0)
	{
		double modPhi = SourcePhi(cellptr->MsQ().Phi);

		(cellptr->MsQ().Phi_x) = modPhi*(cellptr->MsQ().Phi_x)/L
								 + cellptr->MsQ().dPhiU()/(RT*dt);

		(cellptr->MsQ().Phi_y) = modPhi*(cellptr->MsQ().Phi_y)/L
								 + cellptr->MsQ().dPhiV()/(RT*dt);
	}
	else
	{
		(cellptr->MsQ().Phi_x) = cellptr->MsQ().dPhiU()/(RT*dt);
		(cellptr->MsQ().Phi_y) = cellptr->MsQ().dPhiV()/(RT*dt);
	}
	//
	cellptr->MsQ().prevPhiU = (cellptr->MsQ().Phi)*(cellptr->MsQ().U);
	cellptr->MsQ().prevPhiV = (cellptr->MsQ().Phi)*(cellptr->MsQ().V);
}