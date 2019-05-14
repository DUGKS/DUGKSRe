#include <iostream>
#include <string>
#include "DUGKSDeclaration.h"
#include <random>
//
using std::string;
using std::cout;
using std::endl;

string caseName;

//--------------------------------------------------------------------------------------
//------------------------------------------Initialization------------------------------
void UniformFlow()
{
		caseName = __func__;
		MacroQuantity Initial(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
		LoopPS(Cells)
		{
			CellArray[n].MsQ() = Initial;
			CellArray[n].MsQ().Fx = 0.0;
			CellArray[n].MsQ().Fy = 0.0;
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.tau = Tau0;
			#endif
			CellArray[n].Factor();

			Update_DVDF_Eq(CellArray[n]);

			LoopVS(Q)
			{
				//!momentum
				#ifdef _ARK_MOMENTUM_FLIP
				CellArray[n].f.Tilde[k]  = CellArray[n].f.Eq[k];
				#endif
//isothermal
				#ifdef _ARK_THERMAL_FLIP
				CellArray[n].g.Tilde[k]  = CellArray[n].g.Eq[k];
				#endif
			}
		}
		LoopPS(Faces)
		{
			FaceArray[n].MsQh() = Initial;
			FaceArray[n].MsQh().Fx = 0.0;
			FaceArray[n].MsQh().Fy = 0.0;
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			FaceArray[n].f.tauh = Tau0;
			#endif
			FaceArray[n].Factor();
		}
		#ifdef _P_INLET_4_BCS_FLIP
		LoopPSB(P_InletFaceNum)
		{
			LoopVS(Q)
			{
				//!momentum
				#ifdef _ARK_MOMENTUM_FLIP
				P_InletShadowCA[nB].f.BarP[k] = P_InletShadowCA[nB].Cell_C[0]->f.Eq[k];
				P_InletShadowCA[nB].f.BarP_x[k] = 0;
				P_InletShadowCA[nB].f.BarP_y[k] = 0;
				#endif
				//!isothermal
				#ifdef _ARK_THERMAL_FLIP
				P_InletShadowCA[nB].g.BarP[k] = P_InletShadowCA[nB].Cell_C[0]->g.Eq[k];
				P_InletShadowCA[nB].g.BarP_x[k] = 0;
				P_InletShadowCA[nB].g.BarP_y[k] = 0;
				#endif			
			}
		}
		#endif
		#ifdef _P_OUTLET_5_BCS_FLIP
		LoopPSB(P_OutletFaceNum)
		{
			LoopVS(Q)
			{
				//!momentum
				#ifdef _ARK_MOMENTUM_FLIP
				P_OutletShadowCA[nB].f.BarP[k] = P_OutletShadowCA[nB].Cell_C[0]->f.Eq[k];
				P_OutletShadowCA[nB].f.BarP_x[k] = 0;
				P_OutletShadowCA[nB].f.BarP_y[k] = 0;
				#endif
				//!isothermal
				#ifdef _ARK_THERMAL_FLIP
				P_OutletShadowCA[nB].g.BarP[k] = P_OutletShadowCA[nB].Cell_C[0]->gEq[k];
				P_OutletShadowCA[nB].g.BarP_x[k] = 0;
				P_OutletShadowCA[nB].g.BarP_y[k] = 0;
				#endif					
			}
		}
		#endif
}
void ShockStructure()
{
	caseName = __func__;
	double Mu_R = Mu0*pow(T_R/T0,Omega0),Tau_R = Mu_R/p_R;
//
	MacroQuantity left(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	MacroQuantity right(Rho_R,U_R,V_R,p_R,T_R,Lambda_R,Mu_R);
	LoopPS(Cells)
	{
		if(CellArray[n].xc <= 0)
		{
			CellArray[n].MsQ() = left;
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.tau = Tau0;
			#endif
			//!energy
			#ifdef _ARK_THERMAL_FLIP
			CellArray[n].g.tau = Tau0;
			#endif
		}
		else
		{
			CellArray[n].MsQ() = right;
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.tau = Tau_R;
			#endif
			//!Energy
			#ifdef _ARK_THERMAL_FLIP
			CellArray[n].g.tau = Tau_R;
			#endif
		}		
		CellArray[n].Factor();
		//
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
			//isothermal
			#ifdef _ARK_THERMAL_FLIP
			CellArray[n].g.Tilde[k] = CellArray[n].g.Eq[k];
			#endif
		}
	}
	#ifdef _P_INLET_4_BCS_FLIP
	LoopPSB(P_InletFaceNum)
	{
		P_InletFaceA[nB]->MsQh() = left;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		P_InletFaceA[nB]->f.tauh = Tau0;
		#endif
		//!Energy
		#ifdef _ARK_THERMAL_FLIP
		P_InletFaceA[nB]->g.tauh = Tau0;
		#endif
		//
		P_InletFaceA[nB]->Factor();
		LoopVS(Q)
		{
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			P_InletShadowCA[nB].f.BarP[k] = P_InletFaceA[nB]->owner->f.Eq[k];
			P_InletShadowCA[nB].f.BarP_x[k] = 0;		
			P_InletShadowCA[nB].f.BarP_y[k] = 0;
			#endif
			//!isothermal
			#ifdef _ARK_THERMAL_FLIP
			P_InletShadowCA[nB].g.BarP[k] = P_InletFaceA[nB]->owner->gEq[k];
			P_InletShadowCA[nB].g.BarP_x[k] = 0;
			P_InletShadowCA[nB].g.BarP_y[k] = 0;
			#endif
		}
	}
	#endif
	//
	#ifdef _P_OUTLET_5_BCS_FLIP
	LoopPS(P_OutletFaceNum)
	{
		P_OutletFaceA[nB]->MsQh() = right;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		P_OutletFaceA[nB]->f.tauh = Tau_R;
		#endif
		//!Energy
		#ifdef _ARK_THERMAL_FLIP
		P_OutletFaceA[nB]->g.tauh = Tau_R;
		#endif
		//
		P_OutletFaceA[nB]->Factor();
		LoopVS(Q)
		{
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			P_OutletShadowCA[nB].f.BarP[k] = P_OutletFaceA[nB]->owner->f.Eq[k];
			P_OutletShadowCA[nB].f.BarP_x[k] = 0;
			P_OutletShadowCA[nB].f.BarP_y[k] = 0;
			#endif
			//!isothermal
			#ifdef _ARK_THERMAL_FLIP
			P_OutletShadowCA[nB].g.BarP[k] = P_OutletFaceA[nB]->owner->g.Eq[k];
			P_OutletShadowCA[nB].g.BarP_x[k] = 0;
			P_OutletShadowCA[nB].g.BarP_y[k] = 0;
			#endif
		}
	}
	#endif
}
//
void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p)
{
	/*u = -U0*cos(2.0*PI*x)*sin(2.0*PI*y)*exp(-8.0*PI*PI*nu*t);
	v =  U0*sin(2.0*PI*x)*cos(2.0*PI*y)*exp(-8.0*PI*PI*nu*t);
	p = -0.25*U0*U0*(cos(4.0*PI*x) + cos(4.0*PI*y))*exp(-16.0*PI*PI*nu*t);*/
	u = -U0*cos(2.0*PI*x)*sin(2.0*PI*y)*exp(-8.0*PI*PI*Nu0*t)/(2.0*PI);
	v =  U0*sin(2.0*PI*x)*cos(2.0*PI*y)*exp(-8.0*PI*PI*Nu0*t)/(2.0*PI);
	p = -0.25*U0*U0*(cos(4.0*PI*x) + cos(4.0*PI*y))*exp(-16.0*PI*PI*Nu0*t)/(4.0*PI*PI);
}
//-----------------------------------ForceDrivenTG-------------------------------
void AnalyticalForceDrivenTG(double x,double y,double &u_A, double &v_A,double &p_A)
{
	u_A = U0*sin(2.0*PI*x/ChLength)*sin(2.0*PI*y/ChLength);
	v_A = U0*cos(2.0*PI*x/ChLength)*cos(2.0*PI*y/ChLength);
	p_A = 0.25*Rho0*U0*U0*(cos(4.0*PI*x/ChLength)-cos(4.0*PI*y/ChLength));
}

void ForceDrivenTG()
{
	LoopPS(Cells)
	{	
		AnalyticalForceDrivenTG
		(
			CellArray[n].xc,CellArray[n].yc,
			CellArray[n].MsQ().U,CellArray[n].MsQ().V,CellArray[n].MsQ().p
		);
// CellArray[n].MsQ().U = U0*sin(2.0*PI*CellArray[n].xc)*sin(2.0*PI*CellArray[n].yc);
// CellArray[n].MsQ().V = U0*cos(2.0*PI*CellArray[n].xc)*cos(2.0*PI*CellArray[n].yc);
// CellArray[n].MsQ().p = 0.25*Rho0*U0*U0*(cos(4.0*PI*CellArray[n].xc)-cos(4.0*PI*CellArray[n].yc));
		CellArray[n].MsQ().Rho = Rho0 + 2.0*CellArray[n].MsQ().p*Lambda0;
		CellArray[n].MsQ().Fx = 8.0*PI*PI*Nu0*CellArray[n].MsQ().U/ChLength/ChLength;
		CellArray[n].MsQ().Fy = 8.0*PI*PI*Nu0*CellArray[n].MsQ().V/ChLength/ChLength;
		CellArray[n].MsQ().T = T0;
		CellArray[n].MsQ().Mu = Mu0;
		CellArray[n].MsQ().Lambda = Lambda0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = Tau0;
		#endif
		CellArray[n].Factor();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}
	}
	LoopPSB(PeriodicFaceNum)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		PeriodicShadowCA[nB].f.tau = Tau0;
		#endif
		PeriodicShadowCA[nB].Factor();
	}
	LoopPS(Faces)
	{
		double x = FaceArray[n].xf, y = FaceArray[n].yf;
		FaceArray[n].MsQh().Fx = 8.0*PI*PI*Nu0*U0*sin(2.0*PI*x/ChLength)*sin(2.0*PI*y/ChLength);
		FaceArray[n].MsQh().Fy = 8.0*PI*PI*Nu0*U0*cos(2.0*PI*x/ChLength)*cos(2.0*PI*y/ChLength);
		FaceArray[n].MsQh().Fx /= (ChLength*ChLength);
		FaceArray[n].MsQh().Fy /= (ChLength*ChLength);
		FaceArray[n].MsQh().T = T0;
		FaceArray[n].MsQh().Mu  = Mu0;
		FaceArray[n].MsQh().Lambda = Lambda0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = Tau0;
		#endif
		FaceArray[n].Factor();
	}
}
void SCMP_FlatInterface()
{
	using namespace PseudoPotentialSC;
	MacroQuantity Initial(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		CellArray[n].MsQ() = Initial;
//------------------flat interface-----------------------
		CellArray[n].MsQ().Rho = 		RhoV + 0.5*(RhoL-RhoV)
								*(
								  tanh(0.4*(CellArray[n].yc-0.25*ChLength))
								  -
								  tanh(0.4*(CellArray[n].yc-0.75*ChLength))
								 );
		CellArray[n].MsQ().p = CellArray[n].MsQ().Rho*RT;
		CellArray[n].MsQ().Fx = 0.0;
		CellArray[n].MsQ().Fy = 0.0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = Tau0;
		#endif
		CellArray[n].Factor();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}
	}
	LoopPS(Faces)
	{
		#ifdef _ARK_FORCE_FLIP
		FaceArray[n].MsQh().Fx = 0.5*((FaceArray[n].neigh->msq->Fx) +(FaceArray[n].owner->msq->Fx));
		FaceArray[n].MsQh().Fy = 0.5*((FaceArray[n].neigh->msq->Fy) +(FaceArray[n].owner->msq->Fy));
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = Tau0;
		#endif
		FaceArray[n].Factor();
	}
}
void SCMP_Drop()
{
	caseName = __func__;
	using namespace PseudoPotentialSC;
	MacroQuantity Initial(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		CellArray[n].MsQ() = Initial;
//------------------flat interface-----------------------
		// CellArray[n].MsQ().Rho = 		RhoV + 0.5*(RhoL-RhoV)
		// 						*(
		// 						  tanh(0.4*(CellArray[n].yc-0.25*ChLength))
		// 						  -
		// 						  tanh(0.4*(CellArray[n].yc-0.75*ChLength))
		// 						 );
//------------------Drop with transition region----------------		
		CellArray[n].MsQ().Rho = 0.5*(RhoL + RhoV) - 0.5*(-RhoL + RhoV)*
		 			tanh
		 			(
		 				2*(sqrt
		 					(
								(CellArray[n].xc-0.5*ChLength)*(CellArray[n].xc-0.5*ChLength)
							+	(CellArray[n].yc-0.5*ChLength)*(CellArray[n].yc-0.5*ChLength)
		 					)-radius
		 					)/wI
		 			);
//------------------Drop with transition region--------------
		// CellArray[n].MsQ().Rho = 0.5*(RhoL + RhoV) - 0.5*(RhoL - RhoV)*
		// 					tanh(2*(sqrt(
		// 						(CellArray[n].yc-0.5*ChLength)*(CellArray[n].yc-0.5*ChLength)
		// 					+	(CellArray[n].xc-0.5*ChLength)*(CellArray[n].xc-0.5*ChLength)
		// 								)- 0.25*ChLength
		// 							)/(0.1*ChLength) 
		// 						);

		CellArray[n].MsQ().p = CellArray[n].MsQ().Rho*RT;
		CellArray[n].MsQ().Fx = 0.0;
		CellArray[n].MsQ().Fy = 0.0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = Tau0;
		#endif
		CellArray[n].Factor();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}
	}
	LoopPS(Faces)
	{
		#ifdef _ARK_FORCE_FLIP
		FaceArray[n].MsQh().Fx = 0.5*((FaceArray[n].neigh->msq->Fx) +(FaceArray[n].owner->msq->Fx));
		FaceArray[n].MsQh().Fy = 0.5*((FaceArray[n].neigh->msq->Fy) +(FaceArray[n].owner->msq->Fy));
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = Tau0;
		#endif
		FaceArray[n].Factor();
	}
}
double AnalyticalPhiAC_Drop(double const x,double const y)
{
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::centerX;
	using PhaseFieldAC::centerY;
	using PhaseFieldAC::wI;
	using PhaseFieldAC::radius;

	double phi = 
	 0.5*(PhiL+PhiV) + 0.5*(PhiL - PhiV)*		//bubble
	//0.5*(PhiL+PhiV) + 0.5*(-PhiL + PhiV)*		//drop
	tanh
	(
	(sqrt((x-centerX)*(x-centerX) + (y-centerY)*(y-centerY))-radius)
	*2/wI
	);

	// int slotL = centerX - 10;
	// int slotR = centerX + 10;
	// double phi = 0;
	// if(sqrt((x-centerX)*(x-centerX) + (y-centerY)*(y-centerY)) > radius)
	// {
	// 	phi = PhiV;
	// }
	// else if(x > slotL && x < slotR && y < centerY)
	// {
	// 	phi = PhiV;
	// }
	// else
	// {
	// 	phi = PhiL;
	// }

	return phi;
}
void AC_ZalesakDisk()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::centerX;
	using PhaseFieldAC::centerY;
	using PhaseFieldAC::wI;
	using PhaseFieldAC::radius;
	#endif
	

	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);

	LoopPS(Cells)
	{
		CellArray[n].MsQ() = init;

		#ifdef _ARK_ALLENCAHN_FLIP
		double const &xc = CellArray[n].xc, &yc = CellArray[n].yc;
		CellArray[n].MsQ().Phi = AnalyticalPhiAC_Drop(xc,yc);
		CellArray[n].MsQ().Rho = aPhi*(CellArray[n].MsQ().Phi) + bPhi;
		CellArray[n].MsQ().U = -U0*PI/ChLength*(yc-centerY);
		CellArray[n].MsQ().V = U0*PI/ChLength*(xc-centerX);
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		#endif

		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = CellArray[n].MsQ().calcTau();
		#endif
		CellArray[n].Factor();
		CellArray[n].FactorAC();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}		
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;

		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].MsQh().U = -U0*PI/ChLength*(FaceArray[n].yf-centerY);
		FaceArray[n].MsQh().V = U0*PI/ChLength*(FaceArray[n].xf-centerX);
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = FaceArray[n].MsQh().calcTau();
		#endif
		FaceArray[n].Factor();
		FaceArray[n].FactorAC();
	}
}
void AC_Smoothed()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::wI;
	using PhaseFieldAC::radius;
	#endif
	//
	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		double &xc = CellArray[n].xc, &yc = CellArray[n].yc;
		CellArray[n].MsQ() = init;

		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].MsQ().Phi = AnalyticalPhiAC_Drop(xc,yc);
		CellArray[n].MsQ().Rho = aPhi*(CellArray[n].MsQ().Phi) + bPhi;
		#endif

		CellArray[n].MsQ().U = 
		-U0*sin(4*PI*xc/ChLength)*sin(4*PI*yc/ChLength);
		CellArray[n].MsQ().V = 
		-U0*cos(4*PI*xc/ChLength)*cos(4*PI*yc/ChLength);
		
	//
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = CellArray[n].MsQ().calcTau();
		#endif
		CellArray[n].Factor();
		CellArray[n].FactorAC();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}		
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = FaceArray[n].MsQh().calcTau();
		#endif
		FaceArray[n].Factor();
		FaceArray[n].FactorAC();
	}
}
void AC_InterfaceElongation()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::wI;
	using PhaseFieldAC::radius;
	#endif

	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		double const &xc = CellArray[n].xc, &yc = CellArray[n].yc;
		CellArray[n].MsQ() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].MsQ().Phi = AnalyticalPhiAC_Drop(CellArray[n].xc,CellArray[n].yc);
		CellArray[n].MsQ().Rho = aPhi*(CellArray[n].MsQ().Phi) + bPhi;
		#endif

		// CellArray[n].MsQ().U = 
		// U0*PI*sin(PI*xc/ChLength)*cos(PI*yc/ChLength);
		// CellArray[n].MsQ().V = 
		// -U0*PI*cos(PI*xc/ChLength)*sin(PI*yc/ChLength);

		CellArray[n].MsQ().U = 
		U0*sin(PI*xc/ChLength)*sin(PI*xc/ChLength)*sin(2*PI*yc/ChLength);
		CellArray[n].MsQ().V = 
		-U0*sin(PI*yc/ChLength)*sin(PI*yc/ChLength)*sin(2*PI*xc/ChLength);
		

		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = CellArray[n].MsQ().calcTau();
		#endif
		CellArray[n].Factor();
		CellArray[n].FactorAC();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}		
	}
	LoopPS(Faces)
	{
		double const &xf = FaceArray[n].xf, &yf = FaceArray[n].yf;
		FaceArray[n].MsQh() = init;

		FaceArray[n].MsQh().U = 
		U0*PI*sin(PI*xf/ChLength)*cos(PI*yf/ChLength);
		FaceArray[n].MsQh().V = 
		-U0*PI*cos(PI*xf/ChLength)*sin(PI*yf/ChLength);

		// FaceArray[n].MsQh().U = 
		// U0*sin(PI*xf/ChLength)*sin(PI*xf/ChLength)*sin(2*PI*yf/ChLength);
		// FaceArray[n].MsQh().V =
		// -U0*sin(PI*yf/ChLength)*sin(PI*yf/ChLength)*sin(2*PI*xf/ChLength);

		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = FaceArray[n].MsQh().calcTau();
		#endif
		FaceArray[n].Factor();
		FaceArray[n].FactorAC();
	}
}
void AC_DiagTrans()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::centerX;
	using PhaseFieldAC::centerY;
	using PhaseFieldAC::wI;
	using PhaseFieldAC::radius;
	#endif

	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		CellArray[n].MsQ() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].MsQ().Phi = AnalyticalPhiAC_Drop(CellArray[n].xc,CellArray[n].yc);
		CellArray[n].MsQ().Rho = aPhi*(CellArray[n].MsQ().Phi) + bPhi;
		#endif
		CellArray[n].MsQ().U = U0;
		CellArray[n].MsQ().V = U0;
		

		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = CellArray[n].MsQ().calcTau();
		#endif
		CellArray[n].Factor();
		CellArray[n].FactorAC();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}		
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;
		FaceArray[n].MsQh().U = U0;
		FaceArray[n].MsQh().V = U0;
		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = FaceArray[n].MsQh().calcTau();
		#endif
		FaceArray[n].Factor();
		FaceArray[n].FactorAC();
	}
}
void AC_SpinodalDecomposition()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP

	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::wI;

	#endif

	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		CellArray[n].MsQ() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].MsQ().Phi = 0.6 + 0.001*(random()%10);
		CellArray[n].MsQ().Rho = aPhi*(CellArray[n].MsQ().Phi) + bPhi;
		#endif
		
//
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = CellArray[n].MsQ().calcTau();
		#endif
		CellArray[n].Factor();
		CellArray[n].FactorAC();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}		
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = FaceArray[n].MsQh().calcTau();
		#endif
		FaceArray[n].Factor();
		FaceArray[n].FactorAC();
	}
}
void AC_Drop()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP

	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::wI;
	using PhaseFieldAC::radius;

	#endif
//
	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		CellArray[n].MsQ() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].MsQ().Phi = AnalyticalPhiAC_Drop(CellArray[n].xc,CellArray[n].yc);
		CellArray[n].MsQ().Rho = aPhi*(CellArray[n].MsQ().Phi) + bPhi;
		CellArray[n].MsQ().Fx = 0.0;
		CellArray[n].MsQ().Fy = 0.0;
		#endif
		
//
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = CellArray[n].MsQ().calcTau();
		#endif
		CellArray[n].Factor();
		CellArray[n].FactorAC();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}		
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = FaceArray[n].MsQh().calcTau();
		#endif
		FaceArray[n].Factor();
		FaceArray[n].FactorAC();
	}
}
void AC_RisingBubble()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP

	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::wI;
	using PhaseFieldAC::centerX;
	using PhaseFieldAC::centerY;
	using PhaseFieldAC::radius;

	#endif

	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		CellArray[n].MsQ() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].MsQ().Phi = 0.5*(PhiL+PhiV) + 0.5*(PhiL - PhiV)*
		tanh
		(
			2*(sqrt(
					(CellArray[n].xc-centerX)*(CellArray[n].xc-centerX)
				 +	(CellArray[n].yc-centerY)*(CellArray[n].yc-centerY)
				 )-radius
			)/wI
		);
		CellArray[n].MsQ().Rho = aPhi*(CellArray[n].MsQ().Phi) + bPhi;
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		#else
		cout <<"idiocy!!!"<<endl;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = CellArray[n].MsQ().calcTau();
		#endif
		CellArray[n].Factor();
		CellArray[n].FactorAC();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}	
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = FaceArray[n].MsQh().calcTau();
		#endif
		FaceArray[n].Factor();
		FaceArray[n].FactorAC();
	}
}
void LayeredPoiseuilleAnalytical(double const xyz,double &u_A)
{
	using PhaseFieldAC::MuL;
	using PhaseFieldAC::MuV;
	using PhaseFieldAC::Gx;
	double const H = Ly/2;
	if(xyz > 0)
	{
		u_A = Gx*H*H*(-(xyz/H)*(xyz/H)-xyz/H*(MuV-MuL)/(MuV+MuL)+2*MuV/(MuV+MuL))/(2*MuV);
	}
	else
	{
		u_A = Gx*H*H*(-(xyz/H)*(xyz/H)-xyz/H*(MuV-MuL)/(MuV+MuL)+2*MuL/(MuV+MuL))/(2*MuL);
	}
}
void AC_LayeredPoiseuille()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::RhoL;
	using PhaseFieldAC::RhoV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::wI;
	double const MidY = 0;
	#endif

	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Cells)
	{
		CellArray[n].MsQ() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].MsQ().Phi = 0.5*(PhiL+PhiV) + 0.5*(PhiL-PhiV)*tanh((MidY-2*CellArray[n].yc)/wI);
		CellArray[n].MsQ().Rho = 0.5*(RhoL+RhoV) + 0.5*(RhoL-RhoV)*tanh((MidY-2*CellArray[n].yc)/wI);		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		CellArray[n].FactorAC();
		#else
		cout <<"idiocy!!!"<<endl;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = Tau0;
		#endif
		CellArray[n].Factor();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}	
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		FaceArray[n].FactorAC();
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = Tau0;
		#endif
		FaceArray[n].Factor();
	}
	#ifdef _Wall_3_BCs_FLIP
	LoopPSB(WallFaceNum)
	{
		WallShadowCA[nB].MsQ() = init;
		WallShadowCA[nB].MsQ().U = 0;
		WallShadowCA[nB].MsQ().V = 0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		WallShadowCA[nB].f.tau = Tau0;
		#endif
		WallShadowCA[nB].Factor();
	}
	#endif

}
void AC_RayleighTaylor()
{
	caseName = __func__;

	#ifdef _ARK_ALLENCAHN_FLIP
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;
	using PhaseFieldAC::RhoL;
	using PhaseFieldAC::RhoV;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::bPhi;
	using PhaseFieldAC::wI;
	#endif

	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);

	LoopPS(Cells)
	{
		CellArray[n].MsQ() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		double const xc = CellArray[n].xc, yc = CellArray[n].yc;
		double MidY = 2*ChLength + 0.1*ChLength*cos(2*PI*xc/ChLength);
		CellArray[n].MsQ().Phi = 0.5*(PhiL+PhiV) + 0.5*(PhiL-PhiV)*tanh(2*(yc-MidY)/wI);
		CellArray[n].MsQ().Rho = 0.5*(RhoL+RhoV) + 0.5*(RhoL-RhoV)*tanh(2*(yc-MidY)/wI);		
		CellArray[n].h.tau = PhaseFieldAC::TauMass;
		CellArray[n].FactorAC();
		#else
		cout <<"idiocy!!!"<<endl;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = Tau0;
		#endif
		CellArray[n].Factor();

		Update_DVDF_Eq(CellArray[n]);
		
		LoopVS(Q)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			CellArray[n].h.Tilde[k] = CellArray[n].h.Eq[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
		}	
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh() = init;
		#ifdef _ARK_ALLENCAHN_FLIP
		FaceArray[n].h.tauh = PhaseFieldAC::TauMass;
		FaceArray[n].FactorAC();
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = Tau0;
		#endif
		FaceArray[n].Factor();
	}
	#ifdef _Wall_3_BCs_FLIP
	LoopPSB(WallFaceNum)
	{
		WallShadowCA[nB].MsQ() = init;
		WallShadowCA[nB].MsQ().U = 0;
		WallShadowCA[nB].MsQ().V = 0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		WallShadowCA[nB].f.tau = Tau0;
		#endif
		WallShadowCA[nB].Factor();
	}
	#endif
}
/*void TG_Initialization()
{	
	for(int i = 0;i != Cells;++i)
	{
		double &x = CellArray[i].xc, &y = CellArray[i].yc;
		TaylorGreenVortex(0.0,CellArray[i].xc,CellArray[i].yc,
							  CellArray[i].MsQ().U,CellArray[i].v,CellArray[i].MsQ().p);
		CellArray[i].MsQ().Rho = Rho0 + CellArray[i].MsQ().p/RT;
		Update_fEq(CellArray[i]);
//---------------------------------Initialize_fT-------------------------------
		double grad_ux = U0*sin(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vx = U0*cos(2.0*PI*x)*cos(2.0*PI*y);
//
		double grad_uy = -grad_vx;
		double grad_vy = -grad_ux;

		double grad_ut = nu*4.0*PI*U0*cos(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vt = -nu*4.0*PI*U0*sin(2.0*PI*x)*cos(2.0*PI*y);
		for(int k = 0;k < Q;++k)
		{
			double u1 = xi[k].MsQ().U * CellArray[i].MsQ().U + xi[k].v * CellArray[i].v;
			double A_u = (xi[k].MsQ().U + u1*xi[k].MsQ().U/RT - CellArray[i].MsQ().U);
			double A_v = (xi[k].v + u1*xi[k].v/RT - CellArray[i].v);
			double A = omega[k]*Rho0/RT;
			double fEq_t = A * (A_u*grad_ut + A_v*grad_vt);
			double fEq_x = A * (A_u*grad_ux + A_v*grad_vx) * xi[k].MsQ().U;
			double fEq_y = A * (A_u*grad_uy + A_v*grad_vy) * xi[k].v;
			double f = CellArray[i].f.Eqk] - tau*(fEq_t + fEq_x + fEq_y);
			CellArray[i].f.Tilde[k] = f; //- 0.5*dt*(fEq_t + fEq_x + fEq_y);
		}
	}
}*/
void LidDrivenSquare()
{
	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	for(int n = 0;n < Cells;++n)
	{
		CellArray[n].MsQ() = init;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		CellArray[n].f.tau = Tau0;
		#endif
		CellArray[n].Factor();
		Update_DVDF_Eq(CellArray[n]);
		LoopVS(Q)
		{
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			CellArray[n].f.Tilde[k] = CellArray[n].f.Eq[k];
			#endif
			#ifdef _ARK_THERMAL_FLIP
			CellArray[n].g.Tilde[k] = CellArray[n].g.Eq[k];
			#endif
		}
	}
	LoopPS(Faces)
	{
		FaceArray[n].MsQh()=init;	
		if(VelocityZone == FaceArray[n].zone)
		{
			FaceArray[n].MsQh().U = U0;
			FaceArray[n].MsQh().V = V0;
			FaceArray[n].neigh->MsQ().U = U0;
			FaceArray[n].neigh->MsQ().V = V0;
		}
		else
		{
			FaceArray[n].MsQh().U = 0;
			FaceArray[n].MsQh().V = 0;
			FaceArray[n].neigh->MsQ().U = 0;
			FaceArray[n].neigh->MsQ().V = 0;
		}
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		FaceArray[n].f.tauh = Tau0;
		#endif
		FaceArray[n].Factor();
	}
	LoopPSB(WallFaceNum)
	{
		WallShadowCA[nB].MsQ() = init;
		WallShadowCA[nB].MsQ().U = 0;
		WallShadowCA[nB].MsQ().V = 0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		WallShadowCA[nB].f.tau = Tau0;
		#endif
		WallShadowCA[nB].Factor();
	}
}
void TaylorCouetteAnalytical(double x,double y,double &u_A)
{
	const double eta = TC_r/TC_R;
	const double A = -W_i*eta*eta/(1.0 - eta*eta), B = W_i*TC_r*TC_r/(1.0 - eta*eta);
	double r = sqrt(x*x + y*y);
	u_A = A*r + B/r;	
}
//------------------------------------D2GHn------------------------------
void TaylorCouetteInitialization()
{
	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);
	LoopPS(Faces)
	{
		Face_2D &face = FaceArray[n];
		face.MsQh() = init;
		if(4 == face.zone)
		{
			face.MsQh().U = -W_i*TC_r*face.Vy;
			face.MsQh().V =  W_i*TC_r*face.Vx;
			face.neigh->MsQ().U = face.MsQh().U;
			face.neigh->MsQ().V = face.MsQh().V;
		}
		else
		{
			face.MsQh().U = 0;
			face.MsQh().V = 0;
			face.neigh->MsQ().U = 0;
			face.neigh->MsQ().V = 0;
		}
		face.Factor();
	}
	LoopPS(Cells)
	{
		Cell_2D &cell = CellArray[n];
		cell.MsQ() = init;
		cell.MsQ().U = 0;
		cell.MsQ().V = 0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.tau = Tau0;
		#endif
		cell.Factor();
		Update_DVDF_Eq(cell);
		LoopVS(Q)
		{
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			cell.f.Tilde[k] = cell.f.Eq[k];
			#endif
		}
	}
	LoopPSB(WallFaceNum)
	{
		WallShadowCA[nB].MsQ() = init;
		WallShadowCA[nB].MsQ().U = 0;
		WallShadowCA[nB].MsQ().V = 0;
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		WallShadowCA[nB].f.tau = Tau0;
		#endif
		WallShadowCA[nB].Factor();
	}
}
double RBFTestAnalytic(double const x,double const y)
{
	return Rho0*(cos(4.0*PI*x/ChLength)-cos(4.0*PI*y/ChLength));
	// return Rho0*sin(4.0*PI*x/ChLength)*cos(2.0*PI*y/ChLength);
}
void RBFTestFunc()
{
	caseName = __func__;
	MacroQuantity init(Rho0,U0,V0,p0,T0,Lambda0,Mu0);

	LoopPS(Cells)
	{
		Cell_2D &cell = CellArray[n];
		cell.MsQ().Rho = RBFTestAnalytic(cell.xc,cell.yc);

		#ifdef _ARK_ALLENCAHN_FLIP
		CellArray[n].h.BarP[5] = cell.MsQ().Rho;
		#endif

		LoopVS(Q)
		{
			
		}
	}

}
//-----------------------------------D2Q9------------------------
/*void TaylorCouetteInitialization()
{
	for(int i = 0;i != WallFaceNum;++i)
	{
		Face_2D &face = *WallFaceA[i];
		WallFaceA[i]->rhoh = Rho0;
		WallFaceA[i]->neigh->rho = Rho0;
		if(VelocityZone == WallFaceA[i]->zone)
		{
			WallFaceA[i]->uh = -U0*face.Vy;
			WallFaceA[i]->vh = U0*face.Vx;
			WallFaceA[i]->neigh->u = WallFaceA[i]->uh;
			WallFaceA[i]->neigh->v = WallFaceA[i]->vh;
		}
		else
		{
			WallFaceA[i]->uh = 0.0;
			WallFaceA[i]->vh = 0.0;
			WallFaceA[i]->neigh->u = 0.0;
			WallFaceA[i]->neigh->v = 0.0;
		}
	}
	for(int i = 0;i != Cells;++i)
	{
		CellArray[i].MsQ().Rho = Rho0;
		//CellArray[i].MsQ().U = 0.0;
		//CellArray[i].v = 0.0;
		Update_fEq(CellArray[i]);
		for(int k = 0;k < Q;++k)
			CellArray[i].f.Tilde[k] = CellArray[i].f.Eqk];
	}
}*/
// void Riemann2D()
// {
// 	for(int n = 0;n < Cells;++n)
// 	{
// 		if(CellArray[n].xc > 0.5 && CellArray[n].yc > 0.5)
// 		{
// 			CellArray[n].MsQ().Rho = Rho1;
// 			CellArray[n].MsQ().U = U1;
// 			CellArray[n].MsQ().V = V1;
// 			CellArray[n].MsQ().T = T1;
// 			CellArray[n].MsQ().p = Rho1*R0*T1;
// 			CellArray[n].MsQ().Lambda = 0.5/(T1*R0);
// 			CellArray[n].MsQ().Mu = Mu0;
// 		}
// 		else if(CellArray[n].xc <= 0.5 && CellArray[n].yc > 0.5)
// 		{
// 			CellArray[n].MsQ().Rho = Rho2;
// 			CellArray[n].MsQ().U = U2;
// 			CellArray[n].MsQ().V = V2;
// 			CellArray[n].MsQ().T = T2;
// 			CellArray[n].MsQ().p = Rho2*R0*T2;
// 			CellArray[n].MsQ().Lambda = 0.5/(T2*R0);
// 			CellArray[n].MsQ().Mu = Mu0*pow(T2/T0,Omega0);
// 		}
// 		else if(CellArray[n].xc <= 0.5 && CellArray[n].yc <= 0.5)
// 		{
// 			CellArray[n].MsQ().Rho = Rho3;
// 			CellArray[n].MsQ().U = U3;
// 			CellArray[n].MsQ().V = V3;
// 			CellArray[n].MsQ().T = T3;
// 			CellArray[n].MsQ().p = Rho3*R0*T3;
// 			CellArray[n].MsQ().Lambda = 0.5/(T3*R0);
// 			CellArray[n].MsQ().Mu = Mu0*pow(T3/T0,Omega0);
// 		}
// 		else if(CellArray[n].xc > 0.5 && CellArray[n].yc <= 0.5)
// 		{
// 			CellArray[n].MsQ().Rho = Rho4;
// 			CellArray[n].MsQ().U = U4;
// 			CellArray[n].MsQ().V = V4;
// 			CellArray[n].MsQ().T = T4;
// 			CellArray[n].MsQ().p = Rho4*R0*T4;
// 			CellArray[n].MsQ().Lambda = 0.5/(T4*R0);
// 			CellArray[n].MsQ().Mu = Mu0*pow(T4/T0,Omega0);
// 		}
// 		CellArray[n].MsQ().qx = 0;
// 		CellArray[n].MsQ().qy = 0;
// 		CellArray[n].Factor();
// 		Update_DVDF_Eq(CellArray[n]);
// 		for(int i = 0;i < DV_Qu;++i)
// 		for(int j = 0;j < DV_Qv;++j)
// 		{
// 			CellArray[n].f.Tilde[i][j] = CellArray[n].f.Eq[i][j];
// 			CellArray[n].g.Tilde[i][j] = CellArray[n].g.Eq[i][j];
// 		}
// 	}
// #ifdef _P_INLET_4_BCS_FLIP
// 	LoopPSB(P_InletFaceNum)
// 	{
// 		for(int i = 0;i < DV_Qu;++i)
// 		for(int j = 0;j < DV_Qv;++j)
// 		{
// 			P_InletShadowCA[nB].f.BarP[i][j] = P_InletShadowCA[nB].Cell_C[0]->f.Eq[i][j];
// 			P_InletShadowCA[nB].g.BarP[i][j] = P_InletShadowCA[nB].Cell_C[0]->gEq[i][j];
// 			P_InletShadowCA[nB].f.BarP_x[i][j] = 0;
// 			P_InletShadowCA[nB].g.BarP_x[i][j] = 0;
// 			P_InletShadowCA[nB].f.BarP_y[i][j] = 0;
// 			P_InletShadowCA[nB].g.BarP_y[i][j] = 0;
// 		}
// 	}
// #endif
// #ifdef _P_OUTLET_5_BCS_FLIP
	// LoopPSB(P_OutletFaceNum)
// 	{
// 		for(int i = 0;i < DV_Qu;++i)
// 		for(int j = 0;j < DV_Qv;++j)
// 		{
// 			P_OutletShadowCA[nB].f.BarP[i][j] = P_OutletShadowCA[nB].Cell_C[0]->f.Eq[i][j];
// 			P_OutletShadowCA[nB].g.BarP[i][j] = P_OutletShadowCA[nB].Cell_C[0]->gEq[i][j];
// 			P_OutletShadowCA[nB].f.BarP_x[i][j] = 0;
// 			P_OutletShadowCA[nB].g.BarP_x[i][j] = 0;
// 			P_OutletShadowCA[nB].f.BarP_y[i][j] = 0;
// 			P_OutletShadowCA[nB].g.BarP_y[i][j] = 0;
// 		}
// 	}
// #endif
// }
void SelfCheck()
{
//
	#if defined _ARK_FORCE_FLIP && defined _ARK_STRANGSPLIT_FLIP
	#error "Fatal Error : Source model collision"
	#endif
//
	#if defined _ARK_ALLENCAHN_FLIP && defined _ARK_PSEUDOPSI_FLIP
	#error "Fatal Error : multiphase model collision"
	#endif
//
	#if defined _FLUX_SCHEME_CD_ARK && defined _FLUX_SCHEME_UW_ARK
	#error "Fatal Error : Flux Scheme CD-UW collision"
	#endif

	#if defined _FLUX_SCHEME_CD_ARK && defined _FLUX_SCHEME_RBF_ARK
	#error "Fatal Error : Flux Scheme CD-RBF collision"
	#endif

	#if defined _FLUX_SCHEME_UW_ARK && defined _FLUX_SCHEME_RBF_ARK
	#error "Fatal Error : Flux Scheme UW-RBF collision"
	#endif

	#if defined _GRAD_SCHEME_9Points && defined _GRAD_SCHEME_5Points
	#error "Fatal Error : Laplacian Scheme Collision"
	#endif
//
	#ifdef _Wall_3_BCs_NEE
	string _bc_ark = _BC_ARK;
	if("NEE" != _bc_ark)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"\"NEE\" != _BC_ARK"<<endl;
		getchar();
	}
	#endif
	#ifdef _Wall_3_BCs_DS
	string _bc_ark = _BC_ARK;
	if("DS" != _BC_ARK)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"\"DS\" != _BC_ARK"<<endl;
		getchar();
	}
	#endif
//
	#if defined _FLUX_SCHEME_CD_ARK
	string flux_scheme = _FLUX_SCHEME_ARK;
		if("CD" != flux_scheme)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"\"CD\" != _FLUX_SCHEME_ARK"<<endl;
			getchar();
		}
	#elif defined _FLUX_SCHEME_UW_ARK
	string flux_scheme = _FLUX_SCHEME_ARK;
		if("UW" != flux_scheme)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"\"UW\" != _FLUX_SCHEME_ARK"<<endl;
			getchar();
		}
	#elif defined _FLUX_SCHEME_RBF_ARK
	string flux_scheme = _FLUX_SCHEME_ARK;
	if("RBF" != flux_scheme)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"\"RBF\" != _FLUX_SCHEME_ARK"<<endl;
		getchar();
	}
	#elif defined _FLUX_SCHEME_UW3RD_ARK
	string flux_scheme = _FLUX_SCHEME_ARK;
	if("UW3rd" != flux_scheme)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"\"RBF\" != _FLUX_SCHEME_ARK"<<endl;
		getchar();
	}
	#endif
	string meshType = _MESHTYPE_ARK;
	if("Quad" == meshType || "Tri" == meshType)
	{
		if("CD" == flux_scheme)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"\"UW\" != _FLUX_SCHEME_ARK"<<endl;
			getchar();
		}
	}
	if(End_Step != 0 && RESIDUAL < 1E-5)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"End_Step != 0 while RESIDUAL < 1E-5"<<endl;
		exit(-1);
	}
	if(End_Step == 0 && RESIDUAL > 1E-5)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"End_Step == 0 while RESIDUAL > 1E-5"<<endl;
		exit(-1);
	}
}