#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

using PhaseFieldAC::PhiL;
using PhaseFieldAC::PhiV;
using PhaseFieldAC::RhoL;
using PhaseFieldAC::RhoV;
using PhaseFieldAC::aPhi;
using PhaseFieldAC::bPhi;

#include <iostream>
using std::cout;

// namespace D2Q9{

// double const xi_u[Q] = {0,1,1,0,-1,-1,-1,0,1};

// double const xi_v[Q] = {0,0,1,1,1,0,-1,-1,-1};

// int const _BB[Q] = {0,3,4,1,2,7,8,5,6};

// }
extern const char DmQnName[] = "D2Q9AC"; 

extern double * const xi_u = new double[Q];

extern double * const xi_v = new double[Q];

// const double omega[Q]={4.0/9.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0};

double const KForce = 1;

double const Kp = RT*9/5.0;

void MacroQuantity::calcMu()
{
		// Mu = PhaseFieldAC::MuV + 
		// 	(Phi-PhaseFieldAC::PhiV)*(PhaseFieldAC::MuL-PhaseFieldAC::MuV);
		
		// Mu = PhaseFieldAC::MuV*PhaseFieldAC::MuL/
		// 	( 
		// 	  (Phi-PhaseFieldAC::PhiV)*PhaseFieldAC::MuV 
		// 	+ (PhaseFieldAC::PhiL-Phi)*PhaseFieldAC::MuL
		// 	);

		// Mu = PhaseFieldAC::NuV*PhaseFieldAC::NuL/
		// 	( 
		// 	  (Phi-PhaseFieldAC::PhiV)*PhaseFieldAC::NuV 
		// 	+ (PhaseFieldAC::PhiL-Phi)*PhaseFieldAC::NuL
		// 	);
		// Mu *= Rho;

		// Mu = PhaseFieldAC::NuV + 
		// 		(Phi-PhaseFieldAC::PhiV)*(PhaseFieldAC::NuL-PhaseFieldAC::NuV);
		// Mu *= Rho;

		Mu = Rho*RT*Tau0;
}
double MacroQuantity::calcTau()
{
	return Mu/(Rho*RT);
}
void DiscreteVelocityAssign()
{
	LoopVS(Q)
	{
		xi_u[k] = MaSpan*D2Q9::xi_u[k];
		xi_v[k] = MaSpan*D2Q9::xi_v[k];
	}
}
void Update_DVDF_Eq(Cell_2D &cell)
{
	double uu,u1,GAMMA;
	uu = cell.MsQ().SqUV();
	LoopVS(Q)
	{
		u1 = (xi_u[k]*cell.MsQ().U + xi_v[k]*cell.MsQ().V)/RT;
		GAMMA = u1 + 0.5*u1*u1 - uu*Lambda0;
//
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.Eq[k] = D2Q9::omega[k]*(1.0 + u1)*cell.MsQ().Phi;
		cell.h.So[k] = D2Q9::omega[k]*
						(
							cell.MsQ().Phi_x*xi_u[k]
						  + cell.MsQ().Phi_y*xi_v[k]
						);
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.Eq[k] = D2Q9::omega[k]*(cell.MsQ().p/RT + cell.MsQ().Rho*GAMMA)
						 - static_cast<int>(!k)*cell.MsQ().p/RT;
		//
		#ifdef _ARK_FORCE_FLIP

		#ifdef _ARK_FORCE_H_Liang2014
		cell.f.So[k] = 
		D2Q9::omega[k]*
		(
		(1+GAMMA)*((xi_u[k]-cell.MsQ().U)*(cell.MsQ().Fx) + (xi_v[k]-cell.MsQ().V)*(cell.MsQ().Fy))
		/RT
		+
		GAMMA*((xi_u[k]-cell.MsQ().U)*(cell.MsQ().Rho_x) + (xi_v[k]-cell.MsQ().V)*(cell.MsQ().Rho_y))
		);
		#endif

		#ifdef _ARK_FORCE_H_Liang2018
		cell.f.So[k] = 
		D2Q9::omega[k]*
		(
				(xi_u[k]*cell.MsQ().Fx + xi_v[k]*cell.MsQ().Fy)/RT
		+	u1*(xi_u[k]*cell.MsQ().Rho_x + xi_v[k]*cell.MsQ().Rho_y)
		);
		#endif
		//!Force FLip
		#endif
		//!Momentum Flip
		#endif
	}
}
void Update_DVDF_Eqh(Face_2D &face)
{
	double uu,u1,GAMMA;
	uu = face.MsQh().SqUV();
	LoopVS(Q)
	{
		u1 = (xi_u[k]*face.MsQh().U + xi_v[k]*face.MsQh().V)/RT; 
		GAMMA = u1 + 0.5*u1*u1 - uu*Lambda0;

		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.EqhDt[k] = D2Q9::omega[k]*face.MsQh().Phi*(1.0 + u1);
		face.h.SohDt[k] = D2Q9::omega[k]*
							(
								face.MsQh().Phi_x*xi_u[k]
							 	+ 
							 	face.MsQh().Phi_y*xi_v[k]
						    );
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		face.f.EqhDt[k] = D2Q9::omega[k]*(face.MsQh().p/RT + face.MsQh().Rho*GAMMA)
							 - static_cast<int>(!k)*face.MsQh().p/RT;
//
		#ifdef _ARK_FORCE_FLIP

		#ifdef _ARK_FORCE_H_Liang2014
		face.f.SohDt[k] = 
		D2Q9::omega[k]*
		(
		(1+GAMMA)*((xi_u[k]-face.MsQh().U)*face.MsQh().Fx+(xi_v[k]-face.MsQh().V)*(face.MsQh().Fy))
		/RT
		+
		GAMMA*((xi_u[k]-face.MsQh().U)*face.MsQh().Rho_x+(xi_v[k]-face.MsQh().V)*(face.MsQh().Rho_y))
		);
		#endif

		#ifdef _ARK_FORCE_H_Liang2018
		face.f.SohDt[k] = 
		D2Q9::omega[k]*
		(
			(xi_u[k]*face.MsQh().Fx + xi_v[k]*face.MsQh().Fy)/RT
		+	u1*(xi_u[k]*face.MsQh().Rho_x + xi_v[k]*face.MsQh().Rho_y)
		);
		#endif
		//!Force End
		#endif
		//!Momentum End
		#endif
	}
}

void Update_MacroVar_h(Face_2D& face)
{

	face.MsQh().Phi_x = 0.5*((face.owner->msq->Phi_x) + (face.neigh->msq->Phi_x));
	face.MsQh().Phi_y = 0.5*((face.owner->msq->Phi_y) + (face.neigh->msq->Phi_y));
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.MsQh().Rho_x = 0.5*((face.owner->msq->Rho_x) + (face.neigh->msq->Rho_x));
	face.MsQh().Rho_y = 0.5*((face.owner->msq->Rho_y) + (face.neigh->msq->Rho_y));

	face.MsQh().Fx = 0.5*((face.owner->msq->Fx) + (face.neigh->msq->Fx));
	face.MsQh().Fy = 0.5*((face.owner->msq->Fy) + (face.neigh->msq->Fy));
	#endif

	#ifdef _ARK_ALLENCAHN_FLIP
	face.MsQh().Phi  = IntegralGH(Q,face.h.BhDt);
	//resetPhi(face.MsQh().Phi);
	#endif
	//
	face.MsQh().Rho  = aPhi*face.MsQh().Phi + bPhi;
	// double RhoModify = resetRho(face.MsQh().Rho);
	//--------------------------------smoothed flow---------------------------
	// double &xf = face.xf, &yf = face.yf;
	// face.MsQh().U = 
	// -U0*sin(4*PI*xf/ChLength)*sin(4*PI*yf/ChLength)*cos(PI*(step+0.5)*dt/PhaseFieldAC::T);
	// face.MsQh().V = 
	// -U0*cos(4*PI*xf/ChLength)*cos(4*PI*yf/ChLength)*cos(PI*(step+0.5)*dt/PhaseFieldAC::T);
	//--------------------------------decayed shear flow------------------------------
	// double &xf = face.xf, &yf = face.yf;
	// face.MsQh().U = 
	// U0*sin(PI*xf/ChLength)*sin(PI*xf/ChLength)*sin(2*PI*yf/ChLength)
	// *cos(PI*(step+0.5)*dt/PhaseFieldAC::T);
	// face.MsQh().V = 
	// -U0*sin(PI*yf/ChLength)*sin(PI*yf/ChLength)*sin(2*PI*xf/ChLength)
	// *cos(PI*(step+0.5)*dt/PhaseFieldAC::T);
	//---------------------------------non decay---------------------------------------
		// double &xf = face.xf, &yf = face.yf;
	//------------------------------------PRE2014Liang---------------------------
	// if(PhaseFieldAC::iT/2 == step)
	// {
		// face.MsQh().U = 
		// U0*sin(PI*xf/ChLength)*sin(PI*xf/ChLength)*sin(2*PI*yf/ChLength)
		// *cos(PI*(step+0.5)/PhaseFieldAC::iT);

		// face.MsQh().V = 
		// -U0*sin(PI*yf/ChLength)*sin(PI*yf/ChLength)*sin(2*PI*xf/ChLength)
		// *cos(PI*(step+0.5)/PhaseFieldAC::iT);
	// }
	//------------------------------------PRE2016Ren------------------------------
	// if(PhaseFieldAC::iT/2 == step)
	// {
	// 	face.MsQh().U = 
	// 	-U0*PI*sin(PI*xf/ChLength)*cos(PI*yf/ChLength);
	// 	face.MsQh().V = 
	// 	U0*PI*cos(PI*xf/ChLength)*sin(PI*yf/ChLength);
	// }
//
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.MsQh().U    = IntegralGH(Q,face.f.BhDt,xi_u);
	face.MsQh().V    = IntegralGH(Q,face.f.BhDt,xi_v);
	
	#ifdef _ARK_FORCE_FLIP
	face.MsQh().U += 0.5*hDt*face.MsQh().Fx;
	face.MsQh().V += 0.5*hDt*face.MsQh().Fy;
	#endif

	face.MsQh().U /= face.MsQh().Rho;
	face.MsQh().V /= face.MsQh().Rho;

	face.MsQh().p    = IntegralGH(Q,face.f.BhDt)-face.f.BhDt[0];
	face.MsQh().p   += (0.5*hDt*(face.MsQh().RhoXUYV())
				    - D2Q9::omega[0]*face.MsQh().Rho*(face.MsQh().SqUV())*Lambda0);//to be solved +
	face.MsQh().p   *= Kp;
//
	face.MsQh().calcMu();
	face.f.tauh = face.MsQh().calcTau();
	face.Factor();
	#endif
// //
	// face.MsQh().U    = IntegralGH(face.f.BhDt[0],xi_u)/(face.MsQh().Rho*RT);
	// face.MsQh().V    = IntegralGH(face.f.BhDt[0],xi_v)/(face.MsQh().Rho*RT);
	// face.MsQh().p    = IntegralGH(face.f.BhDt[0]);
	// #ifdef _ARK_FORCE_FLIP
	// face.MsQh().U += 0.5*hDt*face.MsQh().Fx/face.MsQh().Rho;
	// face.MsQh().V += 0.5*hDt*face.MsQh().Fy/face.MsQh().Rho;
	// face.MsQh().p += 0.5*hDt*(face.MsQh().U*face.MsQh().Rho_x + face.MsQh().V*face.MsQh().Rho_y);
	// #endif
//-------------------------------------------------------------------

	// face.Mu_h = dynamicViscosity(face.MsQh().Phi);
	// face.Tau_h = face.Mu_h/(face.MsQh().Rho*RT);
	// face.Factor();
	// face.Ms_h   = PhiSource(face.MsQh().Phi);
	//
	// face.MsQh().U    = (IntegralGH(face.f.BhDt[0],xi_u))/Rho0;
	// face.MsQh().V    = (IntegralGH(face.f.BhDt[0],xi_v))/Rho0;
	// face.MsQh().p    = 0.5*(face.MsQh().Rho - Rho0)/Lambda0;
}
void Update_MacroVar(Cell_2D& cell)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	cell.MsQ().Phi  = IntegralGH(Q,cell.h.Tilde);
	#endif
	cell.MsQ().Rho  = aPhi*cell.MsQ().Phi + bPhi;

	//--------------------------------smoothed flow---------------------------
	// double &xc = cell.xc, &yc = cell.yc;
	// cell.MsQ().U = 
	// -U0*sin(4*PI*xc/ChLength)*sin(4*PI*yc/ChLength)*cos(PI*(step+1)*dt/PhaseFieldAC::T);
	// cell.MsQ().V = 
	// -U0*cos(4*PI*xc/ChLength)*cos(4*PI*yc/ChLength)*cos(PI*(step+1)*dt/PhaseFieldAC::T);
	//--------------------------------decayed shear flow-------------------------------
	// double &xc = cell.xc, &yc = cell.yc;
	// cell.MsQ().U = 
	// U0*sin(PI*xc/ChLength)*sin(PI*xc/ChLength)*sin(2*PI*yc/ChLength)
	// *cos(PI*(step+1)*dt/PhaseFieldAC::T);
	// cell.MsQ().V = 
	// -U0*sin(PI*yc/ChLength)*sin(PI*yc/ChLength)*sin(2*PI*xc/ChLength)
	// *cos(PI*(step+1)*dt/PhaseFieldAC::T);
	//------------------Liang2014 interface elongation----------------------
	// double &xc = cell.xc, &yc = cell.yc;

	// cell.MsQ().U = 
	// U0*sin(PI*xc/ChLength)*sin(PI*xc/ChLength)*sin(2*PI*yc/ChLength)
	// *cos(PI*(step+1)/PhaseFieldAC::iT);
	// cell.MsQ().V = 
	// -U0*sin(PI*yc/ChLength)*sin(PI*yc/ChLength)*sin(2*PI*xc/ChLength)
	// *cos(PI*(step+1)/PhaseFieldAC::iT);

	// if(PhaseFieldAC::iT/2 == step)
	// {
	// 	cell.MsQ().U = -U0*PI*sin(PI*xc/ChLength)*cos(PI*yc/ChLength);

	// 	cell.MsQ().V = U0*PI*cos(PI*xc/ChLength)*sin(PI*yc/ChLength);
	// }
//
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	cell.MsQ().U    = IntegralGH(Q,cell.f.Tilde,xi_u);
	cell.MsQ().V    = IntegralGH(Q,cell.f.Tilde,xi_v);
//
	// cell.MsQ().U = 0.0;
	// cell.MsQ().V = 0.0;

	cell.MsQ().p = IntegralGH(Q,cell.f.Tilde) - cell.f.Tilde[0];
//
	cell.MsQ().calcMu();
	cell.f.tau = cell.MsQ().calcTau();
	cell.Factor();
	#endif
//----------------------------------------------------------------
}
//
void IntegralShearStress(){}
void Update_DVDF_Source(Cell_2D&){}
