#ifndef _ZERO_INFORMATION_H_
#define _ZERO_INFORMATION_H_

#include "ZeroReference.h"
#include "D2Q9.h"

//------------------------------Normalized Parameters----------------------------
const double 

CFL = -4,

dt = CFL > 0.0 ? CFL*MinL/MaxU: 1.0E-4,

hDt = 0.5*dt;

//---------------------Rankine-Hugoniot-----------------
/*const double

Ma_Outlet = sqrt((Ma*Ma*(Gamma-1.0)+2.0)/(2.0*Gamma*Ma*Ma-(Gamma-1))),

Rho_Outlet = Rho0*(((Gamma+1)*Ma*Ma)/((Gamma-1)*Ma*Ma+2)),

T_Outlet = T0*((1+0.5*(Gamma-1)*Ma*Ma)*(2.0*Gamma/(Gamma-1)*Ma*Ma-1)/
			(Ma*Ma*(2.0*Gamma/(Gamma-1)+0.5*(Gamma-1)))),

U_Outlet = Ma_Outlet*sqrt(Gamma*R0*T_Outlet),

V_Outlet = 0.0;*/

//---------------------2D-Riemann---------------------------------

const double

Rho1 = 0.5313, U1 = 0.0000,V1 = 0.0000,p1 = 0.4,T1 = p1/(Rho1*R0),

Rho2 = 1.0000, U2 = 0.7276,V2 = 0.0000,p2 = 1.0,T2 = p2/(Rho2*R0),

Rho3 = 0.8000, U3 = 0.0000,V3 = 0.0000,p3 = 1.0,T3 = p3/(Rho3*R0),

Rho4 = 1.0000, U4 = 0.0000,V4 = 0.7276,p4 = 1.0,T4 = p4/(Rho4*R0);

//-------------------------Shock Tube---------------------------

const double 

Rho_R = 0.125,

U_R = 0.0,

V_R = 0.0,

T_R = 0.8,

p_R = Rho_R*R0*T_R,

Lambda_R = 1.0/(2*R0*T_R);

namespace RBF
{
	const double c = 35;

	int const DIM = 9;

	void setInverseM();
}

//------------------------multiphase------------------------------
namespace PhaseFieldAC
{
	const double
	
	M_Phi = 0.1,
	
	Cn = 5.0/ChLength,
	
	Pe = U0*ChLength/M_Phi,
	
	wI = ChLength*Cn,

	Sigma = 1E-3,

	Beta = 12.0*Sigma/wI,

	Kappa = 3.0*Sigma*wI/2;

	const double 
	
	PhiL = 1.0,	PhiV = 0.0,
	
	RhoL = 10.0, RhoV = 1.0,

	// MuL = Mu0*RhoL,	MuV = Mu0,

	// NuL  = MuL/RhoL, NuV  = MuV/RhoV,

	MuL = Nu0*RhoL,	MuV = Nu0*RhoV,

	NuL = Nu0, NuV  = Nu0,
	
	aPhi = (RhoL-RhoV)/(PhiL-PhiV), bPhi = RhoV-aPhi*PhiV,

	Uc = 0.0,//1E-4,

	Gx = 0.0,//4.0*Uc*(MuL+MuV)/(Ly*Ly), 

	Gy = 0.0,//-6.25E-6,//0,
	
	TauMass = M_Phi/RT,

	centerX = 0.5*ChLength,

	centerY = 0.5*ChLength,
	
	radius = 0.25*ChLength,

	diameter = 2*radius;

	const double 

	At = (RhoL-RhoV)/(RhoL+RhoV),

	Mo = fabs(Gy)*(RhoL-RhoV)*MuL*MuL*MuL*MuL/(RhoL*RhoL*Sigma*Sigma*Sigma),

	Eo = fabs(Gy)*(RhoL-RhoV)*4*radius*radius/Sigma,

	//ReMP = sqrt(Gy*RhoL*(RhoL-RhoV)*diameter)*diameter/MuL,	//rising bubble

	ReMP = RhoL*ChLength*sqrt(fabs(Gy*ChLength))/MuL,

	La = 0,//Sigma*radius/Nu0/Nu0,

	// ReMP = sqrt(Gy*ChLength)*ChLength/Nu0,

	T = sqrt(ChLength/At/fabs(Gy))/::dt;//(ChLength/U0)/::dt;

	int const 

	iT = 3.0*T;
}

namespace PseudoPotentialSC
{
	double const

	RhoL = 2.4806E-1, RhoV = 4.5435E-2;

	double const 

	K = 1,

	Tr = 0.85,

	wI = 5,

	radius = 0.25*ChLength;
}

namespace D2Q9
{
	double const Alpha = 0.8;
}
/*const int

I_TimeEnd = log(2.0)/(8.0*PI*PI*nu*dt);//TimeEnd/dt;

const double

TimeEnd = I_TimeEnd*dt;*/
//------------------------------Taylor-Couette--------------------------------
const double 

TC_R = 1.5,

TC_r = 0.5,

W_o = 0.0,

W_i = 0.2;

//-----------------------------Output-------------------
const int

VelocityZone = 7,//7 == TC

End_Step = 0,//PhaseFieldAC::iT + 10000,//log(2.0)/(8.0*PI*PI*Nu0*dt),

ZeroDebugControl = 100, //

ConvergenceControl = 100, //SumRho,SumT,independent

ResidualControl = 1000, //print to screen

writeFileControl = 1000; //always >= ResidualControl

double const

RESIDUAL = 1E-6;

#endif
