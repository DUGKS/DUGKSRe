#ifndef _MACROSCOPIC_QUANTITY_ARK_H
#define _MACROSCOPIC_QUANTITY_ARK_H

#include "ZeroNamespace.h"
// declaration
class MacroQuantity;
//-------------friend arithmetic operators-------------------------
MacroQuantity operator-(const MacroQuantity &lhs,const MacroQuantity &rhs);
MacroQuantity operator+(const MacroQuantity &lhs,const MacroQuantity &rhs);
MacroQuantity operator*(double a,const MacroQuantity &rhs);
MacroQuantity operator*(const MacroQuantity &lhs,double a);
//
class MacroQuantity
{
public:
	MacroQuantity(){};
	MacroQuantity
	(
		double Rho0,double U0,double V0,double p0,double T0,
		double Lambda0,double Mu0
	):Rho(Rho0),U(U0),V(V0),p(p0),T(T0),Lambda(Lambda0),Mu(Mu0)
	{}
	//
	friend MacroQuantity operator-(const MacroQuantity &lhs,const MacroQuantity &rhs);
	friend MacroQuantity operator+(const MacroQuantity &lhs,const MacroQuantity &rhs);
	friend MacroQuantity operator*(double a,const MacroQuantity &rhs);
	friend MacroQuantity operator*(const MacroQuantity &lhs,double a);
	//
	MacroQuantity& operator-=(const MacroQuantity &rhs);
	MacroQuantity& operator+=(const MacroQuantity &rhs);
public:
	double Psi = 0.0;
	double Phi = 0.0,Phi_x = 0.0,Phi_y = 0.0;
	double laplacianPhi = 0.0, prevPhiU = 0.0,prevPhiV = 0.0;
	double Rho = 0.0,Rho_x = 0.0,Rho_y = 0.0;
//-------MOMENTUM-----------------
	double U = 0.0,V = 0.0;
	double Fx = 0.0,Fy = 0.0;
//-----------ENERGY------------------
	double p = 0.0,T = 0.0,qx = 0.0,qy = 0.0;
//-----------------------------------
	double Lambda = 0.0,Mu = 0.0,scG = 0;
//
	double Rho_1k = 0.0,U_1k = 0.0,V_1k = 0.0,T_1k = 0.0;
//
	inline double SqUV(){return U*U+V*V;}
	inline double SqPhixPhiy(){return Phi_x*Phi_x+Phi_y*Phi_y;}
	inline double dPhiU(){return Phi*U-prevPhiU;}
	inline double dPhiV(){return Phi*V-prevPhiV;}
	inline double RhoXUYV(){return Rho_x*U+Rho_y*V;}
	void calcMu();
	double calcTau();
	//
	inline void calcRho_xRho_y()
	{
		Rho_x = PhaseFieldAC::aPhi*Phi_x;
		Rho_y = PhaseFieldAC::aPhi*Phi_y;
	}
	inline void calcFxFy(double F)
	{
		Fx = F*Phi_x;
		Fy = F*Phi_y;
	}
	inline void calcPsi()
	{
		scG = (p > Rho*RT) - (p < Rho*RT);
		Psi = sqrt(2*(p - Rho*RT)/scG);
	}
};

#endif