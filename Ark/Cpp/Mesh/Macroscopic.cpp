#include "Macroscopic.h"
inline namespace{
using MQ = MacroQuantity;
}
//
MacroQuantity operator-(const MacroQuantity &lhs,const MacroQuantity &rhs)
{
	MQ tmp = lhs;
	tmp -= rhs;
//
	return tmp;
}
MacroQuantity operator+(const MacroQuantity &lhs,const MacroQuantity &rhs)
{
	MQ tmp = lhs;
	tmp += rhs;
//
	return tmp;
}
MacroQuantity operator*(double a,const MacroQuantity &rhs)
{
	MQ tmp;
	tmp.Psi = a*rhs.Psi;
	tmp.Phi = a*rhs.Phi;
	tmp.Phi_x = a*rhs.Phi_x;
	tmp.Phi_y = a*rhs.Phi_y;
	tmp.laplacianPhi = a*rhs.laplacianPhi;
	tmp.prevPhiU = a*rhs.prevPhiU;
	tmp.prevPhiV = a*rhs.prevPhiV;
	tmp.Rho = a*rhs.Rho;
	tmp.Rho_x = a*rhs.Rho_x;
	tmp.Rho_y = a*rhs.Rho_y;
	tmp.U = a*rhs.U;
	tmp.V = a*rhs.V;
	tmp.Fx = a*rhs.Fx;
	tmp.Fy = a*rhs.Fy;
	tmp.p = a*rhs.p;
	tmp.T = a*rhs.T;
	tmp.qx = a*rhs.qx;
	tmp.qy = a*rhs.qy;
//
	tmp.Lambda = rhs.Lambda;
	tmp.Mu = rhs.Mu;
	tmp.scG = rhs.scG;
//
	return tmp;
}
MacroQuantity operator*(const MacroQuantity &lhs,double a)
{
	return operator*(a,lhs);
}
MacroQuantity& MacroQuantity::operator-=(const MacroQuantity &rhs)
{
	Psi -= rhs.Psi;
	Phi -= rhs.Phi;
	Phi_x -= rhs.Phi_x;
	Phi_y -= rhs.Phi_y;
	laplacianPhi -= rhs.laplacianPhi;
	prevPhiU -= rhs.prevPhiU;
	prevPhiV -= rhs.prevPhiV;
	Rho -= rhs.Rho;
	Rho_x -= rhs.Rho_x;
	Rho_y -= rhs.Rho_y;
	U -= rhs.U;
	V -= rhs.V;
	Fx -= rhs.Fx;
	Fy -= rhs.Fy;
	p -= rhs.p;
	T -= rhs.T;
//
	Lambda = rhs.Lambda;
	Mu = rhs.Mu;
	scG = rhs.scG;
//
	return *this;
}
MacroQuantity& MacroQuantity::operator+=(const MacroQuantity &rhs)
{
	Psi += rhs.Psi;
	Phi += rhs.Phi;
	Phi_x += rhs.Phi_x;
	Phi_y += rhs.Phi_y;
	laplacianPhi += rhs.laplacianPhi;
	prevPhiU += rhs.prevPhiU;
	prevPhiV += rhs.prevPhiV;
	Rho += rhs.Rho;
	Rho_x += rhs.Rho_x;
	Rho_y += rhs.Rho_y;
	U += rhs.U;
	V += rhs.V;
	Fx += rhs.Fx;
	Fy += rhs.Fy;
	p += rhs.p;
	T += rhs.T;
//
	Lambda = rhs.Lambda;
	Mu = rhs.Mu;
	scG = rhs.scG;
//
	return *this;
}