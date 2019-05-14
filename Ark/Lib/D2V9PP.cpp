#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

#include <iostream>
using std::cout;

extern const char DmQnName[] = "D2V9";

extern double * const xi_u = new double[Q];

extern double * const xi_v = new double[Q];

extern void Update_PseudoPsi(Cell_2D &cell);

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
	double const &Rho = cell.MsQ().Rho;
	double const &U = cell.MsQ().U;
	double const &V = cell.MsQ().V;

	double mEq[Q] = {0.0};

	mEq[0] = Rho;
	mEq[1] = (-2.0 + 3.0*(U*U + V*V))*Rho;
	mEq[2] = (D2Q9::Alpha - 3.0*(U*U + V*V))*Rho;
	mEq[3] =  Rho*U;
	mEq[4] = -Rho*U;
	mEq[5] =  Rho*V;
	mEq[6] = -Rho*V;
	mEq[7] = Rho*(U*U - V*V);
	mEq[8] = Rho*U*V;

	LoopVS(Q)
	{
		cell.f.Eq[k] = IntegralGH(Q,D2Q9::IM[k],mEq);
	}
}
void Update_DVDF_Eqh(Face_2D &face)
{
	double const &Rho = face.MsQh().Rho;
	double const &U = face.MsQh().U;
	double const &V = face.MsQh().V;

	double mEq[Q] = {0.0};

	mEq[0] = Rho;
	mEq[1] = (-2.0 + 3.0*(U*U + V*V))*Rho;
	mEq[2] = (D2Q9::Alpha - 3.0*(U*U + V*V))*Rho;
	mEq[3] =  Rho*U;
	mEq[4] = -Rho*U;
	mEq[5] =  Rho*V;
	mEq[6] = -Rho*V;
	mEq[7] = Rho*(U*U - V*V);
	mEq[8] = Rho*U*V;

	LoopVS(Q)
	{
		face.f.EqhDt[k] = IntegralGH(Q,D2Q9::IM[k],mEq);
	}
}
void Update_MacroVar(Cell_2D& cell)
{
	cell.MsQ().Rho  = IntegralGH(Q,cell.f.Tilde);
	cell.MsQ().U    = IntegralGH(Q,cell.f.Tilde,xi_u)/cell.MsQ().Rho;
	cell.MsQ().V    = IntegralGH(Q,cell.f.Tilde,xi_v)/cell.MsQ().Rho;	
	//
	Update_PseudoPsi(cell);
}
void Update_MacroVar_h(Face_2D& face)
{
	face.MsQh().Rho  = IntegralGH(Q,face.f.BhDt);
	face.MsQh().U    = IntegralGH(Q,face.f.BhDt,xi_u)/face.MsQh().Rho;
	face.MsQh().V    = IntegralGH(Q,face.f.BhDt,xi_v)/face.MsQh().Rho;
}
void Update_DVDF_Source(Cell_2D &cell)
{
	// double const &Rho = cell.MsQ().Rho;
	double const &U = cell.MsQ().U;
	double const &V = cell.MsQ().V;
	double const &Fx = cell.MsQ().Fx;
	double const &Fy = cell.MsQ().Fy;

	double mSo[Q] = {0.0};

	mSo[0] = 0;
	mSo[1] = 6.0*(Fx*U + Fy*V);
	mSo[2] = -mSo[1];
	mSo[3] = Fx;
	mSo[4] = -Fx;
	mSo[5] = Fy;
	mSo[6] = -Fy;
	mSo[7] = 2.0*(Fx*U - Fy*V);
	mSo[8] = (Fx*V + Fy*U);

	LoopVS(Q)
	{
		cell.f.So[k] = IntegralGH(Q,D2Q9::IM[k],mSo);
	}
}
void IntegralShearStress(){}
void unsteadyTaylorGreen(){}