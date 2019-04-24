#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

extern const char DmQnName[] = "D2Q9PP"; 

extern double * const xi_u = new double[Q];

extern double * const xi_v = new double[Q];

extern void Update_PseudoPsi(Cell_2D &cell);

double const KForce = 1;

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
	double uu,u1;
	uu = cell.MsQ().U*cell.MsQ().U + cell.MsQ().V*cell.MsQ().V;
	LoopVS(Q)
	{
		u1 = (xi_u[k]*cell.MsQ().U + xi_v[k]*cell.MsQ().V)/RT;
		//cell.f.Eq[i][j] = D2Q9::omega[j]*(cell.MsQ().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		cell.f.Eq[k] = D2Q9::omega[k]*cell.MsQ().Rho*(1.0 + u1 + 0.5*u1*u1 - uu*Lambda0);
	}
}
void Update_DVDF_Eqh(Face_2D &face)
{
	double uu,u1;
	uu = face.MsQh().U*face.MsQh().U + face.MsQh().V*face.MsQh().V;
	LoopVS(Q)
	{
		u1 = (xi_u[k]*face.MsQh().U + xi_v[k]*face.MsQh().V)/RT; 
		//face.f.EqhDt[i][j] = D2Q9::omega[j]*(face.MsQh().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		face.f.EqhDt[k] = D2Q9::omega[k]*face.MsQh().Rho*(1.0 + u1 + 0.5*u1*u1 - uu*Lambda0);
	}
}
void Update_MacroVar(Cell_2D& cell)
{
	cell.MsQ().Rho  = IntegralGH(Q,cell.f.Tilde);
	cell.MsQ().U    = (IntegralGH(Q,cell.f.Tilde,xi_u))/cell.MsQ().Rho;
	cell.MsQ().V    = (IntegralGH(Q,cell.f.Tilde,xi_v))/cell.MsQ().Rho;	
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
	LoopVS(Q)
	{
		cell.f.So[k] =
		(cell.MsQ().Fx*(xi_u[k]-cell.MsQ().U) + cell.MsQ().Fy*(xi_v[k]-cell.MsQ().V))
		*cell.f.Eq[k]/(cell.MsQ().Rho*RT);
	}
}
//------------------------------------LGA------------------------------------
// void Update_DVDF_Source(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		cell.f.So[i][j] = (cell.MsQ().Fx*(xi_u[QuIndex]) + cell.MsQ().Fy*xi_v[k])/RT*D2Q9::omega[j];
// 	}
// }
// void Update_DVDF_Source_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		face.f.SohDt[i][j] = (face.MsQh().Fx*xi_u[QuIndex] + face.MsQh().Fy*xi_v[j])/RT*D2Q9::omega[j];
// 	}
// }
// //------------------------------He-Shan-Doolen---------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		cell.f.So[i][j] = (*cell.Fx) * (xi_u[QuIndex]-cell.MsQ().U) + (*cell.Fy) * (xi_v[j]-cell.MsQ().V);
// 		cell.f.So[i][j] *= cell.f.Eq[i][j]/(cell.MsQ().Rho*RT);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		face.forceh[i][j] = face.Fx_h*(xi_u[QuIndex]-face.MsQh().U) + face.Fy_h*(xi_v[j]-face.MsQh().V);
// 		face.forceh[i][j] *= face.fEqh[i][j]/(face.MsQh().Rho*RT);
// 	}
// }
//---------------------------------Guo + Qing Li------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double Rho = cell.MsQ().Rho, U = cell.MsQ().U, V = cell.MsQ().V, Fx = (*cell.Fx), Fy = (*cell.Fy);
// 		double ex = xi_u[QuIndex], ey = xi_v[j];
// //
// 		double xiDotF = ex*Fx +ey*Fy;
// 		double UFx = (2*U*Fx + KForce*Fx*Fx/Rho)* (ex*ex-RT);
// 		double UFyVFx = (U*Fy + V*Fx + KForce*Fx*Fy/Rho) * (ex*ey);
// 		double VFy = (2*V*Fy + KForce*Fy*Fy/Rho) * (ey*ey-RT);
// 		cell.f.So[i][j] = 
// 		D2Q9::omega[j]*
// 		(
// 			xiDotF/RT + (UFx + UFyVFx + VFy)/(2*RT*RT)
// 		);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double Rho = face.MsQh().Rho, U = face.MsQh().U, V = face.MsQh().V, Fx = face.Fx_h, Fy = face.Fy_h;
// 		double ex = xi_u[QuIndex], ey = xi_v[j];
// //
// 		double xiDotF = ex*Fx +ey*Fy;
// 		double UFx = (2*U*Fx + KForce*Fx*Fx/Rho)* (ex*ex-RT);
// 		double UFyVFx = (U*Fy + V*Fx + KForce*Fx*Fy/Rho) * (ex*ey);
// 		double VFy = (2*V*Fy + KForce*Fy*Fy/Rho) * (ey*ey-RT);
// 		face.forceh[i][j] = 
// 		D2Q9::omega[j]*
// 		(
// 			xiDotF/RT + (UFx + UFyVFx + VFy)/(2*RT*RT)
// 		);
// 	}
// }
//----------------------------------Lishi-Luo force----------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double xiDotu = (xi_u[j]*cell.MsQ().U + xi_v[j]*cell.MsQ().V);
// 		cell.f.So[i][j] = 
// 		D2Q9::omega[j]*
// 		(
// 			((xi_u[j]-cell.MsQ().U)/RT + xiDotu*xi_u[j]/(RT*RT))*(*cell.Fx)
// 		+	((xi_v[j]-cell.MsQ().V)/RT + xiDotu*xi_v[j]/(RT*RT))*(*cell.Fy)
// 		);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double xiDotu = (xi_u[j]*face.MsQh().U + xi_v[j]*face.MsQh().V);
// 		face.forceh[i][j] = 
// 		D2Q9::omega[j]*
// 		(
// 			((xi_u[j]-face.MsQh().U)/RT + xiDotu*xi_u[j]/(RT*RT))*(face.Fx_h)
// 		+	((xi_v[j]-face.MsQh().V)/RT + xiDotu*xi_v[j]/(RT*RT))*(face.Fy_h)
// 		);
// 	}
// }
//-----------------------------------------------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	cell.f.So[i][j] = (cell.Fx*(xi_u[j]-cell.MsQ().U) + cell.Fy*(xi_v[j]-cell.MsQ().V))
// 						*2*Lambda0*cell.f.Eq[i][j]/Rho0;
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	face.forceh[i][j] = (face.Fx_h*(xi_u[j]-face.MsQh().U)+face.Fy_h*(xi_v[j]-face.MsQh().V))
// 						*2*Lambda0*face.fEqh[i][j]/Rho0;
// 	}
// }
// void Wall_3_Boundary(Face_2D &face)
// {

// }
void IntegralShearStress(){}
//
//-----------------------------Inc unsteady Taylor Green vortex------------------
extern void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p);
void unsteadyTaylorGreen()
{	
	for(int n = 0;n != Cells;++n)
	{
		double x = CellArray[n].xc, y = CellArray[n].yc;
		TaylorGreenVortex(0.0,x,y,CellArray[n].MsQ().U,CellArray[n].MsQ().V,CellArray[n].MsQ().p);
		CellArray[n].MsQ().Rho = Rho0 + CellArray[n].MsQ().p/RT;
		CellArray[n].MsQ().T = T0;
		CellArray[n].MsQ().Lambda = Lambda0;
		CellArray[n].MsQ().Mu = Mu0;
		CellArray[n].f.tau = 2.0*Nu0*Lambda0;
		CellArray[n].Factor();
		Update_DVDF_Eq(CellArray[n]);
//---------------------------------Initialize_fT-------------------------------
		double grad_ux = U0*sin(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vx = U0*cos(2.0*PI*x)*cos(2.0*PI*y);
//
		double grad_uy = -grad_vx;
		double grad_vy = -grad_ux;

		double grad_ut = Nu0*4.0*PI*U0*cos(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vt = -Nu0*4.0*PI*U0*sin(2.0*PI*x)*cos(2.0*PI*y);
		LoopVS(Q)
		{
			double u1 = (xi_u[k]) * CellArray[n].MsQ().U + (xi_v[k]) * CellArray[n].MsQ().V;
			double A_u = ((xi_u[k]) + u1*(xi_u[k])/RT - CellArray[n].MsQ().U);
			double A_v = ((xi_v[k]) + u1*(xi_v[k])/RT - CellArray[n].MsQ().V);
			double A = D2Q9::omega[k]*Rho0/RT;
			double fEq_t = A * (A_u*grad_ut + A_v*grad_vt);
			double fEq_x = A * (A_u*grad_ux + A_v*grad_vx) * (xi_u[k]);
			double fEq_y = A * (A_u*grad_uy + A_v*grad_vy) * (xi_v[k]);
			double f = CellArray[n].f.Eq[k] - CellArray[n].f.tau*(fEq_t + fEq_x + fEq_y);
			CellArray[n].f.Tilde[k] = f; //- 0.5*dt*(fEq_t + fEq_x + fEq_y);
		}
	}
	for(int n = 0;n < Faces;++n)
	{
		FaceArray[n].MsQh().T = T0;
		FaceArray[n].MsQh().Mu = Mu0;
		FaceArray[n].MsQh().Lambda = Lambda0;
		FaceArray[n].f.tauh = 2.0*Nu0*Lambda0;
		FaceArray[n].Factor();
	}
}