#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

extern double * const xi_u = new double[Q];

extern double * const xi_v = new double[Q];

extern const char DmQnName[] = "D2Q9"; 

void MacroQuantity::calcMu()
{
		// Mu=PhaseFieldAC::MuV + (Phi-PhaseFieldAC::PhiV)*(PhaseFieldAC::MuL-PhaseFieldAC::MuV);
		
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
		// Mu = PhaseFieldAC::NuV + (Phi-PhaseFieldAC::PhiV)*(PhaseFieldAC::NuL-PhaseFieldAC::NuV);
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
	double const &Rho = cell.MsQ().Rho, &U = cell.MsQ().U, &V = cell.MsQ().V;
	double const &Fx = cell.MsQ().Fx, &Fy = cell.MsQ().Fy;
	//
	double uu,u1;
	uu = U*U + V*V;
	LoopVS(Q)
	{
		u1 = (xi_u[k]*U + xi_v[k]*V)/RT;
		//cell.f.Eq[i][j] = D2Q9::omega[k]*(cell.MsQ().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.Eq[k] = D2Q9::omega[k]*(Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		cell.f.So[k] = ((xi_u[k] - U)*Fx + (xi_v[k] - V)*Fy)*cell.f.Eq[k]/(Rho0*RT);
		#endif
	}
}
void Update_DVDF_Eqh(Face_2D &face)
{
	double const &Rhoh = face.MsQh().Rho, &Uh = face.MsQh().U, &Vh = face.MsQh().V;
	double const &Fxh = face.MsQh().Fx, &Fyh = face.MsQh().Fy;
	//
	double uu,u1;
	uu = Uh*Uh + Vh*Vh;
	LoopVS(Q)
	{
		u1 = (xi_u[k]*Uh + xi_v[k]*Vh)/RT; 
		//face.f.EqhDt[i][j] = D2Q9::omega[k]*(face.MsQh().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		#ifdef _ARK_MOMENTUM_FLIP
		face.f.EqhDt[k] = D2Q9::omega[k]*(Rhoh + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		face.f.SohDt[k] = ((xi_u[k] - Uh)*Fxh + (xi_v[k] - Vh)*Fyh)*face.f.EqhDt[k]/(Rho0*RT);
		#endif
	}
}
void Update_MacroVar(Cell_2D& cell)
{
	cell.MsQ().Rho  = IntegralGH(Q,cell.f.Tilde);
	//
	cell.MsQ().U    = IntegralGH(Q,cell.f.Tilde,xi_u) + 0.5*::dt*cell.MsQ().Fx;
	cell.MsQ().U    /= Rho0;
	//
	cell.MsQ().V    = IntegralGH(Q,cell.f.Tilde,xi_v) + 0.5*::dt*cell.MsQ().Fy;
	cell.MsQ().V    /= Rho0;
	cell.MsQ().p = (cell.MsQ().Rho - Rho0)*RT;
}
void Update_MacroVar_h(Face_2D& face)
{
	face.MsQh().Rho  = IntegralGH(Q,face.f.BhDt);
	//
	face.MsQh().U    = IntegralGH(Q,face.f.BhDt,xi_u) + 0.5*::hDt*face.MsQh().Fx;
	face.MsQh().U    /= Rho0;
	//
	face.MsQh().V    = IntegralGH(Q,face.f.BhDt,xi_v) + 0.5*::hDt*face.MsQh().Fy;
	face.MsQh().V    /= Rho0;
	face.MsQh().p = (face.MsQh().Rho - Rho0)*RT;
}
//------------------------------------LGA------------------------------------
// void Update_DVDF_Source(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		cell.f.So[i][j] = (cell.MsQ().Fx*(xi_u[QuIndex]) + cell.MsQ().Fy*xi_v[k])/RT*omega[k];
// 	}
// }
// void Update_DVDF_Source_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		face.f.SohDt[i][j] = (face.MsQh().Fx*xi_u[QuIndex] + face.MsQh().Fy*xi_v[k])/RT*omega[k];
// 	}
// }
// //------------------------------He-Shan-Doolen---------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		cell.f.So[i][j] = (*cell.Fx) * (xi_u[QuIndex]-cell.MsQ().U) + (*cell.Fy) * (xi_v[k]-cell.MsQ().V);
// 		cell.f.So[i][j] *= cell.f.Eq[i][j]/(cell.MsQ().Rho*RT);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		face.forceh[i][j] = face.Fx_h*(xi_u[QuIndex]-face.MsQh().U) + face.Fy_h*(xi_v[k]-face.MsQh().V);
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
// 		double ex = xi_u[QuIndex], ey = xi_v[k];
// //
// 		double xiDotF = ex*Fx +ey*Fy;
// 		double UFx = (2*U*Fx + KForce*Fx*Fx/Rho)* (ex*ex-RT);
// 		double UFyVFx = (U*Fy + V*Fx + KForce*Fx*Fy/Rho) * (ex*ey);
// 		double VFy = (2*V*Fy + KForce*Fy*Fy/Rho) * (ey*ey-RT);
// 		cell.f.So[i][j] = 
// 		omega[k]*
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
// 		double ex = xi_u[QuIndex], ey = xi_v[k];
// //
// 		double xiDotF = ex*Fx +ey*Fy;
// 		double UFx = (2*U*Fx + KForce*Fx*Fx/Rho)* (ex*ex-RT);
// 		double UFyVFx = (U*Fy + V*Fx + KForce*Fx*Fy/Rho) * (ex*ey);
// 		double VFy = (2*V*Fy + KForce*Fy*Fy/Rho) * (ey*ey-RT);
// 		face.forceh[i][j] = 
// 		omega[k]*
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
// 		double xiDotu = (xi_u[k]*cell.MsQ().U + xi_v[k]*cell.MsQ().V);
// 		cell.f.So[i][j] = 
// 		omega[k]*
// 		(
// 			((xi_u[k]-cell.MsQ().U)/RT + xiDotu*xi_u[k]/(RT*RT))*(*cell.Fx)
// 		+	((xi_v[k]-cell.MsQ().V)/RT + xiDotu*xi_v[k]/(RT*RT))*(*cell.Fy)
// 		);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double xiDotu = (xi_u[k]*face.MsQh().U + xi_v[k]*face.MsQh().V);
// 		face.forceh[i][j] = 
// 		omega[k]*
// 		(
// 			((xi_u[k]-face.MsQh().U)/RT + xiDotu*xi_u[k]/(RT*RT))*(face.Fx_h)
// 		+	((xi_v[k]-face.MsQh().V)/RT + xiDotu*xi_v[k]/(RT*RT))*(face.Fy_h)
// 		);
// 	}
// }
//-----------------------------------------------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	cell.f.So[i][j] = (cell.Fx*(xi_u[k]-cell.MsQ().U) + cell.Fy*(xi_v[k]-cell.MsQ().V))
// 						*2*Lambda0*cell.f.Eq[i][j]/Rho0;
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	face.forceh[i][j] = (face.Fx_h*(xi_u[k]-face.MsQh().U)+face.Fy_h*(xi_v[k]-face.MsQh().V))
// 						*2*Lambda0*face.fEqh[i][j]/Rho0;
// 	}
// }
// void Wall_3_Boundary(Face_2D &face)
// {

// }
void IntegralShearStress(){}
void Update_DVDF_Source(Cell_2D &cell){}
//
//-----------------------------Inc unsteady Taylor Green vortex------------------
extern void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p);
void unsteadyTaylorGreen()
{	
	LoopPS(Cells)
	{
		double x = CellArray[n].xc, y = CellArray[n].yc;
		TaylorGreenVortex(0.0,x,y,CellArray[n].MsQ().U,CellArray[n].MsQ().V,CellArray[n].MsQ().p);
		CellArray[n].MsQ().Rho = Rho0 + CellArray[n].MsQ().p/RT;
		CellArray[n].MsQ().T = T0;
		CellArray[n].MsQ().Lambda = Lambda0;
		CellArray[n].MsQ().Mu = Mu0;
		CellArray[n].f.tau = 2.0*Nu0*Lambda0;
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
	LoopPS(Faces)
	{
		FaceArray[n].MsQh().T = T0;
		FaceArray[n].MsQh().Mu = Mu0;
		FaceArray[n].MsQh().Lambda = Lambda0;
		FaceArray[n].f.tauh = 2.0*Nu0*Lambda0;
	}
}
/*void DiffusiveScatter(Face_2D &face)
{
	double uu = face.uh*face.uh + face.vh*face.vh;
	double u1,xi_dot_n;
	double rhohWall_A = 0.0, rhohWall_b = 0.0;
	for(int k = 0;k < Q;++k)
	{
		xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		u1 = xi[k].u*face.uh + xi[k].v*face.vh;
		if(xi_dot_n > 0)
		{
			Solve_Interior_fBh(face,face.owner,k);
			rhohWall_b += xi_dot_n*(face.ah*face.f.BhDt[k] + face.bh*ConstanceInfEq(u1,uu,k));
			rhohWall_A += xi_dot_n*bh*omega[k];
		}
		else if(xi_dot_n < 0)
		{
			rhohWall_b += xi_dot_n*ConstanceInfEq(u1,uu,k);
			rhohWall_A += xi_dot_n*omega[k];
		}
	}
	face.rhoh = -rhohWall_b/rhohWall_A;
	Update_fEqh(face);
	double out = 0.0,in = 0.0;
	for(int k = 0;k < Q;++k)
	{
		xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		if(xi_dot_n > 0)
		{
			face.fh[k] = face.ah*face.f.BhDt[k] + face.bh*face.fEqh[k];
			out += face.fh[k]*xi_dot_n;
		}
		else if(xi_dot_n < 0)
		{
			face.fh[k] = face.fEqh[k];
			in += face.fh[k]*xi_dot_n;
		}
	}
	if(in + out > 1.0E-16 || in + out < -1.0E-16)
	{
		cout <<in + out<<endl;
		getchar();
	}
}

void NonEquilibriumExtrapolation(Face_2D &face)
{
	Cell_2D &cell = *face.owner;
//
	face.rhoh = cell.rho;
//
	Update_fEqh(face);
	for(int k = 0;k != Q;++k)
	{
		face.fh[k] =  face.fEqh[k] + cell.aNEq*(cell.f.Tilde[k] - cell.fEq[k]);
	}
}
void BounceBack(Face_2D &face)
{
	face.rhoh = face.owner->rho;
	Update_fEqh(face);
	for(int k = 0;k < Q;++k)
	{
		double xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		if(xi_dot_n > 0.0)
		{
			Solve_Interior_fBh(face,face.owner,k);
		}
	}
	for(int k = 0;k < Q;++k)
	{
		double xi_dot_n = face.Vx * xi[k].u + face.Vy * xi[k].v;
		if(xi_dot_n < 0.0)
		{
			double xi_dot_u = face.uh * xi[k].u + face.vh * xi[k].v;
			if(k < 5)
			{	
				face.f.BhDt[k] = face.f.BhDt[k + 4] + 2.0*face.rhoh*omega[k]*xi_dot_u/RT;
			}
			else
			{
				face.f.BhDt[k] = face.f.BhDt[k - 4] + 2.0*face.rhoh*omega[k]*xi_dot_u/RT;
			}

		}
	}
	Update_fh(face);
}*/