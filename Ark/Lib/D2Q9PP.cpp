#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

extern const char DmQnName[] = "D2Q9PP"; 

extern double * const xi_u = new double[DV_Qv];

extern double * const xi_v = new double[DV_Qv];

const double omega[DV_Qv]={4.0/9.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0};

double const KForce = 1;

extern void Update_PseudoPsi(Cell_2D &cell);

extern void Update_PseudoPsi(Face_2D &face);

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
	for(int j = 0;j < DV_Qv;++j)
	{
		xi_u[j] = MaSpan*D2Q9::xi_u[j];
		xi_v[j] = MaSpan*D2Q9::xi_v[j];
	}
}
void Update_phi_Eq(Cell_2D &cell)
{
	double uu,u1;
	uu = cell.MsQ().U*cell.MsQ().U + cell.MsQ().V*cell.MsQ().V;
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		u1 = (xi_u[j]*cell.MsQ().U + xi_v[j]*cell.MsQ().V)/RT;
		//cell.f.Eq[i][j] = omega[j]*(cell.MsQ().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		cell.f.Eq[i][j] = omega[j]*cell.MsQ().Rho*(1.0 + u1 + 0.5*u1*u1 - uu*Lambda0);
	}
}
void Update_phi_Eqh(Face_2D &face)
{
	double uu,u1;
	uu = face.MsQh().U*face.MsQh().U + face.MsQh().V*face.MsQh().V;
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		u1 = (xi_u[j]*face.MsQh().U + xi_v[j]*face.MsQh().V)/RT; 
		//face.f.EqhDt[i][j] = omega[j]*(face.MsQh().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		face.f.EqhDt[i][j] = omega[j]*face.MsQh().Rho*(1.0 + u1 + 0.5*u1*u1 - uu*Lambda0);
	}
}
// double vectorDot(const double *const first,const double *const second)
// {
// 	double sum = 0;
// 	for(int j = 0;j < DV_Qv;++j)
// 		sum += first[j]*second[j];
// 	return sum;
// }
// double vectorDot(const double *const first)
// {
// 	double sum = 0;
// 	for(int j = 0;j < DV_Qv;++j)
// 		sum += first[j];
// 	return sum;
// }
void Update_MacroVar(Cell_2D& cell)
{
	cell.MsQ().Rho  = IntegralGH(DV_Qv,cell.f.Tilde[0]);
	cell.MsQ().U    = (IntegralGH(DV_Qv,cell.f.Tilde[0],xi_u))/cell.MsQ().Rho;
	cell.MsQ().V    = (IntegralGH(DV_Qv,cell.f.Tilde[0],xi_v))/cell.MsQ().Rho;	
	#ifdef _ARK_FORCE_FLIP
	cell.MsQ().U += hDt*cell.MsQ().Fx/cell.MsQ().Rho;
	cell.MsQ().V += hDt*cell.MsQ().Fy/cell.MsQ().Rho;
	#endif

	#ifdef _ARK_PSEUDOPSI_FLIP
	Update_PseudoPsi(cell);
	#endif
	// cell.MsQ().U    = (IntegralGH(DV_Qv,cell.f.Tilde[0],xi_u))/Rho0;
	// cell.MsQ().V    = (IntegralGH(DV_Qv,cell.f.Tilde[0],xi_v))/Rho0;	
	// #ifdef _ARK_FORCE_FLIP
	// cell.MsQ().U += hDt*cell.MsQ().Fx/Rho0;
	// cell.MsQ().V += hDt*cell.MsQ().Fy/Rho0;
	// #endif
	//cell.p    = 0.5*(cell.MsQ().Rho - Rho0)/Lambda0;
}
void Update_MacroVar_h(Face_2D& face)
{
	face.MsQh().Rho  = IntegralGH(DV_Qv,face.f.BhDt[0]);
	face.MsQh().U    = IntegralGH(DV_Qv,face.f.BhDt[0],xi_u)/face.MsQh().Rho;
	face.MsQh().V    = IntegralGH(DV_Qv,face.f.BhDt[0],xi_v)/face.MsQh().Rho;
	//
	#ifdef _ARK_FORCE_FLIP
	// face.MsQh().Fx = 0.5*((face.owner->msq->Fx) + (face.neigh->msq->Fx));
	// face.MsQh().Fy = 0.5*((face.owner->msq->Fy) + (face.neigh->msq->Fy));
	face.MsQh().U += 0.5*hDt*face.MsQh().Fx/face.MsQh().Rho;
	face.MsQh().V += 0.5*hDt*face.MsQh().Fy/face.MsQh().Rho;
	#endif

	#if defined _ARK_PSEUDOPSI_FLIP && defined _ARK_FORCE_FLIP
	Update_PseudoPsi(face);
	#endif
	// face.MsQh().U    = IntegralGH(DV_Qv,face.f.BhDt[0],xi_u)/Rho0;
	// face.MsQh().V    = IntegralGH(DV_Qv,face.f.BhDt[0],xi_v)/Rho0;
	// #ifdef _ARK_FORCE_FLIP
	// face.MsQh().Fx = 0.5*((face.owner->msq->Fx) + (face.neigh->msq->Fx));
	// face.MsQh().Fy = 0.5*((face.owner->msq->Fy) + (face.neigh->msq->Fy));
	// //
	// face.MsQh().U += 0.5*hDt*face.MsQh().Fx/Rho0;
	// face.MsQh().V += 0.5*hDt*face.MsQh().Fy/Rho0;
	// #endif
	// face.p_h    = 0.5*(face.MsQh().Rho - Rho0)/Lambda0;
}
void Update_DVDF_Source(Cell_2D &cell)
{
	LoopVS(DV_Qu,DV_Qv)
	{
		cell.f.So[i][j] =
		(cell.MsQ().Fx*(xi_u[j]-cell.MsQ().U)+cell.MsQ().Fy*(xi_v[j]-cell.MsQ().V))
		*cell.f.Eq[i][j]/(cell.MsQ().Rho*RT);
	}
}
void Update_DVDF_Source_h(Face_2D &face)
{
	LoopVS(DV_Qu,DV_Qv)
	{
		face.f.SohDt[i][j] =
		(face.MsQh().Fx*(xi_u[j]-face.MsQh().U) + face.MsQh().Fy*(xi_v[j]-face.MsQh().V))
		*face.f.EqhDt[i][j]/(face.MsQh().Rho*RT);
	}
}
//------------------------------------LGA------------------------------------
// void Update_DVDF_Source(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		cell.f.So[i][j] = (cell.MsQ().Fx*(xi_u[QuIndex]) + cell.MsQ().Fy*xi_v[j])/RT*omega[j];
// 	}
// }
// void Update_DVDF_Source_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		face.f.SohDt[i][j] = (face.MsQh().Fx*xi_u[QuIndex] + face.MsQh().Fy*xi_v[j])/RT*omega[j];
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
// 		omega[j]*
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
// 		omega[j]*
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
// 		omega[j]*
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
// 		omega[j]*
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
		Update_phi_Eq(CellArray[n]);
//---------------------------------Initialize_fT-------------------------------
		double grad_ux = U0*sin(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vx = U0*cos(2.0*PI*x)*cos(2.0*PI*y);
//
		double grad_uy = -grad_vx;
		double grad_vy = -grad_ux;

		double grad_ut = Nu0*4.0*PI*U0*cos(2.0*PI*x)*sin(2.0*PI*y);
		double grad_vt = -Nu0*4.0*PI*U0*sin(2.0*PI*x)*cos(2.0*PI*y);
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			double u1 = (xi_u[j]) * CellArray[n].MsQ().U + (xi_v[j]) * CellArray[n].MsQ().V;
			double A_u = ((xi_u[j]) + u1*(xi_u[j])/RT - CellArray[n].MsQ().U);
			double A_v = ((xi_v[j]) + u1*(xi_v[j])/RT - CellArray[n].MsQ().V);
			double A = omega[j]*Rho0/RT;
			double fEq_t = A * (A_u*grad_ut + A_v*grad_vt);
			double fEq_x = A * (A_u*grad_ux + A_v*grad_vx) * (xi_u[j]);
			double fEq_y = A * (A_u*grad_uy + A_v*grad_vy) * (xi_v[j]);
			double f = CellArray[n].f.Eq[i][j] - CellArray[n].f.tau*(fEq_t + fEq_x + fEq_y);
			CellArray[n].f.Tilde[i][j] = f; //- 0.5*dt*(fEq_t + fEq_x + fEq_y);
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