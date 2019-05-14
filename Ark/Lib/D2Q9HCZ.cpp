#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

using namespace PhaseFieldAC;

#include <iostream>
using std::cout;

extern const char DmQnName[] = "D2Q9HCZ"; 

extern double * const xi_u = new double[Q];

extern double * const xi_v = new double[Q];

double const KForce = 1;

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
		// Mu = PhaseFieldAC::NuV + (Phi-PhaseFieldAC::PhiV)*(PhaseFieldAC::NuL-PhaseFieldAC::NuV);
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
		cell.f.Eq[k] = D2Q9::omega[k]*(cell.MsQ().p + cell.MsQ().Rho*GAMMA*RT);
		//
		#ifdef _ARK_FORCE_FLIP
		cell.f.So[k] = 
		D2Q9::omega[k]*
		(
		(1+GAMMA)*((xi_u[k]-cell.MsQ().U)*(cell.MsQ().Fx) + (xi_v[k]-cell.MsQ().V)*(cell.MsQ().Fy))
		+
		GAMMA*RT*((xi_u[k]-cell.MsQ().U)*(cell.MsQ().Rho_x) + (xi_v[k]-cell.MsQ().V)*(cell.MsQ().Rho_y))
		);
		#endif
		//!momentum end
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
		//
		#ifdef _ARK_MOMENTUM_FLIP
		face.f.EqhDt[k] = D2Q9::omega[k]*(face.MsQh().p + face.MsQh().Rho*RT*GAMMA);
//
		#ifdef _ARK_FORCE_FLIP
		face.f.SohDt[k] = 
		D2Q9::omega[k]*
		(
		(1+GAMMA)*((xi_u[k]-face.MsQh().U)*face.MsQh().Fx+(xi_v[k]-face.MsQh().V)*(face.MsQh().Fy))
		+
		GAMMA*RT*((xi_u[k]-face.MsQh().U)*face.MsQh().Rho_x+(xi_v[k]-face.MsQh().V)*(face.MsQh().Rho_y))
		);
		#endif
		#endif
	}
}

void Update_MacroVar(Cell_2D& cell)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	cell.MsQ().Phi  = IntegralGH(Q,cell.h.Tilde);
	cell.MsQ().Rho  = aPhi*cell.MsQ().Phi + bPhi;
	#endif
	//
//
	#ifdef _ARK_MOMENTUM_FLIP
	//!momentum now is updated by the newest index parameter;
	//!see SourceAC.cpp for more infromation; So is dynamic pressure;

	cell.MsQ().U    = (IntegralGH(Q,cell.f.Tilde,xi_u))/RT;
	cell.MsQ().V    = (IntegralGH(Q,cell.f.Tilde,xi_v))/RT;

	// #ifdef _ARK_FORCE_FLIP
	// cell.MsQ().U += ::hDt * cell.MsQ().Fx;
	// cell.MsQ().V += ::hDt * cell.MsQ().Fy;
	// #endif

	// cell.MsQ().U /= cell.MsQ().Rho;
	// cell.MsQ().V /= cell.MsQ().Rho;

	cell.MsQ().p = IntegralGH(Q,cell.f.Tilde);// + ::hDt*cell.MsQ().RhoXUYV()*RT;

	// cell.MsQ().p = IntegralGH(DV_Qv,cell.f.Tilde[0]) - cell.f.Tilde[0][0];
	// cell.MsQ().p += (hDt*cell.MsQ().RhoXUYV()
	// 		  			+ cell.MsQ().Rho*omega[0]*cell.MsQ().SqUV()*Lambda0);//to be solved +
	// cell.MsQ().p *= Kp;
//
	cell.MsQ().calcMu();
	cell.f.tau = cell.MsQ().calcTau();
	cell.Factor();

	#endif
//
	// cell.MsQ().U    = (IntegralGH(cell.f.Tilde[0],xi_u))/(cell.MsQ().Rho*RT);
	// cell.MsQ().V    = (IntegralGH(cell.f.Tilde[0],xi_v))/(cell.MsQ().Rho*RT);
	// cell.MsQ().p    = IntegralGH(cell.f.Tilde[0]);
	// #ifdef _ARK_FORCE_FLIP
	// cell.MsQ().U += hDt*(cell.MsQ().Fx)/cell.MsQ().Rho;
	// cell.MsQ().V += hDt*(cell.MsQ().Fy)/cell.MsQ().Rho;
	// cell.MsQ().p += hDt*((cell.MsQ().Rho_x)*cell.MsQ().U + (cell.MsQ().Rho_y)*cell.MsQ().V);
	// #endif
	// cell.Mu = dynamicViscosity(cell.MsQ().Phi);
	// cell.Tau = cell.Mu/(cell.MsQ().Rho*RT);
	// cell.Factor();
//----------------------------------------------------------------
	// cell.MsQ().U    = (IntegralGH(cell.f.Tilde[0],xi_u))/cell.MsQ().Rho;
	// cell.MsQ().V    = (IntegralGH(cell.f.Tilde[0],xi_v))/cell.MsQ().Rho;
	// #ifdef _ARK_FORCE_FLIP
	// cell.MsQ().U += hDt*(cell.MsQ().Fx)/cell.MsQ().Rho;
	// cell.MsQ().V += hDt*(cell.MsQ().Fy)/cell.MsQ().Rho;
	// #endif
	// cell.MsQ().U    = (IntegralGH(cell.f.Tilde[0],xi_u))/Rho0;
	// cell.MsQ().V    = (IntegralGH(cell.f.Tilde[0],xi_v))/Rho0;
	// cell.MsQ().p    = 0.5*(cell.MsQ().Rho - Rho0)/Lambda0;
}
void Update_MacroVar_h(Face_2D& face)
{

//	face.Ms_h = 0.5*(*(face.owner->Ms) + *(face.neigh->Ms));
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
	face.MsQh().Rho  = aPhi*face.MsQh().Phi + bPhi;
	#endif
	//
	#ifdef _ARK_MOMENTUM_FLIP

	face.MsQh().U    = IntegralGH(Q,face.f.BhDt,xi_u)/RT;
	face.MsQh().V    = IntegralGH(Q,face.f.BhDt,xi_v)/RT;

	#ifdef _ARK_FORCE_FLIP
	face.MsQh().U += 0.5*hDt*face.MsQh().Fx;
	face.MsQh().V += 0.5*hDt*face.MsQh().Fy;
	#endif

	face.MsQh().U /= face.MsQh().Rho;
	face.MsQh().V /= face.MsQh().Rho;

	// face.MsQh().U = 0.0;
	// face.MsQh().V = 0.0;
	face.MsQh().p = IntegralGH(Q,face.f.BhDt) + 0.5*hDt*RT*face.MsQh().RhoXUYV();

	// face.MsQh().p    = IntegralGH(DV_Qv,face.f.BhDt[0])-face.f.BhDt[0][0];
	// face.MsQh().p   += (
	// 				0.5*hDt*(face.MsQh().RhoXUYV())
	// 			    + omega[0]*face.MsQh().Rho*(face.MsQh().SqUV())*Lambda0//to be solved +
	// 			  );
	// face.MsQh().p   *= Kp;
//
	face.MsQh().calcMu();
	face.f.tauh = face.MsQh().calcTau();
	face.Factor();
	#endif
}
void Update_DVDF_Source(Cell_2D&){}
//------------------------------------LGA------------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double xiDotu = (xi_u[QuIndex]*cell.MsQ().U + xi_v[j]*cell.MsQ().V);
// 		cell.f.So[i][j] = 
// 		omega[j]*
// 		(
// 			((xi_u[QuIndex]-cell.MsQ().U) + xiDotu*xi_u[QuIndex]/RT)*(cell.MsQ().Fx)
// 		+	((xi_v[j]-cell.MsQ().V) + xiDotu*xi_v[j]/RT)*(cell.MsQ().Fy)
// 		);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double xiDotu = (xi_u[QuIndex]*face.MsQh().U + xi_v[j]*face.MsQh().V);
// 		face.f.SohDt[i][j] = 
// 		omega[j]*
// 		(
// 			((xi_u[QuIndex]-face.MsQh().U) + xiDotu*xi_u[QuIndex]/RT)*(face.MsQh().Fx)
// 		+	((xi_v[j]-face.MsQh().V) + xiDotu*xi_v[j]/RT)*(face.MsQh().Fy)
// 		);
// 	}
// }
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		cell.f.So[i][j] = ((cell.MsQ().Fx)*(xi_u[QuIndex]) + (cell.MsQ().Fy)*xi_v[j])/RT*omega[j];
// 	}
// }
// void Update_force_h(Face_2D &face)
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
// 		cell.f.So[i][j] = (cell.MsQ().Fx) * (xi_u[QuIndex]-cell.MsQ().U) + (cell.MsQ().Fy) * (xi_v[j]-cell.MsQ().V);
// 		cell.f.So[i][j] *= cell.f.Eq[i][j]/(cell.MsQ().Rho*RT);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		face.f.SohDt[i][j] = face.MsQh().Fx*(xi_u[QuIndex]-face.MsQh().U) + face.MsQh().Fy*(xi_v[j]-face.MsQh().V);
// 		face.f.SohDt[i][j] *= face.f.EqhDt[i][j]/(face.MsQh().Rho*RT);
// 	}
// }
//---------------------------------Guo + Qing Li------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double Rho = cell.MsQ().Rho, U = cell.MsQ().U, V = cell.MsQ().V, Fx = (cell.MsQ().Fx), Fy = (cell.MsQ().Fy);
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
// 		double Rho = face.MsQh().Rho, U = face.MsQh().U, V = face.MsQh().V, Fx = face.MsQh().Fx, Fy = face.MsQh().Fy;
// 		double ex = xi_u[QuIndex], ey = xi_v[j];
// //
// 		double xiDotF = ex*Fx +ey*Fy;
// 		double UFx = (2*U*Fx + KForce*Fx*Fx/Rho)* (ex*ex-RT);
// 		double UFyVFx = (U*Fy + V*Fx + KForce*Fx*Fy/Rho) * (ex*ey);
// 		double VFy = (2*V*Fy + KForce*Fy*Fy/Rho) * (ey*ey-RT);
// 		face.f.SohDt[i][j] = 
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
// 		double xiDotu = (xi_u[QuIndex]*cell.MsQ().U + xi_v[j]*cell.MsQ().V);
// 		cell.f.So[i][j] = 
// 		omega[j]*
// 		(
// 			((xi_u[QuIndex]-cell.MsQ().U)/RT + xiDotu*xi_u[QuIndex]/(RT*RT))*(cell.MsQ().Fx)
// 		+	((xi_v[j]-cell.MsQ().V)/RT + xiDotu*xi_v[j]/(RT*RT))*(cell.MsQ().Fy)
// 		);
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 		double xiDotu = (xi_u[QuIndex]*face.MsQh().U + xi_v[j]*face.MsQh().V);
// 		face.f.SohDt[i][j] = 
// 		omega[j]*
// 		(
// 			((xi_u[QuIndex]-face.MsQh().U)/RT + xiDotu*xi_u[QuIndex]/(RT*RT))*(face.MsQh().Fx)
// 		+	((xi_v[j]-face.MsQh().V)/RT + xiDotu*xi_v[j]/(RT*RT))*(face.MsQh().Fy)
// 		);
// 	}
// }
//-----------------------------------------------------------------------
// void Update_force(Cell_2D &cell)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	cell.f.So[i][j] = (cell.Fx*(xi_u[QuIndex]-cell.MsQ().U) + cell.Fy*(xi_v[j]-cell.MsQ().V))
// 						*2*Lambda0*cell.f.Eq[i][j]/Rho0;
// 	}
// }
// void Update_force_h(Face_2D &face)
// {
// 	for(int i = 0;i < DV_Qu;++i)
// 	for(int j = 0;j < DV_Qv;++j)
// 	{
// 	face.f.SohDt[i][j] = (face.MsQh().Fx*(xi_u[QuIndex]-face.MsQh().U)+face.MsQh().Fy*(xi_v[j]-face.MsQh().V))
// 						*2*Lambda0*face.f.EqhDt[i][j]/Rho0;
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