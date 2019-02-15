#include "DUGKSDeclaration.h"
#include "GaussHermite.h"

using namespace PhaseFieldAC;

#include <iostream>
using std::cout;

extern const char DmQnName[] = "D2Q9HCZ"; 

extern double * const xi_u = new double[DV_Qv];

extern double * const xi_v = new double[DV_Qv];

const double omega[DV_Qv]={4.0/9.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0,
						1.0/9.0, 1.0/36.0};

double const KForce = 1;

//double const Kp = RT*9/5.0;
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
	for(int j = 0;j < DV_Qv;++j)
	{
		xi_u[QuIndex] = MaSpan*D2Q9::xi_u[QuIndex];
		xi_v[j] = MaSpan*D2Q9::xi_v[j];
	}
}
void Update_phi_Eq(Cell_2D &cell)
{
	double uu,u1,GAMMA;
	uu = cell.MsQ().SqUV();
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		u1 = (xi_u[QuIndex]*cell.MsQ().U + xi_v[j]*cell.MsQ().V)/RT;
		GAMMA = u1 + 0.5*u1*u1 - uu*Lambda0;
		//cell.f.Eq[i][j] = omega[j]*(cell.MsQ().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		//cell.f.Eq[i][j] = omega[j]*cell.MsQ().Rho*(1.0 + u1 + 0.5*u1*u1 - uu*Lambda0);
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.Eq[i][j] = omega[j]*(1.0 + u1) * (cell.MsQ().Phi);
		#endif
		//
		//cell.f.Eq[i][j] = omega[j]*(cell.MsQ().p + cell.MsQ().Rho*RT*GAMMA);
		cell.f.Eq[i][j] = omega[j]*(cell.MsQ().p + cell.MsQ().Rho*GAMMA*RT);
		//
		#ifdef _ARK_FORCE_FLIP
		cell.f.So[i][j] = 
		omega[j]*
		(
		(1+GAMMA)*((xi_u[j]-cell.MsQ().U)*(cell.MsQ().Fx) + (xi_v[j]-cell.MsQ().V)*(cell.MsQ().Fy))
		+
		GAMMA*RT*((xi_u[j]-cell.MsQ().U)*(cell.MsQ().Rho_x) + (xi_v[j]-cell.MsQ().V)*(cell.MsQ().Rho_y))
		);
		// cell.f.So[i][j] = 
		// omega[j]*
		// (
		// 	// (
		// 	// 	(xi_u[QuIndex]-cell.MsQ().U)*(cell.MsQ().Fx)
		// 	// +   (xi_v[j]-cell.MsQ().V)*(cell.MsQ().Fy)
		// 	// )*(1+GAMMA)/RT
		// 	(xi_u[QuIndex]*(cell.MsQ().Fx) + xi_v[j]*cell.MsQ().Fy)/RT
		// +	u1*(xi_u[QuIndex]*(cell.MsQ().Rho_x) + xi_v[j]*(cell.MsQ().Rho_y))
		// );
		#endif
	}
}
void Update_phi_Eqh(Face_2D &face)
{
	double uu,u1,GAMMA;
	uu = face.MsQh().SqUV();
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		u1 = (xi_u[QuIndex]*face.MsQh().U + xi_v[j]*face.MsQh().V)/RT; 
		GAMMA = u1 + 0.5*u1*u1 - uu*Lambda0;
		//face.f.EqhDt[i][j] = omega[j]*(face.MsQh().Rho + Rho0*(u1 + 0.5*u1*u1 - uu*Lambda0));
		//face.f.EqhDt[i][j] = omega[j]*face.MsQh().Rho*(1.0 + u1 + 0.5*u1*u1 - uu*Lambda0);
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.EqhDt[i][j] = omega[j]*face.MsQh().Phi*(1.0 + u1);
		#endif
		//
		face.f.EqhDt[i][j] = omega[j]*(face.MsQh().p + face.MsQh().Rho*RT*GAMMA);
		//face.f.EqhDt[i][j] = omega[j]*(face.MsQh().p/RT + face.MsQh().Rho*GAMMA)
//
		#ifdef _ARK_FORCE_FLIP
		face.f.SohDt[i][j] = 
		omega[j]*
		(
		(1+GAMMA)*((xi_u[QuIndex]-face.MsQh().U)*face.MsQh().Fx+(xi_v[j]-face.MsQh().V)*(face.MsQh().Fy))
		+
		GAMMA*RT*((xi_u[QuIndex]-face.MsQh().U)*face.MsQh().Rho_x+(xi_v[j]-face.MsQh().V)*(face.MsQh().Rho_y))
		);
		// face.f.SohDt[i][j] = 
		// omega[j]*
		// (
		// 	// (
		// 	//     (xi_u[QuIndex]-face.MsQh().U)*face.MsQh().Fx
		// 	//  +  (xi_v[j]-face.MsQh().V)*face.MsQh().Fy
		// 	// )*(1+GAMMA)/RT
		// 	(xi_u[QuIndex]*face.MsQh().Fx + xi_v[j]*face.MsQh().Fy)/RT
		// +	u1*(xi_u[QuIndex]*face.MsQh().Rho_x + xi_v[j]*face.MsQh().Rho_y)
		// );
		#endif
	}
}

void Update_MacroVar(Cell_2D& cell)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	cell.MsQ().Phi  = IntegralGH(DV_Qv,cell.h.Tilde[0]);
	#endif
	//
	cell.MsQ().Rho  = aPhi*cell.MsQ().Phi + bPhi;
//
	cell.MsQ().U    = (IntegralGH(DV_Qv,cell.f.Tilde[0],xi_u))/(cell.MsQ().Rho*RT);
	cell.MsQ().V    = (IntegralGH(DV_Qv,cell.f.Tilde[0],xi_v))/(cell.MsQ().Rho*RT);
	#ifdef _ARK_FORCE_FLIP
	cell.MsQ().U += hDt*(cell.MsQ().Fx)/(cell.MsQ().Rho);
	cell.MsQ().V += hDt*(cell.MsQ().Fy)/(cell.MsQ().Rho);
	#endif
//
	// cell.MsQ().U = 0.0;
	// cell.MsQ().V = 0.0;
	cell.MsQ().p = IntegralGH(DV_Qv,cell.f.Tilde[0]) + ::hDt*cell.MsQ().RhoXUYV()*RT;

	// cell.MsQ().p = IntegralGH(DV_Qv,cell.f.Tilde[0]) - cell.f.Tilde[0][0];
	// cell.MsQ().p += (hDt*cell.MsQ().RhoXUYV()
	// 		  			+ cell.MsQ().Rho*omega[0]*cell.MsQ().SqUV()*Lambda0);//to be solved +
	// cell.MsQ().p *= Kp;
//
	cell.MsQ().calcMu();
	cell.f.tau = cell.MsQ().calcTau();
	cell.Factor();
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
//MsQh().
	face.MsQh().Rho_x = 0.5*((face.owner->msq->Rho_x) + (face.neigh->msq->Rho_x));
	face.MsQh().Rho_y = 0.5*((face.owner->msq->Rho_y) + (face.neigh->msq->Rho_y));

	face.MsQh().Fx = 0.5*((face.owner->msq->Fx) + (face.neigh->msq->Fx));
	face.MsQh().Fy = 0.5*((face.owner->msq->Fy) + (face.neigh->msq->Fy));

	#ifdef _ARK_ALLENCAHN_FLIP
	face.MsQh().Phi  = IntegralGH(DV_Qv,face.h.BhDt[0]);
	#endif
	//
	face.MsQh().Rho  = aPhi*face.MsQh().Phi + bPhi;
//
	face.MsQh().U    = IntegralGH(DV_Qv,face.f.BhDt[0],xi_u)/(face.MsQh().Rho*RT);
	face.MsQh().V    = IntegralGH(DV_Qv,face.f.BhDt[0],xi_v)/(face.MsQh().Rho*RT);
	#ifdef _ARK_FORCE_FLIP
	face.MsQh().U += 0.5*hDt*face.MsQh().Fx/(face.MsQh().Rho);
	face.MsQh().V += 0.5*hDt*face.MsQh().Fy/(face.MsQh().Rho);
	#endif

	// face.MsQh().U = 0.0;
	// face.MsQh().V = 0.0;
	face.MsQh().p = IntegralGH(DV_Qv,face.f.BhDt[0]) + 0.5*hDt*RT*face.MsQh().RhoXUYV();

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
//------------------------------------LGA------------------------------------
void Update_DVDF_Source(Cell_2D &cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.So[i][j] = ((cell.MsQ().Phi_x)*(xi_u[QuIndex]) + (cell.MsQ().Phi_y)*xi_v[j])
						  *omega[j];
		#endif
		//cell.f.So[i][j] = cell.h.So[i][j];
		// ((cell.MsQ().Phi_x)*(xi_u[QuIndex]) + (cell.MsQ().Phi_y)*xi_v[j])
		// 				* (*cell.Ms)*omega[j];
	}
}
void Update_DVDF_Source_h(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.SohDt[i][j] = (face.MsQh().Phi_x*xi_u[QuIndex] + face.MsQh().Phi_y*xi_v[j])
						    *omega[j];
		#endif
		//face.f.SohDt[i][j] = face.h.SohDt[i][j];
		// (face.MsQh().Phi_x*xi_u[QuIndex] + face.MsQh().Phi_y*xi_v[j])
		// 				   * (face.Ms_h)*omega[j];
	}
}
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
		CellArray[n].f.tau = Tau0;
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
			double u1 = (xi_u[QuIndex]) * CellArray[n].MsQ().U + (xi_v[j]) * CellArray[n].MsQ().V;
			double A_u = ((xi_u[QuIndex]) + u1*(xi_u[QuIndex])/RT - CellArray[n].MsQ().U);
			double A_v = ((xi_v[j]) + u1*(xi_v[j])/RT - CellArray[n].MsQ().V);
			double A = omega[j]*Rho0/RT;
			double fEq_t = A * (A_u*grad_ut + A_v*grad_vt);
			double fEq_x = A * (A_u*grad_ux + A_v*grad_vx) * (xi_u[QuIndex]);
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
		FaceArray[n].f.tauh = Tau0;
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
			face.f.hDt[k] = face.ah*face.f.BhDt[k] + face.bh*face.f.EqhDt[k];
			out += face.f.hDt[k]*xi_dot_n;
		}
		else if(xi_dot_n < 0)
		{
			face.f.hDt[k] = face.f.EqhDt[k];
			in += face.f.hDt[k]*xi_dot_n;
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
		face.f.hDt[k] =  face.f.EqhDt[k] + cell.aNEq*(cell.fT[k] - cell.fEq[k]);
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