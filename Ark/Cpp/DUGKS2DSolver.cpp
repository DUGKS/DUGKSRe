#include <iostream>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include "DUGKSDeclaration.h"
using std::ios;
using std::setiosflags;
using std::setprecision;
using std::cout;
using std::endl;

double const aTP = 4.0/3.0, bTP = 1.0/3.0;

int const ThreadNum = 1;//omp_get_max_threads();

double ResidualPer1k = 1.0;

double const PRECISION = 1.0E-16;

double const XDEBUG = 13.25,YDEBUG = 0.25;


//--------------------------DEBUG-------------------------------

extern void Output_L2Norm(double const &t,double &L2_uv, double &L2_p);

extern void Output_Flowfield(double const &t,int step);

extern void Output_SumRho(double t);

extern void Output_SumEk(double t);

extern void Output_RTIlocation(double t,int const X);

extern void Output_Residual(double t,double Residual);

extern void Output_MidX(int step);

extern void Output_MidY(int step);

extern void Output_Diag(int step);

//------------------------Boundary.cpp----------------------

extern void P_Inlet_4_Boundary();

extern void P_Outlet_5_Boundary();

extern void WallShadowC_fBP(Cell_2D &shadowCell);

extern void Wall_3_DS(Face_2D &face);

extern void Wall_3_NEE(Face_2D &face);

extern void Wall_3_BB(Face_2D &face);

extern void fluxCheck(Face_2D const* faceptr);

//----------------------------------------------------------

extern void Update_force(Cell_2D &cell);

extern void Update_force_h(Face_2D &face);

extern void MacroSource(Cell_2D *cellptr);

extern void Update_PseudoForce(Cell_2D &cell);

extern void Update_PseudoForce(Face_2D &face);

extern void Output_xcyc();

//----------------------------------DEBUG---------------------------------------
#ifdef _ZERO_NDEBUG_FLIP

extern void Output_UVP(double const &t);

extern void Output_fBh(Face_2D& face,double t);

extern void Output_gBh(Face_2D& face,double t);

extern void Output_f_hdt(Face_2D& face,double t);

extern void Output_g_hdt(Face_2D& face,double t);

extern void Output_fT(Cell_2D& face,double t);

extern void Output_gT(Cell_2D& face,double t);

extern void Output_fT_Append(Cell_2D &cell,double dt);

extern void Output_gT_Append(Cell_2D &cell,double dt);

extern void Output_f_hdt_Append(Face_2D &face,double dt);

extern void Output_g_hdt_Append(Face_2D &face,double dt);

// extern void Output_fFlux_Append(Cell_2D &cell,double dt);

// extern void Output_gFlux_Append(Cell_2D &cell,double dt);

extern void Output_phi_Bh(Face_2D &face,double t);

#endif

//-------------------------------------GradScheme.cpp---------------------------
extern void LeastSquareDebug();

extern void Grad_VS_LS(Cell_2D *center);

extern void Grad_VS_6points(Cell_2D *center);

extern void Grad_VS_4points(Cell_2D *center);
//------------------------------------------------------------------------------
void Forecast_DVDF_Tilde(Cell_2D &cell);

void StrangSplitting_Source(Cell_2D &cell);

void ZeroSource(Cell_2D &cell);

void Update_DVDF_BP(Cell_2D& cell);

void Update_DVDF_h(Face_2D& face);

void Update_DVDF_Flux_h(Face_2D& face);

void Update_Flux(Face_2D &face);

void Update_BoundFlux(Face_2D &face);

void Update_DVDF_Tilde(Cell_2D& cell);

void Zero_GradBarPlus(Cell_2D& cell);

void Flux_2D(Face_2D &face);

void Flux_2D_Limiter(Face_2D &face);

void Update_Residual(int step);

void UpdateL2Error(int step);

void Update_SumRho(int step);

void Update_SumEk(int step);
//--------------------------------DmVn.cpp--------------------------
extern void Update_DVDF_Source(Cell_2D &cell);

extern void Update_DVDF_Source_h(Face_2D &face);

extern void IntegralShearStress();

//-------------------------------------------------------------------
auto startLoop = std::chrono::system_clock::now();
auto endLoop = std::chrono::system_clock::now();

auto startA = std::chrono::system_clock::now();
auto endA = std::chrono::system_clock::now();

auto startB = std::chrono::system_clock::now();
auto endB = std::chrono::system_clock::now();

auto startC = std::chrono::system_clock::now();
auto endC = std::chrono::system_clock::now();

void DUGKS2DSolver()
{
	step = 0;
	Output_Flowfield(0.0,0);
	// Output_MidX(step);
	// Output_RTIlocation(0,1);
	// Output_RTIlocation(0,Nx/2);
	//
	printSplitLine();
	cout << "iteration start    ThreadNum : "<<ThreadNum<<endl;
	printSplitLine(nl);
	//
	omp_set_num_threads(ThreadNum);
#pragma omp parallel
{
//!------------------------loop selection---------------------------

	// while(step < End_Step)
//
	while(ResidualPer1k > RESIDUAL)
{
//---------------------------------------------------------
	#ifdef _ARK_STRANGSPLIT_FLIP
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		StrangSplitting_Source(CellArray[n]);
		ZeroSource(CellArray[n]);
	}
	#endif
	//
	// #pragma omp single
	// if(step%10 == 0)
	// {
	// startA = std::chrono::system_clock::now();
	// }
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		Update_DVDF_Eq(CellArray[n]);
		Update_DVDF_BP(CellArray[n]);
	}
	// #pragma omp single
	// if(step%10 == 0)
	// {
	// endA = std::chrono::system_clock::now();
	// cout <<step<<fs<<"update_BarPlus"<<fs
	// <<std::chrono::duration_cast<std::chrono::milliseconds>(endA-startA).count()
	// <<'\n';
	// }
//-------------------Update-shadow fBP--------------------------------
	#ifdef _Wall_3_BCs_FLIP
	#pragma omp for schedule(guided)
	for(int n = 0;n < WallFaceNum;++n)
		WallShadowC_fBP(WallShadowCA[n]);
	#endif
	#ifdef _P_INLET_4_BCS_FLIP	
		P_Inlet_4_Boundary();
	#endif
	#ifdef _P_OUTLET_5_BCS_FLIP
		P_Outlet_5_Boundary();
	#endif
//-------------------------------Update Grad DVDF_BarPlus-------------------------------
	#ifdef _FLUX_SCHEME_UW_ARK
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		//Grad_VS_LS(&CellArray[n]);
		Grad_VS_6points(&CellArray[n]);
		//Grad_VS_4points(&CellArray[n]);
		//Zero_GradBarPlus(CellArray[n]);
	}
	#endif
//-------------------------------Flux-------------------------------------
//-------------------------------Interior Face-----------------------------
	// #pragma omp single
	// if(step%10 == 0)
	// {
	// startB = std::chrono::system_clock::now();
	// }

	#ifdef _ARK_LIMITER_FLIP
		#pragma omp for schedule(guided)
		for(int n = 0;n < InteriorFaceNum;++n)
		{
			Flux_2D_Limiter(*InteriorFaceA[n]);
		}
		#pragma omp for schedule(guided)
		for(int n = 0;n < BoundFaceNum;++n)
			Flux_2D(*BoundFaceA[n]);
	#else
		#pragma omp for schedule(guided)
		for(int n = 0;n < InteriorFaceNum;++n)
		{
			Flux_2D(*InteriorFaceA[n]);
		}
		#ifdef _PERIODIC_12_8_BCs_FLIP
		#pragma omp for schedule(guided)
		LoopPSB(PeriodicFaceNum)
		{
			Flux_2D(*PeriodicFaceA[nB]);
		}
		#endif
	#endif
//----------------------Wall Face---------------------------
	#ifdef _Wall_3_BCs_FLIP
		#pragma omp for schedule(guided)
		for(int n = 0;n < WallFaceNum;++n)
		{
			#ifdef _Wall_3_BCs_DS
			Wall_3_DS(*WallFaceA[n]);
			#endif
//
			#ifdef _Wall_3_BCs_NEE
			Wall_3_NEE(*WallFaceA[n]);
			#endif
//
			#ifdef _Wall_3_BCs_BB
			Flux_2D(*WallFaceA[n]);
			Wall_3_BB(*WallFaceA[n]);
			#endif
		}
	#endif
// 
	// #pragma omp single
	// if(step%10 == 0)
	// {
	// endB = std::chrono::system_clock::now();
	// cout <<step<<fs<<"update_Flux"<<fs
	// <<std::chrono::duration_cast<std::chrono::milliseconds>(endB-startB).count()
	// <<'\n';
	// }

	#pragma omp for schedule(guided) 
	LoopPS(Faces)
	{
		Update_DVDF_Flux_h(FaceArray[n]);
	}
	#pragma omp for schedule(guided)
	for(int n = 0;n < WallFaceNum;++n)
	{
		fluxCheck(WallFaceA[n]);
	}

	// #pragma omp single
	// if(step%10 == 0)
	// {
	// startC = std::chrono::system_clock::now();
	// }
//--------------------auxilary DDF : DVDF tilde----------------------------
	#pragma omp for schedule(guided)
	LoopPS(Cells)
		Update_DVDF_Tilde(CellArray[n]);
//---------------------------------macro variables---------------------
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		Update_MacroVar(CellArray[n]);
	}
	#ifdef _ARK_ALLENCAHN_FLIP
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		MacroSource(&CellArray[n]);
	}
	#endif
//!-------------------------pseudopotential force-------------------------
	#ifdef _ARK_PSEUDOPSI_FLIP
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		Update_PseudoForce(CellArray[n]);
	}
	#endif
	//
	#ifdef _ARK_STRANGSPLIT_FLIP
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		StrangSplitting_Source(CellArray[n]);
	}
	#endif

	// #pragma omp single
	// if(step%10 == 0)
	// {
	// endC = std::chrono::system_clock::now();
	// cout <<step<<fs<<"update_Tilde & MacroVars"<<fs
	// <<std::chrono::duration_cast<std::chrono::milliseconds>(endC-startC).count()
	// <<'\n';
	// }

	//
	#pragma omp single
	{
		++step;
		if(step%ConvergenceControl == 0)
		{
			Update_SumRho(step);
			//Update_SumEk(step);
			// Output_RTIlocation(step,1);
			// Output_RTIlocation(step,Nx/2);
			#ifdef _OUTPUT_L2NORM_ERROR_FLIP
			UpdateL2Error(step);
			#endif
		}
		if(step%ResidualControl == 0)
		{
			Update_Residual(step);
		}
	}
  }
}
	#ifdef _ARK_THERMAL_FLIP
	IntegralShearStress();
	#endif

	Output_Flowfield((step)*::dt,step);
	// Output_Diag(step);
	// Output_MidX(step);
//
	printSplitLine();
	cout <<"iteration End"<<nl;
	printSplitLine(nl);
}

void Zero_GradBarPlus(Cell_2D &cell)
{
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.BarP_x[k] = 0.0;
		cell.h.BarP_y[k] = 0.0;
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.BarP_x[k] = 0.0;
		cell.f.BarP_y[k] = 0.0;
		#endif
//isothermal flip
		#ifdef _ARK_THERMAL_FLIP
		cell.g.BarP_x[k] = 0.0;
		cell.g.BarP_y[k] = 0.0;
		#endif
	}
}
void ZeroSource(Cell_2D &cell)
{
	cell.MsQ().Fx = 0.0;
	cell.MsQ().Fy = 0.0;

	LoopVS(Q)
	{
		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.So[k] = 0.0;
		#endif
	}
}
void Forecast_DVDF_Tilde(Cell_2D &cell)
{
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.Tilde[k] += ::hDt*cell.h.So[k];
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.Tilde[k] += ::hDt*cell.f.So[k];
		#endif
	}
}
void StrangSplitting_Source(Cell_2D &cell)
{
	Update_DVDF_Eq(cell);
	Update_DVDF_Source(cell);
	Forecast_DVDF_Tilde(cell);

	#ifdef _ARK_MOMENTUM_FLIP
	cell.MsQ().U += ::hDt*cell.MsQ().Fx/cell.MsQ().Rho;
	cell.MsQ().V += ::hDt*cell.MsQ().Fy/cell.MsQ().Rho;
	#endif
}
void Update_DVDF_BP(Cell_2D& cell)
{
	// for(int i = 0;i < DV_Qu;++i)
	// for(int j = 0;j < DV_Qv;++j)
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.BarP[k] = cell.h.aBP*cell.h.Tilde[k] + cell.h.bBP*cell.h.Eq[k]
							+ cell.h.cBP*cell.h.So[k];
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.BarP[k] = cell.f.aBP*cell.f.Tilde[k] + cell.f.bBP*cell.f.Eq[k]
						+ cell.f.cBP*cell.f.So[k];
		#endif
		//isothermal flip
		#ifdef _ARK_THERMAL_FLIP
		cell.g.BarP[k] = cell.g.aBP*cell.g.Tilde[k] + cell.g.bBP*cell.g.Eq[k];
		#endif
	}
}
void Update_DVDF_h(Face_2D& face)
{
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.hDt[k] = face.h.ah*face.h.BhDt[k] + face.h.bh*face.h.EqhDt[k]
						+ face.h.ch*face.h.SohDt[k];
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		face.f.hDt[k] = face.f.ah*face.f.BhDt[k] + face.f.bh*face.f.EqhDt[k]
						+ face.f.ch*face.f.SohDt[k];
		#endif
		//isothermal flip
		#ifdef _ARK_THERMAL_FLIP
		face.g.hDt[k] = face.g.ah*face.g.BhDt[k] + face.g.bh*face.g.EqhDt[k];
		#endif
	}
}
void Update_DVDF_Flux_h(Face_2D &face)
{
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.hDt[k] *= face.xi_n_dS[k];
		#endif
		//momentum
		#ifdef _ARK_MOMENTUM_FLIP
		face.f.hDt[k] *= face.xi_n_dS[k];
		#endif
		//isothermal flip
		#ifdef _ARK_THERMAL_FLIP
		face.g.hDt[k] *= face.xi_n_dS[k];
		#endif
	}
}
inline
double VenkatakrishnanExpression(double a,double b)
{
	return (a*a + 2*a*b + Kh3)/(a*a + 2*b*b + a*b + Kh3);
} 
void VenkatakrishnanFluxLimiter(Cell_2D &cell,int const k)
{
	#ifdef _ARK_MOMENTUM_FLIP
	double GradfBPDotDelta[4] = {0},LfBP[4] = {0};
	double MaxfBP = cell.f.BarP[k],MinfBP = cell.f.BarP[k];
	double MinfBPLimiter = 1;
	#endif
//
	// double MaxfBP = -1E+5,MinfBP = -MaxfBP, MaxgBP = -1E+5,MingBP = -MaxgBP;
	#ifdef _ARK_THERMAL_FLIP
	double GradgBPDotDelta[4] = {0},LgBP[4] = {0};
	double MaxgBP = cell.g.BarP[k],MingBP = cell.g.BarP[k];
	double MingBPLimiter = 1;
	#endif
//	
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		if(cell.Cell_C[iFace]->f.BarP[k] > MaxfBP)
			MaxfBP = cell.Cell_C[iFace]->f.BarP[k];
		if(cell.Cell_C[iFace]->f.BarP[k] < MinfBP)
			MinfBP = cell.Cell_C[iFace]->f.BarP[k];
		#endif
		//isothermal
		#ifdef _ARK_THERMAL_FLIP
		if(cell.Cell_C[iFace]->g.BarP[k] > MaxgBP)
			MaxgBP = cell.Cell_C[iFace]->g.BarP[k];
		if(cell.Cell_C[iFace]->g.BarP[k] < MingBP)
			MingBP = cell.Cell_C[iFace]->g.BarP[k];
		#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		GradfBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.f.BarP_x[k]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.f.BarP_y[k];					        
		if(GradfBPDotDelta[iFace] > 0)
		{
			LfBP[iFace]
			=
			VenkatakrishnanExpression(MaxfBP-cell.f.BarP[k],GradfBPDotDelta[iFace]);
		}
		else if(GradfBPDotDelta[iFace] < 0)
		{
			LfBP[iFace]
			=
			VenkatakrishnanExpression(MinfBP-cell.f.BarP[k],GradfBPDotDelta[iFace]);
		}
		else
			LfBP[iFace] = 1;
		#endif
		//isothermal
		#ifdef _ARK_THERMAL_FLIP
		GradgBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.g.BarP_x[k]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.g.BarP_y[k];
		if(GradgBPDotDelta[iFace] > 0)
		{
			LgBP[iFace]
			=  
			VenkatakrishnanExpression(MaxgBP-cell.g.BarP[k],GradgBPDotDelta[iFace]);
		}
		else if(GradgBPDotDelta[iFace] < 0)
		{
			LgBP[iFace]
			=
			VenkatakrishnanExpression(MingBP-cell.g.BarP[k],GradgBPDotDelta[iFace]);
		}
		else
			LgBP[iFace] = 1;
		#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		if(LfBP[iFace] < MinfBPLimiter) MinfBPLimiter = LfBP[iFace];
		#endif
		//isothermal
		#ifdef _ARK_THERMAL_FLIP
		if(LgBP[iFace] < MingBPLimiter) MingBPLimiter = LgBP[iFace];
		#endif
	}
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	cell.fBPLimiter = MinfBPLimiter;
	#endif
	//isothermal
	#ifdef _ARK_THERMAL_FLIP
	cell.gBPLimiter = MingBPLimiter;
	#endif
}
void UW_Interior_DVDF_Bh_Limiter(Face_2D& face,Cell_2D* ptr_C,int const k)
{
	double dx = face.xf - hDt*xi_u[k] - ptr_C->xc;
	double dy = face.yf - hDt*xi_v[k] - ptr_C->yc;
	VenkatakrishnanFluxLimiter(*ptr_C,k);
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = ptr_C->h.BarP[k]
				   + (dx*ptr_C->h.BarP_x[k] + dy*ptr_C->h.BarP_y[k]);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.f.BhDt[k] = ptr_C->f.BarP[k]
	               + ptr_C->fBPLimiter*(dx*ptr_C->f.BarP_x[k] + dy*ptr_C->f.BarP_y[k]);
	#endif
	//isothermal flip
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = ptr_C->g.BarP[k]
				   + ptr_C->gBPLimiter*(dx*ptr_C->g.BarP_x[k] + dy*ptr_C->g.BarP_y[k]);
	#endif
}
void UW_Interior_DVDF_Bh(Face_2D& face,Cell_2D* ptr_C,int const k)
{
	double dx = face.xf - hDt*xi_u[k] - ptr_C->xc;
	double dy = face.yf - hDt*xi_v[k] - ptr_C->yc;
	//
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = ptr_C->h.BarP[k]
					  + (dx*ptr_C->h.BarP_x[k] + dy*ptr_C->h.BarP_y[k]);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP				  
	face.f.BhDt[k] = ptr_C->f.BarP[k]
					  + (dx*ptr_C->f.BarP_x[k] + dy*ptr_C->f.BarP_y[k]);
	#endif
	//isothermal flip
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = ptr_C->g.BarP[k]
					  + (dx*ptr_C->g.BarP_x[k] + dy*ptr_C->g.BarP_y[k]);
	#endif
}
void UW_Interior_DVDF_Bh(Face_2D& face,int const k)
{
	Cell_2D const *ow = face.owner, *ne = face.neigh; 
	double dxow = face.xf - hDt*xi_u[k] - ow->xc;
	double dyow = face.yf - hDt*xi_v[k] - ow->yc;

	double dxne = face.xf - hDt*xi_u[k] - ne->xc;
	double dyne = face.yf - hDt*xi_v[k] - ne->yc;
	#ifdef _ARK_ALLENCAHN_FLIP
	// face.h.BhDt[k] = 0.5*(face.owner->h.BarP[k]+face.neigh->h.BarP[k]);
	face.h.BhDt[k] = (ow->h.BarP[k] + (dxow*ow->h.BarP_x[k] + dyow*ow->h.BarP_y[k])
	+ ne->h.BarP[k] + (dxne*ne->h.BarP_x[k] + dyne*ne->h.BarP_y[k]))/2;
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP			  
	face.f.BhDt[k] = (ow->f.BarP[k] + (dxow*ow->f.BarP_x[k] + dyow*ow->f.BarP_y[k])
	+ ne->f.BarP[k] + (dxne*ne->f.BarP_x[k] + dyne*ne->f.BarP_y[k]))/2;
	// face.f.BhDt[k] = 0.5*(face.owner->f.BarP[k]+face.neigh->f.BarP[k]);
	#endif
	//isothermal flip
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = 0.5*(face.owner->g.BarP[k] + face.neigh->g.BarP[k]);
	#endif
}
void CD_Interior_DVDF_Bh(Face_2D &face,int const k)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	double _hBP_xF = 1000,_hBP_yF = 1000;
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	double _fBP_xF = 1000,_fBP_yF = 1000;
	#endif
	//!isothermal
	#ifdef _ARK_THERMAL_FLIP
	double _gBP_xF = 1000,_gBP_yF = 1000;
	#endif
	//!central difference
	if(0.0 == face.Vx)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		//_hBP_xF = 0.5*(face.owner->h.BarP_x[k] + face.neigh->h.BarP_x[k]);
		_hBP_xF = (0.5*(face.faceCells[3]->h.BarP[k] - face.faceCells[1]->h.BarP[k])
				+ 0.5*(face.faceCells[2]->h.BarP[k] - face.faceCells[0]->h.BarP[k]))
				/face._dx;
		_hBP_yF = (face.owner->h.BarP[k] - face.neigh->h.BarP[k])/face._dy;
		#endif
		//
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		_fBP_xF = (0.5*(face.faceCells[3]->f.BarP[k] - face.faceCells[1]->f.BarP[k])
				+ 0.5*(face.faceCells[2]->f.BarP[k] - face.faceCells[0]->f.BarP[k]))
				/face._dx;
		_fBP_yF = (face.owner->f.BarP[k] - face.neigh->f.BarP[k])/face._dy;
		#endif
		//!isothemal
		#ifdef _ARK_THERMAL_FLIP
		//_gBP_xF =  0.5*(face.owner->g.BarP_x[k] + face.neigh->g.BarP_x[k]);
		_gBP_xF = (0.5*(face.faceCells[3]->g.BarP[k] - face.faceCells[1]->g.BarP[k])
				+ 0.5*(face.faceCells[2]->g.BarP[k] - face.faceCells[0]->g.BarP[k]))
				/face._dx;
		_gBP_yF = (face.owner->g.BarP[k] - face.neigh->g.BarP[k])/face._dy;
		#endif
	}
	else if(0.0 == face.Vy)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		//_hBP_yF = 0.5*(face.owner->h.BarP_y[k] + face.neigh->h.BarP_y[k]);
		_hBP_yF = (0.5*(face.faceCells[3]->h.BarP[k] - face.faceCells[1]->h.BarP[k])
				+ 0.5*(face.faceCells[2]->h.BarP[k] - face.faceCells[0]->h.BarP[k]))
				/face._dy;
		_hBP_xF = (face.owner->h.BarP[k] - face.neigh->h.BarP[k])/face._dx;
		#endif
		//
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		_fBP_yF = (0.5*(face.faceCells[3]->f.BarP[k] - face.faceCells[1]->f.BarP[k])
				+ 0.5*(face.faceCells[2]->f.BarP[k] - face.faceCells[0]->f.BarP[k]))
				/face._dy;
		_fBP_xF = (face.owner->f.BarP[k] - face.neigh->f.BarP[k])/face._dx;
		#endif
		//
		#ifdef _ARK_THERMAL_FLIP
		_gBP_yF = (0.5*(face.faceCells[3]->g.BarP[k] - face.faceCells[1]->g.BarP[k])
				+ 0.5*(face.faceCells[2]->g.BarP[k] - face.faceCells[0]->g.BarP[k]))
				/face._dy;
		_gBP_xF = (face.owner->g.BarP[k] - face.neigh->g.BarP[k])/face._dx;
		#endif
	}
	else
	{
		cout <<"wrong interface"<<endl;
		getchar();
	}
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[k] = 0.5*(face.owner->h.BarP[k] + face.neigh->h.BarP[k])
			- hDt*(_hBP_xF*xi_u[k] + _hBP_yF*xi_v[k]);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	face.f.BhDt[k] = 0.5*(face.owner->f.BarP[k] + face.neigh->f.BarP[k])
			- hDt*(_fBP_xF*xi_u[k] + _fBP_yF*xi_v[k]);
	#endif
	//!isothemal
	#ifdef _ARK_THERMAL_FLIP
	face.g.BhDt[k] = 0.5*(face.owner->g.BarP[k] + face.neigh->g.BarP[k])
			- hDt*(_gBP_xF*xi_u[k] + _gBP_yF*xi_v[k]);
	#endif	
}
void Update_DVDF_Bh(Face_2D &face)
{
	LoopVS(Q)
	{
		#if defined _FLUX_SCHEME_CD_ARK
			CD_Interior_DVDF_Bh(face,k);
		#elif defined _FLUX_SCHEME_UW_ARK
		if(face.xi_n_dS[k] > 0)
			UW_Interior_DVDF_Bh(face,face.owner,k);
		else if(face.xi_n_dS[k] < 0)
			UW_Interior_DVDF_Bh(face,face.neigh,k);
		else
			UW_Interior_DVDF_Bh(face,k);
		#else 
			exit(-1);
		#endif
	}
}
void Update_DVDF_Bh_Limiter(Face_2D &face)
{
	LoopVS(Q)
	{
		#if defined _FLUX_SCHEME_CD_ARK
			CD_Interior_DVDF_Bh(face,k);
		#elif defined _FLUX_SCHEME_UW_ARK	
		if(face.xi_n_dS[k] >= 0)
			UW_Interior_DVDF_Bh_Limiter(face,face.owner,k);
		else
			UW_Interior_DVDF_Bh_Limiter(face,face.neigh,k);
		#else 
			exit(-1);
		#endif
	}
}
void Flux_2D_Limiter(Face_2D &face)
{	
	Update_DVDF_Bh_Limiter(face);
	Update_MacroVar_h(face);
	Update_DVDF_Eqh(face);
	Update_DVDF_h(face);
//------------------------------------------DEBUG---------------------------------
}
void Flux_2D(Face_2D &face)
{	
	Update_DVDF_Bh(face);
	Update_MacroVar_h(face);
	Update_DVDF_Eqh(face);
	Update_DVDF_h(face);
}
//-------------------------------------------------------------------------------
void Update_DVDF_Tilde(Cell_2D &cell)
{
	LoopVS(Q)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		double hFluxSum = 0.0;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		double fFluxSum = 0;
		#endif
		//isothermal flip
		#ifdef _ARK_THERMAL_FLIP
		gFluxSum = 0;
		#endif
		for(int iF = 0;iF < cell.celltype;++iF)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			hFluxSum += cell.signFlux[iF]*cell.Face_C[iF]->h.hDt[k];
			#endif
			//!momentum
			#ifdef _ARK_MOMENTUM_FLIP
			fFluxSum += cell.signFlux[iF]*cell.Face_C[iF]->f.hDt[k];
			#endif
			//isothermal flip
			#ifdef _ARK_THERMAL_FLIP
			gFluxSum += cell.signFlux[iF]*cell.Face_C[iF]->g.hDt[k];
			#endif
		}
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.Tilde[k] = aTP*cell.h.BarP[k] - bTP*cell.h.Tilde[k]
							+ cell.DtSlashVolume*hFluxSum;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		cell.f.Tilde[k] = aTP*cell.f.BarP[k] - bTP*cell.f.Tilde[k]
							+ cell.DtSlashVolume*fFluxSum;
		#endif
//isothermal flip
		#ifdef _ARK_THERMAL_FLIP
		cell.g.Tilde[k] = aTP*cell.g.BarP[k] - bTP*cell.g.Tilde[k]
							+ cell.DtSlashVolume*gFluxSum;
		#endif
	}
}
void Update_SumRho(int step)
{
	::SumRho = 0.0;
	::SumT = 0.0;
	LoopPS(Cells)
	{
		::SumRho += CellArray[n].MsQ().Phi;
		#ifdef _ARK_THERMAL_FLIP 
		::SumT += CellArray[n].MsQ().T;
		#endif
	}
//	::SumRho -= CarCellArray[64][65]->MsQ().Rho;
	Output_SumRho(step*dt);
}
void Update_SumEk(int step)
{
	::SumEk = 0.0;
	LoopPS(Cells)
	{
		if(::SumEk < CellArray[n].MsQ().SqUV())
		   ::SumEk = CellArray[n].MsQ().SqUV();
	}
	Output_SumEk(step*::dt);
}
void Update_Residual(int step)
{
#ifndef _ARK_ENDTIME_FLIP
//---------------------------density residual-------------------------
	double SumRho = 0.0,SumdRho = 0.0;
	double dRho = 0.0;
	LoopPS(Cells)
	{
		dRho = CellArray[n].MsQ().Rho - CellArray[n].MsQ().Rho_1k;
		SumdRho += dRho*dRho;
		SumRho += CellArray[n].MsQ().Rho*CellArray[n].MsQ().Rho;
		CellArray[n].MsQ().Rho_1k = CellArray[n].MsQ().Rho;
	}
	ResidualPer1k = sqrt(SumdRho/(SumRho + 1.0E-30));
//---------------------------velocity residual--------------------------
	// double SumUV = 0.0,Sumdudv = 0.0;
	// double du = 0.0, dv = 0.0;
	// for(int i = 0;i < Cells;++i)
	// {
	// 	du = CellArray[i].MsQ().U - CellArray[i].MsQ().U_1k;
	// 	dv = CellArray[i].MsQ().V - CellArray[i].MsQ().V_1k;
	// 	Sumdudv += du*du + dv*dv;
	// 	SumUV += CellArray[i].MsQ().SqUV();
	// 	CellArray[i].MsQ().U_1k = CellArray[i].MsQ().U;
	// 	CellArray[i].MsQ().V_1k = CellArray[i].MsQ().V;
	// }
	// ResidualPer1k = sqrt(Sumdudv/(SumUV + 1.0E-30));
//----------------------------------------------------------------------

	Output_Residual(step*dt,ResidualPer1k);
	if(step%writeFileControl == 0)
	{
		Output_Flowfield(step*dt, step);
		// Output_Diag(step);
		// Output_MidX(step);
	}
#endif

#if defined _ARK_ENDTIME_FLIP
	if
	(
		PhaseFieldAC::iT == step 
	||  PhaseFieldAC::iT/8 == step
	||  PhaseFieldAC::iT/4 == step
	||  PhaseFieldAC::iT*3/8 == step
	||  PhaseFieldAC::iT/2 == step
	||  PhaseFieldAC::iT*5/8 == step
	||  PhaseFieldAC::iT*3/4 == step
	||  PhaseFieldAC::iT*7/8 == step
	)
	{
		Output_Flowfield(step*dt, step);
		// Output_MidX(step);
	}
#endif
	//!suppress print to screen on server;
#ifndef _ARK_NOHUP_FLIP
	endLoop = std::chrono::system_clock::now();
	cout << setiosflags(ios::scientific) << setprecision(12);
	cout <<step <<"    "<<step*dt<<"    "<<SumRho<<"    "<<::SumT<<"    "<<ResidualPer1k<<fs
	<<std::chrono::duration_cast<std::chrono::milliseconds>(endLoop-startLoop).count()
	<<'\n';
	startLoop = std::chrono::system_clock::now();
#endif
}
void UpdateL2Error(int step)
{
	double L2_uv, L2_p=0;
	Output_L2Norm(step*dt,L2_uv,L2_p);
}