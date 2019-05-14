#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <string>
#include "DUGKSDeclaration.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::map;
using std::ios;
using std::setprecision;

int const Out_precision = 12;

extern char const DmQnName[];

extern string caseName;

extern void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p);

extern void TaylorCouetteAnalytical(double x,double y,double &u_A);

extern void AnalyticalForceDrivenTG(double x,double y,double &u_A, double &v_A,double &p_A);

extern void LayeredPoiseuilleAnalytical(double const xyz,double &u_A);

extern double AnalyticalPhiAC_Drop(double const x,double const y);


void oss_XXX(ostringstream& oss_r,const string &folder,const string &suffix,double const &t)
{
	oss_r << "../FlowField/"<<folder<<"/Time" <<t<<"."<<suffix;
}
void FileOpen(ofstream &OutFile_XXX,ostringstream &oss_XXX,string const &s)
{
	OutFile_XXX.open(oss_XXX.str().c_str());
	if(!OutFile_XXX)
	{
		cout <<"  "<< ("OutFile_" + s + " open failed") << endl; 
		printErrorMessage(__LINE__,__FILE__,__func__);
		getchar();
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
void FileOpenAppend(ofstream &OutFile_XXX,ostringstream &oss_XXX,string const &s)
{
	OutFile_XXX.open(oss_XXX.str().c_str(),ios::app);
	if(!OutFile_XXX)
	{
		cout <<"  "<< ("OutFile_" + s + " open failed") << endl; 
		printErrorMessage(__LINE__,__FILE__,__func__);
		getchar();
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
void OutputCase()
{
	ofstream OutFile_Case("../FlowField/Case.ark");
	if(!OutFile_Case)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"Case.ark open failed."<<endl;
		getchar();
		return;
	}
	OutFile_Case <<caseName<<"    =    Case Name"<<endl
				 <<_MESHFILE_NAME_ARK<<"    =    Mesh File"<<endl
				 <<_MESHTYPE_ARK<<"    =    Mesh Type"<<endl
				 <<_FLUX_SCHEME_ARK<<"    =    Flux Scheme"<<endl
				 <<_BC_ARK<<"    =    Boundary Condition"<<endl
				 #ifdef _Wall_3_BCs_FLIP
				 <<_Wall_3_BCs_MODEL<<"    =    Wall BC Type"<<endl
				 #endif
				 << DmQnName<<"    =    velocity model"<<endl
				 #ifdef _ARK_FORCE_FLIP
				 <<_ARK_FORCE_MODEL<<"    =    Force model"<<endl
				 <<_ARK_FORCE_H_Liang<<"    =    Liang Model"<<endl
				 #endif
				 <<_EOS_MODEL_ARK<<"    =    EoS model"<<endl
				 // <<left<<fs<<"="<<fs<<"left"<<endl
				 // <<right<<fs<<"="<<fs<<"right"<<endl
				 // <<top<<fs<<"="<<fs<<"top"<<endl
				 // <<bottom<<fs<<"="<<fs<<"bottom"<<endl
				 <<"#-------------------------------------------"<<endl<<endl
				 #ifdef _ARK_THERMAL_FLIP
				 <<"#--------------thermodynamics------------"<<endl
				 <<Omega0<<"    =    Omega0//sutherland power"<<endl
				 <<Pr<<"    =    Prandtl"<<endl
				 <<nK<<"    =    nK//internal degrees of freedom"<<endl
				 <<Cv<<"    =    specific heat of constant volume"<<endl
				 <<Gamma<<"    =    Gamma//specific heat ratio"<<endl
				 <<"#-------------------------------------------"<<endl<<endl
				 #endif
				 <<"#----------------Reference----------------"<<endl
				 <<Q<<"    =    Velocity Space"<<endl
				 <<NL<<"    =    NL"<<endl
				 <<ChLength<<"    =    Characteristic Length"<<endl
				 <<MinL<<"    =    MinL"<<endl
				 <<Kn<<"    =    Knudsen"<<endl
				 <<Ma<<"    =    Mach"<<endl
				 <<Re<<"    =    Re"<<endl
				 <<U0<<"    =    U0"<<endl
				 <<R0<<"    =    R0"<<endl
				 <<T0<<"    =    T0"<<endl
				 <<Mu0<<"    =    Mu0(dynamic viscosity)"<<endl
				 <<Nu0<<"    =    Nu0(kinetic viscosity)"<<endl
				 <<Tau0<<"    =    tau0"<<endl
				 <<Lambda0<<"    =    Lambda0"<<endl
				 <<CFL<<"    =    CFL"<<endl
				 <<dt<<"    =    dt"<<endl				 
				 <<dt/Tau0<<"    =    dt/tau0"<<endl
				 <<MaSpan<<"    =    MaSpan"<<endl
				 <<Eta<<"    =    Eta// T of D2V16"<<endl
				 <<"#-------------------------------------------"<<endl<<endl
				 <<"---------------write step----------------"<<endl
				 <<::End_Step<<"    =    End_Step"<<endl
				 <<::ConvergenceControl<<"    =    ConvergenceControl"<<endl
				 <<::ResidualControl<<"    =    ResidualControl"<<endl
				 <<::writeFileControl<<"    =    writeFileControl"<<endl
				 <<::RESIDUAL<<"    =    Residual"<<endl
				<<"#-------------------------------------------"<<endl<<endl
				 #ifdef _ARK_ALLENCAHN_FLIP
				 <<"-------------Phase Field--------------"<<endl
				 <<PhaseFieldAC::Cn<<"    =    Cahn number"<<endl
				 <<PhaseFieldAC::Pe<<"    =    Peclet number"<<endl
				 <<PhaseFieldAC::At<<"    =    Atwood number"<<endl
				 <<PhaseFieldAC::Eo<<"    =    Eotvos(Bond) number"<<endl
				 <<PhaseFieldAC::Mo<<"    =    Morton number"<<endl
				 <<PhaseFieldAC::ReMP<<"    =    Reynolds number"<<endl
				 <<PhaseFieldAC::La<<"    =    Laplacian number"<<endl
				 <<PhaseFieldAC::T<<"    =    Reference Time"<<endl
				 <<PhaseFieldAC::iT<<"    =    End_Time"<<endl
				 <<endl
				 <<PhaseFieldAC::M_Phi<<"    =    mobility"<<endl
				 <<PhaseFieldAC::TauMass<<"    =    mobility relaxition time"<<endl
				 <<dt/PhaseFieldAC::TauMass<<"    =    dt/tauMass"<<endl
				 <<PhaseFieldAC::wI<<"    =    inteface width"<<endl
				 <<PhaseFieldAC::Sigma<<"    =    Sigma(surface tension coefficient)"<<endl
				 <<PhaseFieldAC::Beta<<"    =    Beta"<<endl
				 <<PhaseFieldAC::Kappa<<"    =    Kappa"<<endl
				 <<PhaseFieldAC::RhoL<<"    =    liquid density"<<endl
				 <<PhaseFieldAC::RhoV<<"    =    vapor density"<<endl
				 <<PhaseFieldAC::MuL<<"    =    liquid dynamic viscosity"<<endl
				 <<PhaseFieldAC::MuV<<"    =    vapor dynamic viscosity"<<endl
				 <<PhaseFieldAC::NuL<<"    =    liquid kinetic viscosity"<<endl
				 <<PhaseFieldAC::NuV<<"    =    vapor kinetic viscosity"<<endl
				 <<PhaseFieldAC::RhoL/PhaseFieldAC::RhoV<<"    =    "<<"density ratio"<<endl
				 <<PhaseFieldAC::MuL/PhaseFieldAC::MuV<<"    =    "<<"dynamic viscosity ratio"<<endl
				 <<PhaseFieldAC::radius/ChLength<<"    =    radius"<<endl
				 <<_ARK_LAPLACIAN_SCHEME<<"    =    laplacian operator"<<endl
				 <<"#-------------------------------------------"<<endl
				 #endif
//
				 #ifdef _ARK_PSEUDOPSI_FLIP
				 <<"-------------Pseudo Potential--------------"<<endl
				 <<PseudoPotentialSC::RhoL <<"    =    liquid density"<<endl
				 <<PseudoPotentialSC::RhoV <<"    =    vapor density"<<endl
				 <<PseudoPotentialSC::K <<"    =    modified parameter"<<endl
				 <<PseudoPotentialSC::wI<<"    =    initial inteface width"<<endl
				 <<PseudoPotentialSC::radius/ChLength<<"    =    radius/ChLength"<<endl
				 <<D2Q9::Alpha <<"    =    Alpha de D2V9PP Model"<<endl
				 <<"#-------------------------------------------"<<endl
				 #endif
				 <<"#---------------------END--------------------"<<endl<<endl;
	OutFile_Case.close();
//
}
double thermoPressure(Cell_2D const &cell)
{
	using PhaseFieldAC::Kappa;
	using PhaseFieldAC::Beta;
	using PhaseFieldAC::aPhi;
	using PhaseFieldAC::PhiL;
	using PhaseFieldAC::PhiV;

	double Phi = cell.MsQ().Phi;
	double SqGradPhi = (cell.MsQ().Rho_x*cell.MsQ().Rho_x
					 +  cell.MsQ().Rho_y*cell.MsQ().Rho_y)/(aPhi*aPhi);
	double pZero = Phi*2*Beta*(Phi-PhiL)*(Phi-PhiV)*(2*Phi-PhiL-PhiV)
					-Beta*(Phi-PhiL)*(Phi-PhiL)*(Phi-PhiV)*(Phi-PhiV);

	return 
	(pZero - Kappa*Phi*cell.MsQ().laplacianPhi + Kappa*SqGradPhi/2
	+ cell.MsQ().p);
}
void TaylorGreen_L2Norm(double const &t,double &L2_uv, double &L2_p)
{
	double u_A, v_A, p_A, du, dv, dp; 
	double  Sumdudv = 0.0,Sumdp = 0.0,Sumuv_A = 0.0,Sump_A = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A,v_A,p_A);
		du = CellArray[i].MsQ().U - u_A;
		dv = CellArray[i].MsQ().V - v_A;
		dp = CellArray[i].MsQ().p - p_A;
		Sumdudv += du*du + dv*dv;
		Sumuv_A += u_A*u_A + v_A*v_A; 
		Sumdp   += dp*dp;
		Sump_A  += p_A*p_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
	L2_p  = sqrt(Sumdp/Sump_A);
}
void TaylorCouette_L2Norm(double &L2_uv)
{
	double u_A, du;
	double  Sumdudv = 0.0,Sumuv_A = 0.0;
	for(int i =0;i < Cells;++i)
	{
		TaylorCouetteAnalytical(CellArray[i].xc,CellArray[i].yc,u_A);
		du = sqrt(CellArray[i].MsQ().SqUV()) - u_A;
		Sumdudv += du*du;
		Sumuv_A += u_A*u_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
}
void LayeredPoiseuille_L2Norm(double &L2_uv)
{
	double u_A, du;
	double  Sumdudv = 0.0,Sumuv_A = 0.0;
	for(int i =0;i < Cells;++i)
	{
		LayeredPoiseuilleAnalytical(CellArray[i].yc,u_A);
		du = CellArray[i].MsQ().U - u_A;
		Sumdudv += du*du;
		Sumuv_A += u_A*u_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
}
void ForceDrivenTaylorGreen_L2Norm(double &L2_uv, double &L2_p)
{
	double u_A, v_A, p_A, du, dv, dp;
	double  Sumdudv = 0.0,Sumdp = 0.0,Sumuv_A = 0.0,Sump_A = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		AnalyticalForceDrivenTG(CellArray[i].xc,CellArray[i].yc,u_A,v_A,p_A);
		du = CellArray[i].MsQ().U - u_A;
		dv = CellArray[i].MsQ().V - v_A;
		dp = CellArray[i].MsQ().p - p_A;
		Sumdudv += du*du + dv*dv;
		Sumuv_A += u_A*u_A + v_A*v_A; 
		Sumdp   += dp*dp;
		Sump_A  += p_A*p_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
	L2_p  = sqrt(Sumdp/Sump_A);
}
void AC_Drop_L2Norm(double &L2_Phi)
{
	double Phi_A = 0.0,dPhi = 0.0;
	double SumdPhi = 0.0,SumPhi_A = 0.0;
	LoopPS(Cells)
	{
		Phi_A = AnalyticalPhiAC_Drop(CellArray[n].xc,CellArray[n].yc);
		dPhi = CellArray[n].MsQ().Phi - Phi_A;
		SumdPhi += dPhi*dPhi;
		SumPhi_A += Phi_A*Phi_A;
	}
	L2_Phi = sqrt(SumdPhi/SumPhi_A);
}
void Output_L2Norm(double const &t,double &L2_uv, double &L2_p)
{	
	//ForceDrivenTaylorGreen_L2Norm(L2_uv,L2_p);
	//TaylorCouette_L2Norm(L2_uv);
	//LayeredPoiseuille_L2Norm(L2_uv);
	AC_Drop_L2Norm(L2_uv);
	ostringstream oss_L2;
	oss_L2 <<"../FlowField/Convergence/L2_uvp_mu"<<Mu0<<"_Re"<<Re<<"_"<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_L2(oss_L2.str().c_str(),ofstream::app);
	if(!OutFile_L2)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_L2 file open failed"<<endl;
		getchar();
	}
	OutFile_L2 << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_L2 << t <<"    "<<L2_uv<<"    "<<L2_p<<'\n';
	OutFile_L2.close();
}
void Output_SumRho(double t)
{
	ostringstream oss_SumRho;
	oss_SumRho << "../FlowField/Convergence/SumRho_mu"<<Mu0<<"_Re"<<Re<<"_"<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_SumRho(oss_SumRho.str().c_str(),ofstream::app);
	if(!OutFile_SumRho)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_SumRho file open failed"<<endl;
		getchar();
	}
	OutFile_SumRho << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_SumRho <<  t <<"    "<<SumRho<<endl;
	OutFile_SumRho.close();
}
void Output_SumEk(double t)
{
	ostringstream oss_SumEk;
	oss_SumEk <<"../FlowField/Convergence/SumEk_mu"<<Mu0<<"_Re"<<Re<<"_"<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_SumEk(oss_SumEk.str().c_str(),ofstream::app);
	if(!OutFile_SumEk)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_SumEk file open failed"<<endl;
		exit(-1);
	}
	OutFile_SumEk << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_SumEk <<  t <<"    "<<SumEk<<endl;
	OutFile_SumEk.close();
}
double locateInterface(int const X,int &Facek, int &anotherFacek)
{
	double dPhi = 1;

	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D const &cell = *CarCellArray[X][k];
		if(fabs(cell.MsQ().Phi - 0.5) < dPhi)
		{
			dPhi = fabs(cell.MsQ().Phi - 0.5);
			Facek = k;
		}
	}

	Cell_2D const &cellIf = *CarCellArray[X][Facek];
	if(cellIf.MsQ().Phi > 0.5)
	{
		anotherFacek = Facek-1;
	}
	else if(cellIf.MsQ().Phi < 0.5)
	{
		anotherFacek = Facek+1;
	}
	else
	{
		return cellIf.yc;
	}

	double yA = CarCellArray[X][Facek]->yc;
	double yB = CarCellArray[X][anotherFacek]->yc;
	double PhiA = CarCellArray[X][Facek]->MsQ().Phi;
	double PhiB = CarCellArray[X][anotherFacek]->MsQ().Phi;
	return ((0.5-PhiB)/(PhiA-PhiB)*yA + (0.5-PhiA)/(PhiB-PhiA)*yB);

}
void Output_RTIlocation(double step,int const X)
{
	using PhaseFieldAC::ReMP;

	ostringstream oss_RTI;
	oss_RTI <<"../FlowField/Convergence/RTI_X"<<X<<"_mu"<<Mu0<<"_Re"<<ReMP<<"_"<<_MESHFILE_NAME_ARK<<".dat";
	
	ofstream OutFile_RTI;

	FileOpenAppend(OutFile_RTI,oss_RTI,oss_RTI.str());

	int ku,kl;

	locateInterface(X,ku,kl);

	double yc = locateInterface(X,ku,kl);
	double PhiA = CarCellArray[X][ku]->MsQ().Phi;
	double PhiB = CarCellArray[X][kl]->MsQ().Phi;

	yc = (yc/ChLength - 2.0);

	OutFile_RTI << step/PhaseFieldAC::T <<"    "<<yc<<"    "<<PhiA<<"    "<<PhiB<<endl;
	OutFile_RTI.close();
}
void Output_Residual(double t,double Residual)
{
	ostringstream oss_Residual;
	oss_Residual <<"../FlowField/Convergence/Residual_mu"<<Mu0<<"_Re"<<Re<<"_"<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_Residual(oss_Residual.str().c_str(),ofstream::app);
	if(!OutFile_Residual)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_Residual open failed"<<endl;
		getchar();
	}
	OutFile_Residual << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_Residual << t <<"    "<<Residual<<endl;
	OutFile_Residual.close();
}
void writeMidXHead
(
	ofstream &MidFile,
	const string &Name,
	string const &Title,
	double const& MQ
)
{
	MidFile<<"variables = Y,"<<Name<<nl
		   <<"zone T = "<<Title<<MQ<<" "
		   <<"I = "<<Ny<<endl;
}
void writeMidYHead
(
	ofstream &MidFile,
	const string &Name,
	string const &Title,
	double const& MQ
)
{
	MidFile<<"variables = X,"<<Name<<nl
		   <<"zone T = "<<Title<<MQ<<" "
		   <<"I = "<<Nx<<endl;
}
void Output_MidX(int step)
{
	using PhaseFieldAC::RhoL;
	using PhaseFieldAC::RhoV;
	using PhaseFieldAC::M_Phi;
	int const midX = Nx/2;
	//int const mid = 1;
	ostringstream oss_MidX;
	oss_MidX <<"../FlowField/Convergence/"<<"MidX_"<<"Ma"<< Ma<<"_"
					<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"ratio_"<<RhoL/RhoV<<"step"<<step<<".dat";
	ofstream OutFile_MidX(oss_MidX.str().c_str());
	if(!OutFile_MidX)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_MidX open failed" << endl; 
		getchar();
		exit(-1);
	}
	OutFile_MidX<<std::setiosflags(ios::scientific)<<std::setprecision(20);
	//
	writeMidXHead(OutFile_MidX,"Rho","Rho_M",M_Phi);
	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D &cell = *CarCellArray[midX][k];
		OutFile_MidX<<cell.yc<<fs<<cell.MsQ().Rho<<endl;
	}
	
	writeMidXHead(OutFile_MidX,"PhiMid","Phi_M",M_Phi);
	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D &cell = *CarCellArray[midX][k];
		OutFile_MidX<<cell.yc<<fs<<cell.MsQ().Phi<<endl;
	}
	// writeMidXHead(OutFile_MidX,"Phi0","Phi_M",M_Phi);
	// for(int k = 1;k < Nyp1;++k)
	// {
	// 	Cell_2D &cell = *CarCellArray[1][k];
	// 	OutFile_MidX<<cell.yc<<fs<<cell.MsQ().Phi<<endl;
	// }
	//
	// //
	// writeMidXHead(OutFile_MidX,"U","U_M",M_Phi);
	// for(int k = 1;k < Nyp1;++k)
	// {
	// 	Cell_2D &cell = *CarCellArray[midX][k];
	// 	OutFile_MidX<<cell.yc<<fs<<cell.MsQ().U<<endl;
	// }
	//
	// writeMidXHead(OutFile_MidX,"UA","UA_M",M_Phi);
	// for(int k = 1;k < Nyp1;++k)
	// {
	// 	double u_A;
	// 	Cell_2D &cell = *CarCellArray[midX][k];
	// 	LayeredPoiseuilleAnalytical(cell.yc,u_A);
	// 	OutFile_MidX<<cell.yc<<fs<<u_A<<endl;
	// }
	// writeMidXHead(OutFile_MidX,"dU","dU_M",M_Phi);
	// for(int k = 1;k < Nyp1;++k)
	// {
	// 	double u_A;
	// 	Cell_2D &cell = *CarCellArray[midX][k];
	// 	LayeredPoiseuilleAnalytical(cell.yc,u_A);
	// 	OutFile_MidX<<cell.yc<<fs<<fabs(1-cell.MsQ().U/u_A)<<endl;
	// }
	#ifdef _ARK_MOMENTUM_FLIP
	writeMidXHead(OutFile_MidX,"tau","tau_M",M_Phi);
	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D &cell = *CarCellArray[midX][k];
		OutFile_MidX<<cell.yc<<fs<<cell.f.tau<<endl;
	}
	#endif
	writeMidXHead(OutFile_MidX,"DP","DP_M",M_Phi);
	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D &cell = *CarCellArray[midX][k];
		OutFile_MidX<<cell.yc<<fs<<cell.MsQ().p<<endl;
	}
	
	writeMidXHead(OutFile_MidX,"TP","TP_M",M_Phi);
	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D &cell = *CarCellArray[midX][k];
		double TP = thermoPressure(cell);
		OutFile_MidX<<cell.yc<<fs<<TP<<endl;
	}
	
	writeMidXHead(OutFile_MidX,"U","U_M",M_Phi);
	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D &cell = *CarCellArray[midX][k];
		// /double TP = thermoPressure(cell);
		OutFile_MidX<<cell.yc<<fs<<cell.MsQ().U<<endl;
	}
	
	writeMidXHead(OutFile_MidX,"V","V_M",M_Phi);
	for(int k = 1;k < Nyp1;++k)
	{
		Cell_2D &cell = *CarCellArray[midX][k];
		// /double TP = thermoPressure(cell);
		OutFile_MidX<<cell.yc<<fs<<cell.MsQ().V<<endl;
	}
}
void Output_MidX_reverse(int step)
{
	using PhaseFieldAC::RhoL;
	using PhaseFieldAC::RhoV;
	// int const midX = Nx/2;
	//int const mid = 1;
	ostringstream oss_MidX;
	oss_MidX <<"../FlowField/Convergence/"<<"MidXreverse_"<<"Ma"<< Ma<<"_"
					<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"ratio_"<<RhoL/RhoV<<"step"<<step<<".dat";
	ofstream OutFile_MidX(oss_MidX.str().c_str());
	if(!OutFile_MidX)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_MidX open failed" << endl; 
		getchar();
		exit(-1);
	}
	OutFile_MidX<<std::setiosflags(ios::scientific)<<std::setprecision(20);
	// for(int k = Ny;k > 0;--k)
	// {
	// 	double u_A;
	// 	Cell_2D &cell = *CarCellArray[midX][k];
	// 	LayeredPoiseuilleAnalytical(cell.yc,u_A);
	// 	OutFile_MidX<<cell.yc<<fs<<cell.MsQ().Rho<<fs<<cell.MsQ().U<<fs<<u_A<<endl;
	// }
}
void Output_MidY(int step)
{
	using PhaseFieldAC::RhoL;
	using PhaseFieldAC::RhoV;
	using PhaseFieldAC::M_Phi;
	int const midY = Ny/2;
	//int const mid = 1;
	ostringstream oss_MidY;
	oss_MidY <<"../FlowField/Convergence/"<<"MidY_"<<"Ma"<< Ma<<"_"
					<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"ratio_"<<RhoL/RhoV<<"step"<<step<<".dat";
	ofstream OutFile_MidY(oss_MidY.str().c_str());
	if(!OutFile_MidY)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_MidY open failed" << endl; 
		getchar();
		exit(-1);
	}
	OutFile_MidY<<std::setiosflags(ios::scientific)<<std::setprecision(20);
	// for(int k = 1;k < Nxp1;++k)
	// {
	// 	double u_A;
	// 	Cell_2D &cell = *CarCellArray[k][midY];
	// 	LayeredPoiseuilleAnalytical(cell.yc,u_A);
	// 	OutFile_MidY<<cell.xc<<fs<<cell.MsQ().Rho<<fs<<cell.MsQ().V<<fs<<u_A<<endl;
	// }
	writeMidYHead(OutFile_MidY,"Rho","Rho_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		Cell_2D &cell = *CarCellArray[k][midY];
		OutFile_MidY<<cell.xc<<fs<<cell.MsQ().Rho<<endl;
	}
	//
	writeMidYHead(OutFile_MidY,"DP","DP_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		Cell_2D &cell = *CarCellArray[k][midY];
		OutFile_MidY<<cell.xc<<fs<<cell.MsQ().p<<endl;
	}
	//
	writeMidYHead(OutFile_MidY,"TP","TP_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		Cell_2D &cell = *CarCellArray[k][midY];
		OutFile_MidY<<cell.xc<<fs<<thermoPressure(cell)<<endl;
	}
	//
	writeMidYHead(OutFile_MidY,"U","U_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		Cell_2D &cell = *CarCellArray[k][midY];
		OutFile_MidY<<cell.xc<<fs<<cell.MsQ().U<<endl;
	}
	//
	writeMidYHead(OutFile_MidY,"V","V_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		Cell_2D &cell = *CarCellArray[k][midY];
		OutFile_MidY<<cell.xc<<fs<<cell.MsQ().V<<endl;
	}

}
void Output_Diag(int step)
{
	using PhaseFieldAC::RhoL;
	using PhaseFieldAC::RhoV;
	using PhaseFieldAC::M_Phi;
	//
	ostringstream oss_Diag;
	oss_Diag <<"../FlowField/Convergence/"<<"Diag_"<<"Ma"<< Ma<<"_"
					<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"ratio_"<<RhoL/RhoV<<"step"<<step<<".dat";
	ofstream OutFile_Diag(oss_Diag.str().c_str());
	if(!OutFile_Diag)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_Diag open failed" << endl; 
		getchar();
		exit(-1);
	}
	OutFile_Diag<<std::setiosflags(ios::scientific)<<std::setprecision(8);
	//
	writeMidXHead(OutFile_Diag,"Rho","Rho_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		int j = Nxp1 - k;
		Cell_2D &cell = *CarCellArray[k][j];
		OutFile_Diag<<cell.xc<<fs<<cell.MsQ().Rho<<endl;
	}
	//
	writeMidXHead(OutFile_Diag,"DP","DP_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		int j = Nxp1 - k;
		Cell_2D &cell = *CarCellArray[k][j];
		OutFile_Diag<<cell.xc<<fs<<cell.MsQ().p<<endl;
	}
	//
	writeMidXHead(OutFile_Diag,"TP","TP_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		int j = Nxp1 - k;
		Cell_2D &cell = *CarCellArray[k][j];
		OutFile_Diag<<cell.xc<<fs<<thermoPressure(cell)<<endl;
	}
	//
	writeMidXHead(OutFile_Diag,"U","U_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		int j = Nxp1 - k;
		Cell_2D &cell = *CarCellArray[k][j];
		OutFile_Diag<<cell.xc<<fs<<cell.MsQ().U<<endl;
	}
	//
	writeMidXHead(OutFile_Diag,"V","V_M",M_Phi);
	for(int k = 1;k < Nxp1;++k)
	{
		int j = Nxp1 - k;
		Cell_2D &cell = *CarCellArray[k][j];
		OutFile_Diag<<cell.xc<<fs<<cell.MsQ().V<<endl;
	}
	//
}
// void writeHead(ofstream &OutFile,int const subZone)
// {
// 	OutFile<<"(300 ("<<subZone<<" 2 1 0 0 1 "<<Cells<<")("<<endl;
// }
void writeHead(ofstream &OutFile,int const subZone,int const nD = 1)
{
	OutFile<<"(300 ("<<subZone<<" 2 "<<nD<<" 0 0 1 "<<Cells<<")("<<endl;
}
void Output_Flowfield(double const &t,int step)
{
		ostringstream oss_FlowField;
		oss_FlowField <<"../FlowField/global/" << "step" << step <<"Ma"<< Ma<<"_"
						<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"_T"<<t<<".dat";
		ofstream OutFile_FlowField;
		FileOpen(OutFile_FlowField,oss_FlowField,"OutFile_FlowField");

//		int subZone = 700;
		// !-------------------Rho-----------------------
		writeHead(OutFile_FlowField,101);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].MsQ().Rho<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
		//!-------------------U & V------------------------
		writeHead(OutFile_FlowField,2,2);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].MsQ().U<<"  "<<CellArray[n].MsQ().V<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
		//-------------------p---------------------------
		writeHead(OutFile_FlowField,1);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].MsQ().p<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
		//-------------------------T---------------------
		writeHead(OutFile_FlowField,3);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].MsQ().T<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
		//-------------------------Phi---------------------
		writeHead(OutFile_FlowField,700);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].MsQ().Phi<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
		//----------------tau-------------------
		#ifdef _ARK_MOMENTUM_FLIP
		writeHead(OutFile_FlowField,701);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].f.tau<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
		#endif
		//----------------Fx-------------------
		writeHead(OutFile_FlowField,702);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].MsQ().Fx<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
		//----------------Fy-------------------
		writeHead(OutFile_FlowField,703);
		LoopPS(Cells)
		{
			OutFile_FlowField<<CellArray[n].MsQ().Fy<<endl;
		}
		OutFile_FlowField<<"))"<<endl;
}
/* q<sub>x</sub>,q<sub>y</sub>,\
 <Greek>t</Greek><sub>xx</sub>,<Greek>t</Greek><sub>xy</sub>,\
 <Greek>t</Greek><sub>yy</sub>,
*/
//
// void Output_Flowfield(double const &t,int step)
// {
// 	ostringstream oss_FlowField;
// 	oss_FlowField <<"../FlowField/global/" << "step" << step <<"Ma"<< Ma<<"_"
// 					<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"_T"<<t<<".dat";
// 	ofstream OutFile_FlowField(oss_FlowField.str().c_str());
// 	if(!OutFile_FlowField)
// 	{
// 		_PRINT_ERROR_MSG_FLIP
// 		cout <<"  "<<"OutFile_FlowField open failed" << endl; 
// 		getchar();
// 		return;
// 	}
// 	ostringstream VarName,VarLocation,ZoneName,dataNE;
/*
// 	VarName << "VARIABLES = X,Y,<Greek>r</Greek>,U,V,p,T,\
// 	<Greek>f</Greek>,\
// 	<Greek>f</Greek><sub>x</sub>,<Greek>f</Greek><sub>y</sub>,\
// 	Fx,Fy\n";
*/
// 	VarLocation <<"VarLocation=([1-2]=NODAL,[3-12]=CellCentered)\n";
// 	ZoneName<<"ZONE T = Time" << t <<"_"<<"Mu"<<Mu0<<"\n";
// 	dataNE<<"Nodes="<<Nodes<<", Elements="<<Cells<<", ZONETYPE=FEQuadrilateral\n";
// 	string tecformat[5]={VarName.str().c_str(),
// 						ZoneName.str().c_str(),
// 						dataNE.str().c_str(),
// 						"DATAPACKING=BLOCK\n",
// 						VarLocation.str().c_str()};
// 	OutFile_FlowField << tecformat[0]<<tecformat[1]<<tecformat[2]<<tecformat[3]<<tecformat[4];
// 	OutFile_FlowField << setiosflags(ios::scientific) << setprecision(12);
// //	
// 	/*double *u_A = new double[Cells];
// 	double *v_A = new double[Cells];
// 	double *p_A = new double[Cells];
// 	for(int i = 0;i != Cells;++i)
// 		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A[i], v_A[i], p_A[i]);
// 	*/
// //---------------------------------------------Node->X-------------------------------------------
// 	for(int i = 0;i != Nodes;++i)
// 	{
// 		OutFile_FlowField << NodeArray[i].xN <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	OutFile_FlowField << endl;
// //--------------------------------------------Node->Y-------------------------------------------
// 	for(int i = 0;i != Nodes;++i)
// 	{
// 		OutFile_FlowField << NodeArray[i].yN <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	OutFile_FlowField << endl;
// //--------------------------------------------rho-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().Rho <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// //--------------------------------------------u-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().U <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// //--------------------------------------------v-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().V <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// //--------------------------------------------p-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().p <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// //--------------------------------------------T-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().T <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// /*
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().qx <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().qy <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].shearTau[0][0] <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].shearTau[0][1] <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].shearTau[1][1] <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// */
// 	// for(int i = 0;i != Cells;++i)
// 	// {
// 	// 	OutFile_FlowField << (*CellArray[i].pseudopsi) <<"   ";
// 	// 	if((i+1)%16 == 0)
// 	// 		OutFile_FlowField << "\n";
// 	// }
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().Phi <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().Phi_x <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().Phi_y <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().Fx <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << CellArray[i].MsQ().Fy <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// /*//--------------------------------------------u_A-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << u_A[i] <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// //--------------------------------------------v_A-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << v_A[i] <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}
// //--------------------------------------------p_A-------------------------------------------
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << p_A[i] <<"   ";
// 		if((i+1)%16 == 0)
// 			OutFile_FlowField << "\n";
// 	}*/
// //--------------------------------------------relation----------------------------------------	
// 	for(int i = 0;i != Cells;++i)
// 	{
// 		OutFile_FlowField << MeshIndex(CellArray[i].cellNodes[0] , NodeArray)<< " "
// 						  << MeshIndex(CellArray[i].cellNodes[1] , NodeArray)<< " "
// 					 	  << MeshIndex(CellArray[i].cellNodes[2] , NodeArray)<< " " 
// 						  << MeshIndex(CellArray[i].cellNodes[3] , NodeArray)<< endl;
// 	}
// 	OutFile_FlowField.close();
// 	/*delete []u_A;
// 	delete []v_A;
// 	delete []p_A;*/
// }
void Output_xcyc()
{
	ostringstream oss_xcyc;
	oss_xcyc <<"../FlowField/step" <<step<< "Nu"<< Nu0
					<<"_MeshCar"<<NL<<"_CFL"<<CFL<<"_xcyc.dat";
	ofstream OutFile_xcyc(oss_xcyc.str().c_str());
	if(!OutFile_xcyc)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_xcyc open failed" << endl; 
		getchar();
		return;
	}
	ostringstream title,varName,dataIJK,VarLocation;
	title << "title = \"Step"<<step<<"_"<<"Nu"<<Nu0<<"\""<<"\n";
	varName << "variables = x,y,<Greek>r</Greek>,u,v,p,Fx,Fy\n";
	dataIJK << "zone I = "<<Nx<<", J = "<<Ny<<", F = BLOCK\n";
	string tec360Header[3]={title.str().c_str(),varName.str().c_str(),
							dataIJK.str().c_str()};
//
	OutFile_xcyc<<tec360Header[0]<<tec360Header[1]<<tec360Header[2];
	OutFile_xcyc << setiosflags(ios::scientific) << setprecision(Out_precision);
	// for(int i = 0;i != Cells;++i)
	// {
	// 	OutFile_xcyc << CellArray[i].xc<<"  "<<CellArray[i].yc<<"\n";
	// }
	// OutFile_xcyc <<"--------------P_Inlet-----------------"<<'\n';
	// for(int n = 0;n != P_InletFaceNum;++n)
	// 	OutFile_xcyc << P_InletShadowCA[n].xc<<"  "<<P_InletShadowCA[n].yc<<"\n";
	// OutFile_xcyc <<"--------------P_Outlet----------------"<<'\n';
	// for(int n = 0;n != P_OutletFaceNum;++n)
	// 	OutFile_xcyc << P_OutletShadowCA[n].xc<<"  "<<P_OutletShadowCA[n].yc<<"\n";
	// for(int n = 0;n != PeriodicFaceNum;++n)
	// 	OutFile_xcyc << PeriodicShadowCA[n].xc<<"  "<<PeriodicShadowCA[n].yc<<"\n";
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->xc<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//	
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->yc<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->Rho<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->U<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->V<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->p<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	// for(int i = 1;i < Nxp1;++i)
	// for(int j = 1;j < Nyp1;++j)
	// {
	// 	OutFile_xcyc << *(CarCellArray[i][j]->pseudopsi)<<"  ";
	// 	if(j%10 == 0)
	// 		OutFile_xcyc<<"\n";
	// }
	// OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << (CarCellArray[i][j]->msq->Fx)<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << (CarCellArray[i][j]->msq->Fy)<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	OutFile_xcyc.close();
}