#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include "DUGKSDeclaration.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::ios;
using std::setprecision;
using std::setiosflags;

int const Out_precision = 12;

int const IC = 198, JC = 197;
//
void printSplitLine(char c)
{
	cout <<"----------------------------------"<<c<<endl;
}
void printErrorLine(char c)
{
	cout <<"||***?????????????????????????????***||"<<c<<endl;
}
void printErrorMessage(int line,const char *file,const char *func)
{
	cout<<"File : "<<file<<"  Line : "<<line<<"  fun : "<<func<<'\n';
}
#ifdef _ZERO_NDEBUG_FLIP
void printCorners(Cell_2D *cellptr)
{
	cout <<*cellptr->use<<"    "<<cellptr->xc<<"    "<<cellptr->yc<<endl;
}
void Output_UVP(double const &t)
{
	ostringstream oss_rho,oss_u,oss_v,oss_p,oss_uA,oss_vA,oss_pA;
	ofstream OutFile_rho,OutFile_u, OutFile_v, OutFile_p, OutFile_uA, OutFile_vA, OutFile_pA;
//--------------------------------------------------------------------------------
	double *u_A = new double[Cells];
	double *v_A = new double[Cells];
	double *p_A = new double[Cells];
	for(int i = 0;i != Cells;++i)
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A[i], v_A[i], p_A[i]);
//---------------------------------------------------------------------------------
	oss_XXX(oss_rho,"UVP","rho",t);
	oss_XXX(oss_u,"UVP","u",t);
	oss_XXX(oss_v,"UVP","v",t);
	oss_XXX(oss_p,"UVP","p",t);
	oss_XXX(oss_uA,"UVP","uA",t);
	oss_XXX(oss_vA,"UVP","vA",t);
	oss_XXX(oss_pA,"UVP","pA",t);
	FileOpen(OutFile_rho,oss_rho,"rho");
	FileOpen(OutFile_u,oss_u,"u");
	FileOpen(OutFile_v,oss_v,"v");
	FileOpen(OutFile_p,oss_p,"p");
	FileOpen(OutFile_uA,oss_uA,"uA");
	FileOpen(OutFile_vA,oss_vA,"vA");
	FileOpen(OutFile_pA,oss_pA,"pA");
//--------------------------------------------rho-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_rho << CellArray[i].MsQ().Rho <<"\n";
	}
//--------------------------------------------u-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_u << CellArray[i].MsQ().U <<"\n";
	}
//--------------------------------------------v-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_v << CellArray[i].MsQ().V <<"\n";
	}
//--------------------------------------------p-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_p << CellArray[i].MsQ().p <<"\n";
	}
//--------------------------------------------u_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_uA << u_A[i] <<"\n";
	}
//--------------------------------------------v_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_vA << v_A[i] <<"\n";
	}
//--------------------------------------------p_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_pA << p_A[i] <<"\n";
	}
	OutFile_rho.close();
	OutFile_u.close();
	OutFile_v.close();
	OutFile_p.close();
	OutFile_uA.close();
	OutFile_vA.close();
	OutFile_pA.close();
	delete []u_A;
	delete []v_A;
	delete []p_A;
}
void Output_fBP(double const &t,int ii,int jj)
{
	ostringstream oss_fBP;
	oss_fBP <<"../FlowField/fBP/" << "Time" << t <<"_Mu"<< Mu0
					<<"_MeshCar"<<NL<<"-"<<NL<<"_CFL"<<CFL<<"_fBP"<<ii<<"_"<<jj<<".dat";
	ofstream OutFile_fBP(oss_fBP.str().c_str());
	if(!OutFile_fBP)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_fBP open failed" << endl;
		getchar();
		return;
	}
	OutFile_fBP << setiosflags(ios::scientific) << setprecision(12);
	for(int n = 0;n != Cells;++n)
	{
		OutFile_fBP << CellArray[n].f.BarP[ii][jj] <<"\n";
	}
	OutFile_fBP.close();
}
void Output_fBh(Face_2D& face,double t)
{
	ostringstream oss_fBh;
	ofstream OutFile_fBh;
	oss_XXX(oss_fBh,"fBh","fBh",t);
	FileOpen(OutFile_fBh,oss_fBh,"fBh");
	//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_fBh<<face.f.BhDt[i][j]<<'\n';
		// if(face.xi_n_dS[i][j] >= 0)
		// {
		// 	OutFile_fBh <<face.owner->f.BarP[i][j]<<"    "
		// 				<<face.owner->f.Eq[i][j]<<"    "
		// 				<<face.f.BhDt[i][j] - face.owner->f.BarP[i][j]<<"    "
		// 				<<face.owner->f.Eq[i][j] - face.owner->f.BarP[i][j]<<endl;
		// }
		// else
		// {
		// 	OutFile_fBh <<face.neigh->f.BarP[i][j]<<"    "
		// 				<<face.neigh->f.Eq[i][j]<<"    "
		// 				<<face.f.BhDt[i][j] - face.neigh->f.BarP[i][j]<<"    "
		// 				<<face.neigh->f.Eq[i][j] - face.neigh->f.BarP[i][j]<<endl;
		// }
	}
	OutFile_fBh.close();
}
void Output_fh(Face_2D& face,double t)
{
	ostringstream oss_fh;
	ofstream OutFile_fh;
	oss_XXX(oss_fh,"fh","fh",t);
	FileOpen(OutFile_fh,oss_fh,"fh");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_fh << face.f.hDt[i][j] <<"    "
					   << face.f.BhDt[i][j]<<"    "
					   << face.f.EqhDt[i][j]<<'\n';
		}
	OutFile_fh <<"Rho_h : "<<face.MsQh().Rho<<'\n'
			   <<"U_h : "<<face.MsQh().U<<'\n'
			   <<"V_h : "<<face.MsQh().V<<'\n'
			   <<"T_h: "<<face.MsQh().T<<'\n'
			   <<"qx_h : "<<face.MsQh().qx<<'\n'
			   <<"qy_h : "<<face.MsQh().qy<<'\n'
			   <<"xc : "<<face.xf<<'\n'
			   <<"yc : "<<face.yf<<'\n';
	OutFile_fh.close();
}
void Output_fh_Append(Face_2D& face,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_fh;
	ofstream OutFile_fh;
	oss_XXX(oss_fh,"fh","fhAPP",dt);
	FileOpenAppend(OutFile_fh,oss_fh,"fh");
	OutFile_fh 	<< face.f.hDt[I][J] <<"    "<<face.f.ah<<"    "<<face.f.bh<<'\n';
				//<< face.f.BhDt[i][j]<<"    "
				//<< face.f.EqhDt[i][j]<<'\n';
	OutFile_fh.close();				
}
void Output_fT(Cell_2D &cell,double t)
{
	ostringstream oss_fT;
	ofstream OutFile_fT;
	oss_XXX(oss_fT,"fT","fT",t);
	FileOpen(OutFile_fT,oss_fT,"fT");
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_fT <<cell.f.Tilde[i][j]<<"    "
				   <<cell.f.Eq[i][j]<<"    "
//				   <<cell.fFlux[i][j]<<"    "
				   <<cell.f.Tilde[i][j]-cell.f.Eq[i][j]
				   <<'\n';
	}
	OutFile_fT <<"Rho : "<<cell.MsQ().Rho<<'\n'
				<<"U : "<<cell.MsQ().U<<'\n'
				<<"V : "<<cell.MsQ().V<<'\n'
				<<"T: "<<cell.MsQ().T<<'\n'
				<<"Lambda: "<<cell.MsQ().Lambda<<'\n'
				<<"qx : "<<cell.MsQ().qx<<'\n'
				<<"qy : "<<cell.MsQ().qy<<'\n'
				<<"xc : "<<cell.xc<<'\n'
				<<"yc : "<<cell.yc<<'\n';
	OutFile_fT.close(); 
}
void Output_fT_Append(Cell_2D &cell,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_fT;
	ofstream OutFile_fT;
	oss_XXX(oss_fT,"fT","fTAPP",dt);
	FileOpenAppend(OutFile_fT,oss_fT,"fT");
	OutFile_fT  <<cell.f.Tilde[I][J]<<"    "//<<cell.aBP<<"    "<<cell.bBP
				<<cell.Cell_C[0]->f.Tilde[I][J]<<"    "
				<<cell.Cell_C[1]->f.Tilde[I][J]<<"    "
				<<cell.Cell_C[2]->f.Tilde[I][J]<<"    "
				<<cell.Cell_C[3]->f.Tilde[I][J]<<"    "
				<<'\n';
	OutFile_fT.close();
}

// void Output_fFlux_Append(Cell_2D &cell,double dt)
// {
// 	int I = IC,J = JC;
// 	ostringstream oss_fFlux;
// 	ofstream OutFile_fFlux;
// 	oss_XXX(oss_fFlux,"fFlux","fFluxAPP",dt);
// 	FileOpenAppend(OutFile_fFlux,oss_fFlux,"fFlux");
// 	OutFile_fFlux  <<cell.fFlux[I][J]<<"    "
// 				    <<cell.Face_C[0]->fh[I][J]<<"    "
// 				    <<cell.Face_C[1]->fh[I][J]<<"    "
// 				    <<cell.Face_C[2]->fh[I][J]<<"    "
// 				    <<cell.Face_C[3]->fh[I][J]<<"    "
// 				    <<'\n';
// 	OutFile_fFlux.close();
// }
// void Output_gFlux_Append(Cell_2D &cell,double dt)
// {
// 	int I = IC,J = JC;
// 	ostringstream oss_gFlux;
// 	ofstream OutFile_gFlux;
// 	oss_XXX(oss_gFlux,"gFlux","gFluxAPP",dt);
// 	FileOpenAppend(OutFile_gFlux,oss_gFlux,"gFlux");
// 	OutFile_gFlux  <<cell.gFlux[I][J]<<"    "
// 				    <<cell.Face_C[0]->gh[I][J]<<"    "
// 				    <<cell.Face_C[1]->gh[I][J]<<"    "
// 				    <<cell.Face_C[2]->gh[I][J]<<"    "
// 				    <<cell.Face_C[3]->gh[I][J]<<"    "
// 				    <<'\n';
// 	OutFile_gFlux.close();
// }
#ifndef _ARK_ISOTHERMAL_FLIP
void Output_gBh(Face_2D& face,double t)
{
	ostringstream oss_gBh;
	ofstream OutFile_gBh;
	oss_XXX(oss_gBh,"gBh","gBh",t);
	FileOpen(OutFile_gBh,oss_gBh,"gBh");
	//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_gBh<<face.g.BhDt[i][j]<<'\n';
		// if(face.xi_n_dS[i][j] >= 0)
		// {
		// 	OutFile_gBh <<face.owner->g.BarP[i][j]<<"    "
		// 				<<face.owner->gEq[i][j]<<"    "
		// 				<<face.g.BhDt[i][j] - face.owner->g.BarP[i][j]<<"    "
		// 				<<face.owner->gEq[i][j] - face.owner->g.BarP[i][j]<<endl;
		// }
		// else
		// {
		// 	OutFile_gBh <<face.neigh->g.BarP[i][j]<<"    "
		// 				<<face.neigh->gEq[i][j]<<"    "
		// 				<<face.g.BhDt[i][j] - face.neigh->g.BarP[i][j]<<"    "
		// 				<<face.neigh->gEq[i][j] - face.neigh->g.BarP[i][j]<<endl;
		// }
	}
	OutFile_gBh.close();
}
void Output_gT(Cell_2D &cell,double t)
{
	ostringstream oss_gT;
	ofstream OutFile_gT;
	oss_XXX(oss_gT,"gT","gT",t);
	FileOpen(OutFile_gT,oss_gT,"gT");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_gT <<cell.g.Tilde[i][j]<<"    "
//					   <<cell.f.Eq[i][j]<<"    "
//					   <<cell.fFlux[i][j]<<"    "
//					   <<cell.g.Tilde[i][j]-cell.f.Eq[i][j]-cell.DtSlashVolume*cell.fFlux[i][j]
					   <<endl;
		}
	OutFile_gT <<"Rho : "<<cell.Rho<<'\n'
				<<"U : "<<cell.U<<'\n'
				<<"V : "<<cell.V<<'\n'
				<<"T: "<<cell.T<<'\n'
				<<"Lambda: "<<cell.Lambda<<'\n'
				<<"qx : "<<cell.qx<<'\n'
				<<"qy : "<<cell.qy<<'\n'
				<<"xc : "<<cell.xc<<'\n'
				<<"yc : "<<cell.yc<<'\n';
	OutFile_gT.close(); 
}
void Output_gh_Append(Face_2D& face,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_gh;
	ofstream OutFile_gh;
	oss_XXX(oss_gh,"gh","ghAPP",dt);
	FileOpenAppend(OutFile_gh,oss_gh,"gh");
	OutFile_gh 	<< face.g.hdt[I][J] <<"    "<<'\n';
				//<< face.f.BhDt[i][j]<<"    "
				//<< face.f.EqhDt[i][j]<<'\n';
	OutFile_gh.close();				
}
void Output_gh(Face_2D& face,double t)
{
	ostringstream oss_gh;
	ofstream OutFile_gh;
	oss_XXX(oss_gh,"gh","gh",t);
	FileOpen(OutFile_gh,oss_gh,"gh");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_gh << face.g.hDt[i][j] <<"    "
					   << face.g.BhDt[i][j]<<"    "
					   << face.g.EqhDt[i][j]<<'\n';
		}
	OutFile_gh.close();
}
void Output_gT_Append(Cell_2D &cell,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_gT;
	ofstream OutFile_gT;
	oss_XXX(oss_gT,"gT","gTAPP",dt);
	FileOpenAppend(OutFile_gT,oss_gT,"gT");
	OutFile_gT  <<cell.g.Tilde[i][j]<<"    "
				// <<cell.g.BarP[i][j]<<"    "
				// <<cell.g.Eq[i][J]<<"    "
				<<cell.Cell_C[0]->g.Tilde[i][j]<<"    "
				<<cell.Cell_C[1]->g.Tilde[i][j]<<"    "
				<<cell.Cell_C[2]->g.Tilde[i][j]<<"    "
				<<cell.Cell_C[3]->g.Tilde[i][j]
				<<'\n';
	OutFile_gT.close();
}
#endif
void Output_phi_Bh(Face_2D &face,double t)
{
	ostringstream oss_fBh_L;
	ostringstream oss_fBh_R;
	ostringstream oss_gBh_L;
	ostringstream oss_gBh_R;
	ofstream OutFile_fBh_L;
	ofstream OutFile_fBh_R;
	ofstream OutFile_gBh_L;
	ofstream OutFile_gBh_R;
	oss_XXX(oss_fBh_L,"fBh","fBhL",t);
	oss_XXX(oss_fBh_R,"fBh","fBhR",t);
	oss_XXX(oss_gBh_L,"fBh","gBhL",t);
	oss_XXX(oss_gBh_R,"fBh","gBhR",t);
	FileOpen(OutFile_fBh_L,oss_fBh_L,"fBh_L");
	FileOpen(OutFile_fBh_R,oss_fBh_R,"fBh_R");
	FileOpen(OutFile_gBh_L,oss_gBh_L,"gBh_L");
	FileOpen(OutFile_gBh_R,oss_gBh_R,"gBh_R");
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		if(xi_u[QuIndex]>=0)
		{
			OutFile_fBh_R<<face.f.BhDt[i][j]<<'\n';
			#ifndef _ARK_ISOTHERMAL_FLIP
			OutFile_gBh_R<<face.g.BhDt[i][j]<<'\n';
			#endif
		}
		else
		{
			OutFile_fBh_L<<face.f.BhDt[i][j]<<'\n';
			#ifndef _ARK_ISOTHERMAL_FLIP
			OutFile_gBh_L<<face.g.BhDt[i][j]<<'\n';
			#endif
		}
	}
	OutFile_fBh_L.close();
	OutFile_fBh_R.close();
	OutFile_gBh_L.close();
	OutFile_gBh_R.close();
}
#endif