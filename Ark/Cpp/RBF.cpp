#include "DUGKSDeclaration.h"
#include <iostream>
#include <iomanip>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using std::cout;
using std::endl;
using std::ios;
using std::setprecision;
using std::setiosflags;

typedef Matrix<double,9,9> Matrix9d;

extern void Output_Flowfield(double const &t,int step);

namespace RBF{

	inline double kernel(double const Euclid){return sqrt(Euclid/9.0 + RBF::c);}

	const double 
	//if dx == 3 dy == 4
	A = kernel(::dx*::dx), //3
	B = kernel(::dy*::dy), //4
	C = kernel(::dx*::dx + ::dy*::dy), //5
	D = kernel(0.0),
	E = kernel(4 * ::dx*::dx), //6.00083
	F = kernel(4 * ::dy*::dy), //8.00062
	G = kernel(4 * ::dx*::dx + ::dy*::dy), //7.2118
	H = kernel(::dx*::dx + 4 * ::dy*::dy), //8.54459
	I = kernel(4 * ::dx*::dx + 4 * ::dy*::dy); //10.0005

	const double OriginalM[DIM][DIM] =
	{
		{D,A,B,A,B,C,C,C,C},
		{A,D,C,E,C,B,G,G,B},
		{B,C,D,C,F,A,A,H,H},
		{A,E,C,D,C,G,B,B,G},
		{B,C,F,C,D,H,H,A,A},
		{C,B,A,G,H,D,E,I,F},
		{C,G,A,B,H,E,D,F,I},
		{C,G,H,B,A,I,F,D,E},
		{C,B,H,G,A,F,I,E,D}
	};

	double InverseM[DIM][DIM] = {0};
}

void RBF::setInverseM()
{
	Matrix9d MTmp,inverseMTmp;

	for(int i = 0;i < RBF::DIM;++i)
	for(int j = 0;j < RBF::DIM;++j)
	{
		MTmp(i,j) = RBF::OriginalM[i][j];
	}

	inverseMTmp = MTmp.inverse();

	for(int i = 0;i < RBF::DIM;++i)
	for(int j = 0;j < RBF::DIM;++j)
	{
		InverseM[i][j] = inverseMTmp(i,j);
	}
}
//x := deltax, y := deltay

extern Cell_2D* targetCell(Cell_2D* cellptr);

double valueRBF(double const *wRBF, Cell_2D* cellptr,double const x,double const y)
{
	using RBF::kernel;

	cellptr = targetCell(cellptr);
	double 
	tmp = wRBF[0]*kernel(x*x + y*y)
	    + wRBF[1]*kernel((::dx-x)*(::dx-x) + y*y)
	    + wRBF[2]*kernel(x*x + (::dy-y)*(::dy-y))
	    + wRBF[3]*kernel((::dx+x)*(::dx+x) + y*y)
	    + wRBF[4]*kernel(x*x + (::dy+y)*(::dy+y))
	    + wRBF[5]*kernel((::dx-x)*(::dx-x) + (::dy-y)*(::dy-y))
	    + wRBF[6]*kernel((::dx+x)*(::dx+x) + (::dy-y)*(::dy-y))
	    + wRBF[7]*kernel((::dx+x)*(::dx+x) + (::dy+y)*(::dy+y))
	    + wRBF[8]*kernel((::dx-x)*(::dx-x) + (::dy+y)*(::dy+y));
	return tmp;
}

void SetwRBF(double *wRBF,Cell_2D *cellptr,Cell_2D::DVDF Cell_2D::*dvdf,int const k)
{
	using RBF::InverseM;

	cellptr = targetCell(cellptr);
	for(int m = 0;m < RBF::DIM;++m)
	{
		wRBF[m] = 
			InverseM[m][0]*(cellptr->*dvdf).BarP[k]

		+	InverseM[m][1]*(cellptr->Cell_C[0]->*dvdf).BarP[k]
		+	InverseM[m][2]*(cellptr->Cell_C[1]->*dvdf).BarP[k]
		+	InverseM[m][3]*(cellptr->Cell_C[2]->*dvdf).BarP[k]
		+	InverseM[m][4]*(cellptr->Cell_C[3]->*dvdf).BarP[k]

		+	InverseM[m][5]*(cellptr->Cell_Diag[0]->*dvdf).BarP[k]
		+	InverseM[m][6]*(cellptr->Cell_Diag[1]->*dvdf).BarP[k]
		+	InverseM[m][7]*(cellptr->Cell_Diag[2]->*dvdf).BarP[k]
		+	InverseM[m][8]*(cellptr->Cell_Diag[3]->*dvdf).BarP[k];
	}
}
void SetwRBF(double *wRBF, Cell_2D *cellptr,double MacroQuantity::*var)
{
	using RBF::InverseM;

	cellptr = targetCell(cellptr);
	for(int m = 0;m < RBF::DIM;++m)
	{
		wRBF[m] = 
			InverseM[m][0]*(cellptr->MsQ()).*var

		+	InverseM[m][1]*(cellptr->Cell_C[0]->MsQ()).*var
		+	InverseM[m][2]*(cellptr->Cell_C[1]->MsQ()).*var
		+	InverseM[m][3]*(cellptr->Cell_C[2]->MsQ()).*var
		+	InverseM[m][4]*(cellptr->Cell_C[3]->MsQ()).*var

		+	InverseM[m][5]*(cellptr->Cell_Diag[0]->MsQ()).*var
		+	InverseM[m][6]*(cellptr->Cell_Diag[1]->MsQ()).*var
		+	InverseM[m][7]*(cellptr->Cell_Diag[2]->MsQ()).*var
		+	InverseM[m][8]*(cellptr->Cell_Diag[3]->MsQ()).*var;
	}
}
extern double RBFTestAnalytic(double const x,double const y);

void RBFPrecisionTest()
{
	// LoopPS(Cells)
	// {
	// 	Cell_2D &cell = CellArray[n];
	// 	SetwRBF(&cell,&Cell_2D::h,5);
	// }
	// LoopPS(Faces)
	// {
	// 	Face_2D &face = FaceArray[n];
	// 	double 
	// 	dxOwner = face.xf - face.owner->xc,
	// 	dyOwner = face.yf - face.owner->yc,
	// 	dxNeigh = face.xf - face.neigh->xc,
	// 	dyNeigh = face.yf - face.neigh->yc;

	// 	face.h.hDt[5] = 0.0;

	// 	face.h.hDt[5] += valueRBF(face.owner,dxOwner,dyOwner);
	// 	face.h.hDt[5] += valueRBF(face.neigh,dxNeigh,dyNeigh);
	// 	face.h.hDt[5] /= 2;
	// }
	// LoopPS(Faces)
	// {
	// 	Face_2D &face = FaceArray[n];
	// 	Info << (RBFTestAnalytic(face.xf,face.yf))<<fs;
	// 	cout << (face.h.hDt[5])<<fs;
	// 	cout << face.xf << fs << face.yf<<endl;
	// }
	// LoopPS(Cells)
	// {
	// 	Cell_2D &cell = CellArray[n];
	// 	double tmpAN = cell.Cell_C[3]->MsQ().Rho - valueRBF(&cell,0,-::dy);
	// 	if(fabs(tmpAN) > ::infinitesimal)
	// 	Info << tmpAN<<endl;
	// }
	// double ErrorLmax = 0.0,Errortmp = 0.0;
	// double ErrorL2 = 0.0, dVar = 0.0,sumVar = 0.0,sumdVar = 0.0;
	// LoopPS(Faces)
	// {
	// 	Face_2D &face = FaceArray[n];

	// 	Errortmp = fabs(face.h.hDt[5] - RBFTestAnalytic(face.xf,face.yf));
	// 	if(Errortmp > ErrorLmax)
	// 		ErrorLmax = Errortmp;

	// 	dVar = fabs(face.h.hDt[5] - RBFTestAnalytic(face.xf,face.yf));
	// 	sumdVar += dVar*dVar;
	// 	sumVar += RBFTestAnalytic(face.xf,face.yf)*RBFTestAnalytic(face.xf,face.yf);
	// }
	// ErrorL2 = sqrt(sumdVar/(sumVar+1e-30));
	// Info << "ErrorL2 : "<<ErrorL2 << endl;
	// Info << "ErrorLmax : "<<ErrorLmax << endl;
	// Output_Flowfield((step)*::dt,step);
}