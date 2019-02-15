#include "Mesh_2D.h"

double TriArea(double Xbeg, double Ybeg, double X_A, double Y_A,double X_B, double Y_B)
{
	double Area = 0.5 * ((X_A - Xbeg) * (Y_B - Ybeg) - (X_B - Xbeg) * (Y_A - Ybeg));
	return (Area > 0 ? Area : -Area);
}
void allocateARK(double* &f,int const Q)
{
	f = new double[Q];
}
void deallocateARK(double* &f,int const Q)
{
	delete[] f;
	f = nullptr;
}