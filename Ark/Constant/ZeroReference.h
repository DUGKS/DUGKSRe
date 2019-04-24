#ifndef _ZERO_REFERENCE_H
#define _ZERO_REFERENCE_H

#include <cmath>
#include "ZeroFlip.h"

//----------------------Mesh file---------------
const int 

Nx = 128,

Ny = 128,

NL = 128;//name of mesh file

const double

ChLength = 128.0,

MinL = ChLength/NL,//ChLength/NL,//6.037849742228585e-02,//3.079505855617000e-02,//1.555181192035053e-02,
//7.776141069016656e-03,//
X_Beg = 0.0,

X_End = 128,

Y_Beg = 0.0,

Y_End = 128,

Lx = X_End - X_Beg,

Ly = Y_End - Y_Beg,

dx = Lx/Nx,

dy = Ly/Ny,

Kh3 = 1.0E-3*MinL*MinL*MinL; //Limiter


int const 

Nxp1 = Nx + 1, 

Nxp2 = Nx + 2,

Nyp1 = Ny + 1, 

Nyp2 = Ny + 2;

//------------------------------Atomic species-------------------------------------
const double 

Omega0 = 0.5, //hard sphere(HS) = 0.5, variable hard sphere(VHS) = 0.68

Pr = 2.0/3.0, //Prandtl Number

nK = 0,//internal degree 0 = single  2 = double

agEq = (nK + 3.0 - 2.0)/2.0,//used in the Equilibrium of g(x,xi,t)

Cv = (nK + 3.0)/2.0,

Gamma = (nK + 5.0)/(nK + 3.0);

//-----------------------------reference parameters-------------------------

const double 

T0 = 1.0/3,

R0 = 1.0,

RT = R0*T0,

Lambda0 = 1/(2.0*R0*T0),

Rho0 = 1.0,//4.5435E-2,//CS

U0 = 0.0,//W_i*TC_r,

V0 = 0.0,

Ma = U0/sqrt(R0*T0),

p0 = 1,

Re = 10.0,

// Mu0 = U0*ChLength*Rho0/Re,

Mu0 = 0.3,//3.413333E-3,//,//1.0/6.0,//

Nu0 = Mu0/Rho0,//3.413333E-3,//0.0682667,//3.413333E-3,//0.0682667,

// Mu0 = Rho0*U0*ChLength/Re,

// Nu0 = Mu0/Rho0,

Tau0 = Nu0/RT,

Kn = 16.0*Mu0/(5.0*Rho0*R0*T0)*sqrt(1.0/(4.0*Lambda0*PI));

namespace ARK{

int const 

ND = 2,

cellN = 4,

cellF = 4,

faceN = 2,

faceF = 0;

double const

digits = 1e10;

}

#endif
