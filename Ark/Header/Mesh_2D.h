#ifndef _MESH_2D_ARK_H_
#define _MESH_2D_ARK_H_

#include <vector>
#include <cmath>
#include "ZeroNamespace.h"
#include "Macroscopic.h"
using std::vector;

class Node_2D;

class Face_2D;

class Cell_2D;

// enum MME {PHI = 0,MOM,ENE,Num_t};

// enum Bound {ENPTY = 0,RIGHT,LEFT,TOP,BOTTOM};


#include "Node_2D.h"
#include "Face_2D.h"
#include "Cell_2D.h"


inline 
void SetZero(double &d)
{
	if(d < infinitesimal && d > -infinitesimal) d= 0.0;
}
inline 
bool EqualZero(double const &d)
{
	return (d < infinitesimal && d > -infinitesimal) ? true : false;
}
inline 
bool doubleEqual(double const &c,double const &d)
{
	return EqualZero(c-d);
}

extern void allocateARK(double* &f,int const Q);

extern void deallocateARK(double* &f,int const Q);

extern double TriArea(double , double , double , double ,double , double );

#endif