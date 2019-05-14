#ifndef _ARK_D2Q9_H
#define _ARK_D2Q9_H
//
#include <cmath>

int const Q = 9;

double const

MaxU = sqrt(3.0*RT),

Eta = 0,

MaSpan = sqrt(3.0*RT);

// namespace D2Q9{

// double const xi_u[Q] = {0,1,1,0,-1,-1,-1,0,1};

// double const xi_v[Q] = {0,0,1,1,1,0,-1,-1,-1};

// int const _BB[Q] = {0,5,6,7,8,1,2,3,4};

// }

namespace D2Q9{

double const 

A = 1.0/4.0,
B = 1.0/6.0,
C = 1.0/9.0,
D = 1.0/12.0,
E = 1.0/18.0,
F = 1.0/36.0;

double const omega[Q]={4.0/9.0,C,C,C,C,F,F,F,F};

double const xi_u[Q] = {0,1,0,-1,0,1,-1,-1,1};

double const xi_v[Q] = {0,0,1,0,-1,1,1,-1,-1};

int const _BB[Q] = {0,3,4,1,2,7,8,5,6};

double const M[Q][Q] = 
{
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1},
	{-4,-1,-1,-1,-1, 2, 2, 2, 2},
	{ 4,-2,-2,-2,-2, 1, 1, 1, 1},
	{ 0, 1, 0,-1, 0, 1,-1,-1, 1},
	{ 0,-2, 0, 2, 0, 1,-1,-1, 1},
	{ 0, 0, 1, 0,-1, 1, 1,-1,-1},
	{ 0, 0,-2, 0, 2, 1, 1,-1,-1},
	{ 0, 1,-1, 1,-1, 0, 0, 0, 0},
	{ 0, 0, 0, 0, 0, 1,-1, 1,-1}
};

double const IM[Q][Q] = 
{
	{ C,-C, C, 0, 0, 0, 0, 0, 0},
	{ C,-F,-E, B,-B, 0, 0, A, 0},
	{ C,-F,-E, 0, 0, B,-B,-A, 0},
	{ C,-F,-E,-B, B, 0, 0, A, 0},
	{ C,-F,-E, 0, 0,-B, B,-A, 0},
	{ C, E, F, B, D, B, D, 0, A},
	{ C, E, F,-B,-D, B, D, 0,-A},
	{ C, E, F,-B,-D,-B,-D, 0, A},
	{ C, E, F, B, D,-B,-D, 0,-A}
};

}

#endif
// const double RT = 1.0/3.0,s_RT = sqrt(3.0*RT);
// //
// const typexi xi[Q]={{0,0},{s_RT,0},{s_RT,s_RT},{0,s_RT},
// 					{-s_RT,s_RT},{-s_RT,0},{-s_RT,-s_RT},{0,-s_RT},{s_RT,-s_RT}};
// //
// const double omega[Q]={4.0/9.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0};					
// #endif