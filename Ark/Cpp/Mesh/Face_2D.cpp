#include "Mesh_2D.h"

void Face_2D::assign(const Face_2D& rhs)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	h = rhs.h;
	#endif
	//
	#ifdef _ARK_MOMENTUM_FLIP
	f = rhs.f;
	#endif
	//
	#ifndef _ARK_ISOTHERMAL_FLIP
	g = rhs.g;
	#endif
//
	xi_n_dS = rhs.xi_n_dS;
	msqh = rhs.msqh;
}
void Face_2D::allocate()
{
	allocateARK(xi_n_dS,Q);
}
void Face_2D::deallocate()
{
	deallocateARK(xi_n_dS,Q);
	delete msqh;
}
Face_2D::Face_2D():msqh(new MacroQuantity()),use(new int(1))
{
	allocate();
}
Face_2D::Face_2D(const Face_2D& rhs)
{
	assign(rhs);
//
	use = rhs.use;
	++*use;
}
Face_2D& Face_2D::operator=(const Face_2D& rhs)
{
	if(use != rhs.use)
	{
		if(--*use == 0)
		{	
			deallocate();
			delete use;
		}
		assign(rhs);
//
		use = rhs.use;
		++*use;
	}
	return *this;
}
Face_2D::~Face_2D()
{
//
	if(--*use == 0)	
	{
		deallocate();
//
		delete use;
	}
}
//--------------------------------Face::DVDF-----------------------------
void Face_2D::DVDF::allocate()
{
	allocateARK(hDt,Q);
	allocateARK(BhDt,Q);
	allocateARK(EqhDt,Q);
	allocateARK(SohDt,Q);
}
void Face_2D::DVDF::deallocate()
{
	deallocateARK(hDt,Q);
	deallocateARK(BhDt,Q);
	deallocateARK(EqhDt,Q);
	deallocateARK(SohDt,Q);
}
void Face_2D::DVDF::assign(const Face_2D::DVDF& rhs)
{
	hDt = rhs.hDt;
	BhDt = rhs.BhDt;
	EqhDt = rhs.EqhDt;
	SohDt = rhs.SohDt;
}
Face_2D::DVDF::DVDF():token(new int(1))
{
	allocate();
}
Face_2D::DVDF::DVDF(const DVDF& rhs)
{
	assign(rhs);

	token = rhs.token;
	++*token;
}
Face_2D::DVDF& Face_2D::DVDF::operator=(const DVDF& rhs)
{
	if(token != rhs.token)
	{
		if(--*token == 0)
		{
			deallocate();
		//
			delete token;
		}
		assign(rhs);
//
		token = rhs.token;
		++*token;
	}
	return *this;
}
Face_2D::DVDF::~DVDF()
{
	if(--*token == 0)
	{
		deallocate();
	}
}
void Face_2D::DVDF::setxh()
{
	ah = 2.0*tauh/(2.0*tauh+::hDt);
	bh = 1.0 - ah;
	ch = tauh*bh;
}
void Face_2D::Factor()
{
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	f.setxh();
	#endif
	//energy
	#ifndef _ARK_ISOTHERMAL_FLIP
	g.setxh();
	#endif
}
void Face_2D::FactorAC()
{
	#ifdef _ARK_ALLENCAHN_FLIP
	h.setxh();
	#endif
}
void Face_2D::SetArea()
{
	xf = 0.5*(faceNodes[0]->xN + faceNodes[1]->xN);
	yf = 0.5*(faceNodes[0]->yN + faceNodes[1]->yN);
	Area = sqrt(
				(faceNodes[1]->yN - faceNodes[0]->yN) * (faceNodes[1]->yN - faceNodes[0]->yN)
				 +
				(faceNodes[1]->xN - faceNodes[0]->xN) * (faceNodes[1]->xN - faceNodes[0]->xN)
			   );
}
void Face_2D::SetNormalV()
{
	double dy = (faceNodes[1]->yN - faceNodes[0]->yN), 
	 	   dx = (faceNodes[1]->xN - faceNodes[0]->xN);
	SetZero(dx);SetZero(dy); 	 	   
	Vx = dy/Area;Vy = -dx/Area;
}