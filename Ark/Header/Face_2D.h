#ifndef _FACE2D_ARK_H
#define _FACE2D_ARK_H

class Face_2D{
//
public:
//-----------constructors-------------------
	Face_2D();
	Face_2D(const Face_2D &rhs);
	Face_2D& operator=(const Face_2D &rhs);
	~Face_2D();
//
	Face_2D *shadowF = nullptr;
//
	int bc_type = 2;
	unsigned zone = 0;
	Node_2D *faceNodes[2] = {nullptr,nullptr};
	Face_2D *faceFaces[9] = {nullptr};
	Cell_2D *faceCells[4] = {nullptr};
	Cell_2D *owner = nullptr,*neigh = nullptr;
	double xf = 0.0, yf = 0.0, Area = 0.0;
	double Vx = 0.0,Vy = 0.0;
	double _dx = 0.0, _dy = 0.0;
	void SetArea();
	void SetNormalV();
//----------------------------------------------
	double *xi_n_dS = nullptr;
//
	class DVDF	//discrete velocity distribution function
	{
	public:
		DVDF();
		DVDF(const DVDF &rhs);
		DVDF& operator=(const DVDF &rhs);
		~DVDF();
//
		double tauh = 0.0;
		double ah = 0.0,bh = 0.0,ch = 0.0;
		inline void setxh();
//
		double *hDt     = nullptr;	//half time step
		double *BhDt    = nullptr;
		double *EqhDt   = nullptr;
		double *SohDt	= nullptr;	//source term
	private:
		int *token = nullptr;
		void assign(const DVDF& rhs);
		void allocate();
		void deallocate();
	};
//-----------------------discrete velocity distribution function-----------------
	#ifdef _ARK_ALLENCAHN_FLIP
	DVDF h;
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	DVDF f;
	#endif
//
	#ifdef _ARK_THERMAL_FLIP
	DVDF g;
	#endif
//-----------macro variables---------------------
	MacroQuantity * __restrict__ msqh;
	inline MacroQuantity& MsQh(){return *msqh;}
	inline MacroQuantity const & MsQh() const {return *msqh;}
//-----------------Factor------------------------
	void Factor();
	void FactorAC();
//
private:
	int *use = nullptr;
	void assign(const Face_2D& rhs);
	void allocate();
	void deallocate();
};

#endif