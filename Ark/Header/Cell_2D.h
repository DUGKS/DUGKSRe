#ifndef _CELL2D_ARK_H
#define _CELL2D_ARK_H
//
class Cell_2D{
//
public:
	Cell_2D();
	Cell_2D(const Cell_2D &rhs);
	Cell_2D& operator=(const Cell_2D &rhs);
	~Cell_2D();
	static int const nN = 4;
//-----------Mesh------------------
	Cell_2D *ShadowC = nullptr; //always be nullptr for nonshadow cell 
	Face_2D *BoundF = nullptr;
	unsigned zone = 0;
	int celltype = 0;
	double xc = 0.0, yc = 0.0, volume = 0.0, DtSlashVolume = 0.0;
	// double *NodeX_C[4] = {nullptr,nullptr,nullptr,nullptr};
	// double *NodeY_C[4] = {nullptr,nullptr,nullptr,nullptr};
	Node_2D *cellNodes[nN] = {nullptr,nullptr,nullptr,nullptr};
//
	Face_2D *Face_C[nN] = {nullptr,nullptr,nullptr,nullptr};
//
	Cell_2D *Cell_C[nN] = {nullptr,nullptr,nullptr,nullptr};
	Cell_2D *Cell_Diag[nN] = {nullptr,nullptr,nullptr,nullptr};
//
	int signFlux[nN] = {0};
//-------------------------Least Square-------------------------------
	double LS_M[2][2] = {{0.0,0.0},{0.0,0.0}};
	double wdx_C[4] = {0.0};
	double wdy_C[4] = {0.0};
//---------------------Radial basis function--------------------------
	double wRBF[9] = {0.0};// here 9 is RBF::DIM

	void SetVolume();
//
	class DVDF
	{
	public:
		DVDF();
		DVDF(const DVDF &rhs);
		DVDF& operator=(const DVDF &rhs);
		~DVDF();
		//
		double tau = 0.0;
		double aBP = 0.0,bBP = 0.0,cBP = 0.0;
		inline void setxBP();
		//
		static int const nDVDF = 6;

		double *BarP = nullptr;
		double *BarP_x = nullptr;
		double *BarP_y = nullptr;
		double *Tilde = nullptr;
		double *Eq    = nullptr;
		double *So    = nullptr;
	private:
		int *token = nullptr;
		void assign(const DVDF& rhs);
		void allocate();
		void deallocate();
	};
//----------probability density function---------------------
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
//
	
//----------macro variables---------------------
	MacroQuantity * __restrict__ msq;
	inline MacroQuantity& MsQ(){return *msq;}
	inline MacroQuantity const & MsQ() const {return *msq;}
//-------------Factor---------------------
	double fBPLimiter = 1,gBPLimiter = 1;
	void Factor();
	void FactorAC();

private:
	int *use = nullptr;
	void assign(const Cell_2D &rhs);
};

#endif