#ifndef _MESH_CONSTRUCTFUN_H_
#define _MESH_CONSTRUCTFUN_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "DUGKSDeclaration.h"

using std::string;
using std::ios;
using std::ofstream;
using std::istringstream;
using std::ostringstream;
using std::cout;
using std::setiosflags;
using std::setprecision;
using std::endl;

const int MeshPerLine = 5,NumPerCell = ARK::cellN;

const double digits_ = ARK::digits;

void ConstructCells(int* const &ptrHexLine,const int &count);

void Allocate_Mesh(const int &index,int* const &ptrHexLine,int const &count)
{
	int Num = ptrHexLine[2];
	if(index == 10)
	{		
		Nodes = Num;
		NodeArray = new Node_2D[Nodes]();
	}
	else if(index == 12)
	{
		Cells = Num;
		if(0 == ptrHexLine[0])
			CellArray = new Cell_2D[Cells]();
		else
		{
			ConstructCells(ptrHexLine,count);	
		}
	}
	else if(index == 13)
	{
		Faces = Num;
		FaceArray = new Face_2D[Faces]();
	}
	else
	{
		cout << "unknown head type  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<<endl;
	}
}
double validDigits(double &in)
{
	return std::round(in*digits_)/digits_;
}
void ConstructNodes(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	if(ptrHexLine[2] != Nodes)
	{
		printSplitLine();
		cout <<"Nodes number get from body doesn't match with head."<<nl;
		printErrorMessage(__LINE__,__FILE__,__func__);
		printSplitLine(nl);
	}
	InFile_Mesh >> std::dec;
	for(int i = 0;i != Nodes;++i)
	{
		InFile_Mesh >> NodeArray[i].xN >> NodeArray[i].yN;
		NodeArray[i].xN = validDigits(NodeArray[i].xN);
		NodeArray[i].yN = validDigits(NodeArray[i].yN);
		#ifdef _ZERO_NDEBUG_FLIP
		if(i%ZeroDebugControl == 0)
		std::cout <<std::setiosflags(ios::scientific)<<std::setprecision(15)
				  << NodeArray[i].xN << "  " << NodeArray[i].yN <<std::setprecision(6)
				  << resetiosflags(ios::scientific) << std::endl;
		#endif
	}
	cout <<"Nodes Constructed" << endl;
}
void ConstructCells(int* const &ptrHexLine,const int &count)
{
	if(5 != count) 
	{
		std::cerr <<"Error : " <<  __FILE__ <<" : in function "<< __func__
				  <<" at line " <<__LINE__ << endl;
		getchar();
		getchar();
		exit(0);
	}
	if(3 == ptrHexLine[count - 1] )
	{
		for(int i = 0;i != Cells;++i)
			CellArray[i].celltype = 4;
			
	}
	else if(1 == ptrHexLine[count - 1])
	{
		for(int i = 0;i != Cells;++i)
			CellArray[i].celltype = 3;
	}
	else
	{
		std::cerr << "Error : " <<  __FILE__ <<" : in function "<< __func__
				  <<" at line " <<__LINE__ << endl;
		std::cerr << "wrong element type" << endl;
		getchar();
		exit(0);
	}
}
void ConstructCells(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	if(ptrHexLine[2] != Cells)
	{
		printSplitLine();
		cout <<"Cells number get from body doesn't match with head."<<nl;
		printErrorMessage(__LINE__,__FILE__,__func__);
		printSplitLine(nl);
	}	
	int tmp;
	for(int i = 0;i != Cells;++i)
	{
		InFile_Mesh >> tmp;
		if(tmp == 1)
			CellArray[i].celltype = 3;
		else if(tmp == 3)
			CellArray[i].celltype = 4;
		else
		{
			cout << "wrong element index for 2D domain" << endl;
			getchar();
		}
	}
}
void PushNode(Cell_2D &Cell,const int &BegN,const int &EndN)
{
	if(nullptr == Cell.cellNodes[0])
	{
		// Cell.NodeX_C[0] = NodeX+BegN;
		// Cell.NodeX_C[1] = NodeX+EndN;
		// Cell.NodeY_C[0] = NodeY+BegN;
		// Cell.NodeY_C[1] = NodeY+EndN;
		Cell.cellNodes[0] = NodeArray + BegN;
		Cell.cellNodes[1] = NodeArray + EndN;
	}
	else if(nullptr == Cell.cellNodes[2] || nullptr == Cell.cellNodes[3])
	{
		if((NodeArray + BegN) == Cell.cellNodes[1])
		{
			Cell.cellNodes[2] = NodeArray + EndN;
		}
		else if((NodeArray + EndN) == Cell.cellNodes[1])
		{
			Cell.cellNodes[2] = NodeArray + BegN;
		}
		else if((NodeArray + BegN) == Cell.cellNodes[0])
		{
			Cell.cellNodes[3] = NodeArray + EndN;
		}
		else if((NodeArray + EndN) == Cell.cellNodes[0])
		{
			Cell.cellNodes[3] = NodeArray + BegN;
		}
	}
}
void PushFace(Cell_2D &Cell,Face_2D *ptrF)
{
	int i;
	for(i = 0;Cell.Face_C[i] != nullptr;++i);
	Cell.Face_C[i] = ptrF;
}
void PushCell(Cell_2D &Cell,Cell_2D *ptrF)
{
	int i;
	for(i = 0;Cell.Cell_C[i] != nullptr;++i);
	Cell.Cell_C[i] = ptrF;
}
void ConstructFacesCells(const int &FC,const int (&N_C)[MeshPerLine],const int &bc_type,const int &zone);
void ConstructFaces(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	const int &zone = ptrHexLine[0];
	const int &Beg_Face = ptrHexLine[1] -1;
	const int &End_Face = ptrHexLine[2];
	const int &bc_type = ptrHexLine[3];
	const int &face_type = ptrHexLine[4];
	int N_C[MeshPerLine] = {0,0,0,0,0};//set of node index and cell index 

	if(0 == face_type)
	{		
		InFile_Mesh >> std::hex;
		for(int i = Beg_Face;i < End_Face;++i)
		{
			InFile_Mesh >> N_C[0] >> N_C[1] >> N_C[2] >> N_C[3] >> N_C[4];
			ConstructFacesCells(i,N_C,bc_type,zone);
		}
		InFile_Mesh >> std::dec;	
	}
	else if(2 == face_type)
	{
		InFile_Mesh >> std::hex;
		for(int i = Beg_Face;i < End_Face;++i)
		{
			InFile_Mesh >> N_C[1] >> N_C[2] >> N_C[3] >> N_C[4];
			ConstructFacesCells(i,N_C,bc_type,zone);
		}
		InFile_Mesh >> std::dec;
	}
	else
	{
		cout <<__FILE__ <<"  "<<__LINE__<<"  "<<__func__<<endl;
		getchar();
	}
}
void ConstructFacesCells(const int &FC,const int (&N_C)[MeshPerLine],const int &bc_type,const int &zone)
{
	int BegN = N_C[1], EndN = N_C[2], lhsC = N_C[3], rhsC = N_C[4];
	--BegN; --EndN; --lhsC; --rhsC;
	FaceArray[FC].bc_type = bc_type;
	FaceArray[FC].zone = zone;
//
	FaceArray[FC].faceNodes[0] = NodeArray + BegN;
	FaceArray[FC].faceNodes[1] = NodeArray + EndN;
//
	FaceArray[FC].owner = CellArray + lhsC;
//	
	PushNode(CellArray[lhsC],BegN,EndN);
	PushFace(CellArray[lhsC],FaceArray+FC);

	if(bc_type == 2)
	{
		FaceArray[FC].neigh = CellArray + rhsC;
		PushNode(CellArray[rhsC],EndN,BegN);
		PushFace(CellArray[rhsC],FaceArray+FC);
		//PushCell(CellArray[rhsC],CellArray + lhsC);
		//PushCell(CellArray[lhsC],CellArray + rhsC);
		++InteriorFaceNum;
	}
	else if(12 == bc_type || 8 == bc_type)
	{
		++PeriodicFaceNum;
		++BoundFaceNum;
	}
	else if(3 == bc_type)
	{
		++WallFaceNum;
		++BoundFaceNum;
	}
	else if(4 == bc_type)
	{
		++P_InletFaceNum;
		++BoundFaceNum;
	}
	else if(5 == bc_type)
	{
		++P_OutletFaceNum;
		++BoundFaceNum;
	}
	else if(7 == bc_type)
	{
		++SymmetryFaceNum;
		++BoundFaceNum;
	}
	else if(9 == bc_type)
	{
		++P_FarfieldFaceNum;
		++BoundFaceNum;
	}
	else if(10 == bc_type)
	{
		++V_InletFaceNum;
		++BoundFaceNum;
	}
}
int MeshArea()
{
	double dL_Min = 1.0;
	printSplitLine();
	cout << "Traversing Faces : " << endl;
	for(int i = 0;i < Faces;++i)
	{
		FaceArray[i].SetArea();
		FaceArray[i].SetNormalV();
		if(FaceArray[i].Area < dL_Min) dL_Min = FaceArray[i].Area;
		#ifdef _ZERO_NDEBUG_FLIP
		if(i%ZeroDebugControl == 0)
		cout <<"Face : "<<i<<" "<<FaceArray[i].bc_type<<"  "
			 <<std::setiosflags(ios::scientific)<<std::setprecision(6)
			 <<FaceArray[i].xf <<" "<<FaceArray[i].yf<<" "
		     <<FaceArray[i].Area <<" "<<FaceArray[i].Vx<<" "<<FaceArray[i].Vy
		     <<resetiosflags(ios::scientific)<<endl;
		#endif
	}
	if(fabs(dL_Min - MinL) > 1.0E-15)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<std::setiosflags(ios::scientific)<<std::setprecision(15)
		<<"Mesh MinL : "<<dL_Min<<"  Const MinL : "<<MinL<<"    "<<fabs(dL_Min - MinL)<<endl;
		ofstream OutFile_Case("../FlowField/global/dMinL.ark",ios::app);
		if(!OutFile_Case)
		{
			_PRINT_ERROR_MSG_FLIP
			cout <<"dMinL.ark open failed."<<endl;
			getchar();
			return 0;
		}
		OutFile_Case <<std::setiosflags(ios::scientific)<<std::setprecision(15)
					 <<dL_Min<<"    =    dL_Min"<<endl;
		getchar();
	}
	cout << "SetArea Done, SetNormalV Done" <<endl;
	printSplitLine(nl);
	printSplitLine();
	cout << "Traversing Cells : " << endl;
	for(int i = 0;i != Cells;++i)
	{
		CellArray[i].SetVolume();
		#ifdef _ZERO_NDEBUG_FLIP
		if(i%ZeroDebugControl == 0)
		cout <<"Cell : "<<i<<" "<<std::setiosflags(ios::scientific)<<std::setprecision(6)
		 	 <<CellArray[i].xc <<" "<<CellArray[i].yc<<" "<< CellArray[i].volume
		 	 <<resetiosflags(ios::scientific)<<endl;
		#endif
	}
	cout << "SetVolume Done" <<endl;
	printSplitLine(nl);
	return 0;
}
void PeriodFaces(const int &F0,const int &F1)
{
	if(nullptr == FaceArray[F0].neigh && nullptr == FaceArray[F1].neigh)
	{
		FaceArray[F0].shadowF = FaceArray + F1; 
		FaceArray[F1].shadowF = FaceArray + F0; 
	}
	else
	{
		cout << FaceArray[F0].neigh <<"  " << FaceArray[F1].neigh <<endl;
		cout <<__FILE__ <<"  "<<__LINE__<<"  "<<__func__<<endl;
		getchar();
	}
}
void ConstructPeriodicFaces(int* const &ptrHexLine,std::ifstream &InFile_Mesh)
{
	const int &Beg = ptrHexLine[0] -1;
	const int &End = ptrHexLine[1];
	int F0 = 0;
	int F1 = 0;
	InFile_Mesh >> std::hex;
	for(int i = Beg;i < End;++i)
	{
		InFile_Mesh >> F0 >> F1;
		--F0;--F1;
		PeriodFaces(F0,F1);
	}
	InFile_Mesh >> std::dec;
}

#endif
