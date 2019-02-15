#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include "MeshConstructFunc.h"

using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

void EraseParenthese(string &string_line,int &body)
{
	//only process lines whose first char is a left bracket;
	if('(' != string_line.front()) 
	{
		body = -1;
		return;
	}
	//pop space
	//while(' ' == string_line.back() || '\n' == string_line.back() || '\r' == string_line.back())
	while(std::isspace(string_line.back()))
		string_line.pop_back();
	//
	if(string_line.back() == '(') ++body;
	auto It = remove_if(string_line.begin(),string_line.end(),
									[](const char& c){return (c == '(' || c == ')');});
	//
	if(It == string_line.end()) body = -1;
	string_line.erase(It,string_line.end());
}
int StringToHex(int* const &ptrHexLine, istringstream &iss_line,int &count)
{
	for(;iss_line >> std::hex >>ptrHexLine[count];++count);
	return 0;
}
void HeadProcess(const int &index,int* const &ptrHexLine,int const &count)
{
	Allocate_Mesh(index,ptrHexLine,count);
}
void BodyProcess(const int &index,ifstream& InFile_Mesh,int* const &ptrHexLine)
{
	if(index == 10)
	{		
		ConstructNodes(ptrHexLine,InFile_Mesh);
	}
	else if(index == 12)
	{
		ConstructCells(ptrHexLine,InFile_Mesh);
	}
	else if(index == 13)
	{
		ConstructFaces(ptrHexLine,InFile_Mesh);
	}
	else if(index == 18)
	{
		ConstructPeriodicFaces(ptrHexLine,InFile_Mesh);
	}
	else
	{
		cout << "Invalid index during body processing" << endl;
		getchar();
	}
}
int MeshConstruct(const string &s)
{
	printSplitLine();
	cout<<"Mesh File Reading..."<<endl;
	ifstream InFile_Mesh;
	const char *phome = std::getenv("HOME");
	string shome(phome);
	shome += ("/Mesh/" + s +".cas");
	InFile_Mesh.open(shome);
	if(!InFile_Mesh)
	{
		printErrorLine();
		cout << "Mesh file open failed!!!  "<<endl;
		cout <<"name of mesh file : "<<shome<<endl;
		cout <<__FILE__ <<" ; "<<__LINE__<<" ; "<<__func__<<endl; 
		printErrorLine();
		exit(0);
	}
//
	int *ptrHexLine = new int[MeshPerLine];
	int count = 0, index = 0,body = 0;
	string string_line;
	while(getline(InFile_Mesh,string_line))
	{
		if(string_line.size() == 0) continue;
		cout << string_line << endl;
		body = 0;index = 0;	count = 0;
		for(int i = 0;i != MeshPerLine;++i)
			ptrHexLine[i] =  0;
		EraseParenthese(string_line,body);
		if(-1 == body) continue;		
		istringstream iss_line(string_line);
		if(!(iss_line >> index) || !(index == 10 || index == 12 || index == 13 || index == 18))
			continue;
		StringToHex(ptrHexLine,iss_line,count);
		if(0 == body)
			HeadProcess(index,ptrHexLine,count);
		else if(1 == body)
		{
			BodyProcess(index,InFile_Mesh,ptrHexLine);
		}
		else
		{
			cout << "unknown body indicator"<<endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
	}
	delete []ptrHexLine;
	cout <<"Mesh File Reading Done:"<<endl;
	cout << "Faces:" << Faces << " Nodes:" << Nodes << " Cells:" << Cells <<endl;
	printSplitLine(nl);
	return 0;
}

void AllocateFaces(int BoundFaceNum, Face_2D** &ptrBoundFace, vector<Face_2D*> &vecBoundFace)
{
	if (0 == BoundFaceNum) return;
	ptrBoundFace  = new Face_2D*[BoundFaceNum];
	for(int i = 0;i != BoundFaceNum;++i)
		ptrBoundFace[i] = vecBoundFace[i];
}
void FacesClassify()
{
	vector<Face_2D*> InteriorVec,WallVec,PeriodicVec,P_InletVec,P_OutletVec,V_InletVec,
					 SymmetryVec,P_FarfieldVec,BoundVec;//2,3,12-8,4,5,10,7,9,!2;
	for(int i = 0;i != Faces; ++i)
	{
		if(2 == FaceArray[i].bc_type)
			InteriorVec.push_back(&FaceArray[i]);
		else if(12 == FaceArray[i].bc_type || 8 == FaceArray[i].bc_type)
		{
			PeriodicVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(3 == FaceArray[i].bc_type)
		{
			WallVec.push_back(&FaceArray[i]);
			BoundVec.push_back(FaceArray+i);
		}
		else if(4 == FaceArray[i].bc_type)
		{
			P_InletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(5 == FaceArray[i].bc_type)
		{
			P_OutletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(7 == FaceArray[i].bc_type)
		{
			SymmetryVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(9 == FaceArray[i].bc_type)
		{
			P_FarfieldVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(10 == FaceArray[i].bc_type)
		{
			V_InletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else
		{
			cout << "FaceArray : " << i <<" unknown bc_type : " << FaceArray[i].bc_type
			<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
			getchar();
		}
	}
	if(static_cast<unsigned>(InteriorFaceNum) != InteriorVec.size())
	{
		cout << "FaceArray : " << InteriorFaceNum <<" interior Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(static_cast<unsigned>(WallFaceNum) != WallVec.size())
	{
		cout << "FaceArray : " << WallFaceNum <<" Wall Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(static_cast<unsigned>(PeriodicFaceNum) != PeriodicVec.size())
	{
		cout << "FaceArray : " << PeriodicFaceNum <<" Periodic Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(static_cast<unsigned>(P_InletFaceNum) != P_InletVec.size())
	{
		cout << "FaceArray : " << P_InletFaceNum <<" Pressure Inlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(static_cast<unsigned>(P_OutletFaceNum) != P_OutletVec.size())
	{
		cout << "FaceArray : " << P_OutletFaceNum <<" Pressure Outlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(static_cast<unsigned>(V_InletFaceNum) != V_InletVec.size())
	{
		cout << "FaceArray : " << V_InletFaceNum <<" Velocity Inlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(SymmetryVec.size() != static_cast<unsigned>(SymmetryFaceNum))
	{
		cout << "FaceArray : " << SymmetryFaceNum <<" Symmetry Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(P_FarfieldVec.size() != static_cast<unsigned>(P_FarfieldFaceNum))
	{
		cout << "FaceArray : " << P_FarfieldFaceNum <<" Pressure Farfield Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(BoundVec.size() != static_cast<unsigned>(BoundFaceNum))
	{
		cout << "FaceArray : " << BoundFaceNum <<" Boundary Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	printSplitLine();
	cout <<"Faces : "<<Faces<<'\n'
		 <<"BoundFaceNum : "<<BoundFaceNum<<'\n'
		 <<"InFaceNum : "<<InteriorFaceNum<<'\n'
		 <<"WallFaceNum : "<<WallFaceNum<<'\n'
		 <<"PeriodicFaceNum : "<<PeriodicFaceNum<<'\n'
		 <<"P_InletFaceNum : "<<P_InletFaceNum<<'\n'
		 <<"P_OutletFaceNum : "<<P_OutletFaceNum<<'\n'
		 <<"SymmetryFaceNum : "<<SymmetryFaceNum<<'\n'
		 <<"P_FarfieldFaceNum : "<<P_FarfieldFaceNum<<'\n'
		 <<"V_InletFaceNum : "<<V_InletFaceNum<<endl;
	printSplitLine(nl);
	if(InteriorFaceNum + WallFaceNum + PeriodicFaceNum + P_InletFaceNum + P_OutletFaceNum
		+ SymmetryFaceNum + P_FarfieldFaceNum + V_InletFaceNum != Faces)
	{
		_PRINT_ERROR_MSG_FLIP
		cout << "Faces != InteriorFaceNum + WallFaceNum"<<endl;
		getchar();
	}
	if(WallFaceNum + PeriodicFaceNum + P_InletFaceNum + P_OutletFaceNum
		+ SymmetryFaceNum + P_FarfieldFaceNum + V_InletFaceNum != BoundFaceNum)
	{
		_PRINT_ERROR_MSG_FLIP
		cout << "BoundFaceNum != ..."<<endl;
		getchar();
	}
	AllocateFaces(InteriorFaceNum,InteriorFaceA,InteriorVec);
	AllocateFaces(PeriodicFaceNum,PeriodicFaceA,PeriodicVec);
	AllocateFaces(WallFaceNum,WallFaceA,WallVec);
	AllocateFaces(P_InletFaceNum,P_InletFaceA,P_InletVec);
	AllocateFaces(P_OutletFaceNum,P_OutletFaceA,P_OutletVec);
	AllocateFaces(SymmetryFaceNum,SymmetryFaceA,SymmetryVec);
	AllocateFaces(P_FarfieldFaceNum,P_FarfieldFaceA,P_FarfieldVec);
	AllocateFaces(V_InletFaceNum,V_InletFaceA,V_InletVec);
	AllocateFaces(BoundFaceNum,BoundFaceA,BoundVec);
}
//
extern void sortFacesInThisCell(Cell_2D &cell);
void NeighbourCellConstruct()
{
	for(int n = 0;n < Cells;++n)
	{
		#ifdef _CARTESIAN_MESH_FLIP
			sortFacesInThisCell(CellArray[n]);
		#endif
		for(int i = 0;i < CellArray[n].celltype;++i)
		{
		CellArray[n].Cell_C[i] = ((CellArray[n].Face_C[i] -> owner ==  (&CellArray[n])) ? 
						CellArray[n].Face_C[i] -> neigh : CellArray[n].Face_C[i] -> owner);
		CellArray[n].signFlux[i] = ((CellArray[n].Face_C[i] -> owner ==  (&CellArray[n])) ? -1 : 1);
		}
	}
}
void setXiDotdS()
{
	LoopPS(Faces)
	LoopVS(Q)
	{
		FaceArray[n].xi_n_dS[k] = 
		FaceArray[n].Area*(xi_u[k]*FaceArray[n].Vx + xi_v[k]*FaceArray[n].Vy);
	}
}