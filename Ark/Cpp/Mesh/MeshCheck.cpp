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

const int NumPerCell = ARK::cellN;

//------------------------------------Mesh Check--------------------------------------------
void ShadowCellCheck(int BoundFaceNum,const Cell_2D* const &BoundShadowCA,const string& s)
{
	for(int i = 0;i < BoundFaceNum;++i)
	{
		if(nullptr == BoundShadowCA[i].ShadowC)
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<i<<" ShadowC != nullptr"<<endl;
			getchar();
		}
		if(nullptr == BoundShadowCA[i].Cell_C[0])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<i<<" Cell_C[0] == nullptr"<<endl;
			getchar();
		}
		if(nullptr == BoundShadowCA[i].Face_C[0])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<i<<" Face_C[0] == nullptr"<<endl;
			getchar();
		}
		if(nullptr == BoundShadowCA[i].BoundF)
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<i<<" BoundF == nullptr"<<endl;
			getchar();
		}
	}
}
bool nextNeighCheck(Cell_2D const*shadowCell,Cell_2D const*neighb,Cell_2D const*nextNeighb)
{
	double ax = shadowCell->xc - neighb->xc;
	double bx = neighb->xc - nextNeighb->xc;
	double ay = shadowCell->yc - neighb->yc;
	double by = neighb->yc - nextNeighb->yc;
	return (doubleEqual(ax,bx) && doubleEqual(ay,by));
}
// bool nextNeighCheckY(Cell_2D const*shadowCell,Cell_2D const*neighb,Cell_2D const*nextNeighb)
// {
// 	double a = shadowCell->yc - neighb->yc;
// 	double b = neighb->yc - nextNeighb->yc;
// 	return doubleEqual(a,b);
// }
void zoneCheck(const Face_2D &face)
{
	if(X_Beg == face.xf)
	{
		if(left != face.zone || left != face.neigh->zone)
		{
			_PRINT_SPLITLINE_ARK
			cout <<"left doesn't match with face.zone"<<endl;
			_PRINT_ERROR_MSG_FLIP
			exit(-1);
			_PRINT_SPLITLINE_ARK
		}
		if(!nextNeighCheck(face.neigh,face.owner,face.owner->Cell_C[0]))
		{
			_PRINT_SPLITLINE_ARK
			cout <<"wrong next Neighbour cell"<<endl;
			cout <<face.neigh->xc<<fs<<face.owner->xc<<fs<<face.owner->Cell_C[0]->xc<<endl;
			cout <<face.neigh->yc<<fs<<face.owner->yc<<fs<<face.owner->Cell_C[0]->yc<<endl;
			_PRINT_ERROR_MSG_FLIP
			_PRINT_SPLITLINE_ARK
			exit(-1);
		}
	}
	else if(X_End == face.xf)
	{
		if(right != face.zone || right != face.neigh->zone)
		{
			_PRINT_SPLITLINE_ARK
			cout <<"left doesn't match with face.zone"<<endl;
			_PRINT_ERROR_MSG_FLIP
			exit(-1);
			_PRINT_SPLITLINE_ARK
		}
		if(!nextNeighCheck(face.neigh,face.owner,face.owner->Cell_C[2]))
		{
			_PRINT_SPLITLINE_ARK
			cout <<"wrong next Neighbour cell"<<endl;
			cout <<face.neigh->xc<<fs<<face.owner->xc<<fs<<face.owner->Cell_C[2]->xc<<endl;
			cout <<face.neigh->yc<<fs<<face.owner->yc<<fs<<face.owner->Cell_C[2]->yc<<endl;
			_PRINT_ERROR_MSG_FLIP
			_PRINT_SPLITLINE_ARK
			exit(-1);
		}
	}
	else if(Y_Beg == face.yf)
	{
		if(bottom != face.zone || bottom != face.neigh->zone)
		{
			_PRINT_SPLITLINE_ARK
			cout <<"left doesn't match with face.zone"<<endl;
			_PRINT_ERROR_MSG_FLIP
			exit(-1);
			_PRINT_SPLITLINE_ARK
		}
		if(!nextNeighCheck(face.neigh,face.owner,face.owner->Cell_C[1]))
		{
			_PRINT_SPLITLINE_ARK
			cout <<"wrong next Neighbour cell"<<endl;
			cout <<face.neigh->xc<<fs<<face.owner->xc<<fs<<face.owner->Cell_C[1]->xc<<endl;
			cout <<face.neigh->yc<<fs<<face.owner->yc<<fs<<face.owner->Cell_C[1]->yc<<endl;
			_PRINT_ERROR_MSG_FLIP
			_PRINT_SPLITLINE_ARK
			exit(-1);
		}
	}
	else if(Y_End == face.yf)
	{
		if(top != face.zone || top != face.neigh->zone)
		{
			_PRINT_SPLITLINE_ARK
			cout <<"left doesn't match with face.zone"<<endl;
			_PRINT_ERROR_MSG_FLIP
			exit(-1);
			_PRINT_SPLITLINE_ARK
		}
		if(!nextNeighCheck(face.neigh,face.owner,face.owner->Cell_C[3]))
		{
			_PRINT_SPLITLINE_ARK
			cout <<"wrong next Neighbour cell"<<endl;
			cout <<face.neigh->xc<<fs<<face.owner->xc<<fs<<face.owner->Cell_C[3]->xc<<endl;
			cout <<face.neigh->yc<<fs<<face.owner->yc<<fs<<face.owner->Cell_C[3]->yc<<endl;
			_PRINT_ERROR_MSG_FLIP
			_PRINT_SPLITLINE_ARK
			exit(-1);
		}
	}
	else if(3 != face.zone)
	{
		_PRINT_SPLITLINE_ARK
		cout <<"inner faces zone != 3"<<endl;
		_PRINT_ERROR_MSG_FLIP
		exit(-1);
		_PRINT_SPLITLINE_ARK
	}
	if
	(
		left == right 	||	 	 
	 	left == top   	||
	 	left == bottom  ||
	 	right == top 	||
	 	right == bottom ||
	 	top == bottom
	)
	{
		_PRINT_SPLITLINE_ARK
		cout <<"left : "<<left<<fs<<"right : "<<right<<fs
			 <<"bottom : "<<bottom<<fs<<"top : "<<top
			 <<endl;
		_PRINT_ERROR_MSG_FLIP
		exit(-1);
		_PRINT_SPLITLINE_ARK
	}
}
int MeshOutput(const string& s)
{
	cout << "Mesh Output verifing..."<<endl;
	ofstream OutFile_Mesh;
	OutFile_Mesh.open("../MeshOutput/"+ s +".plt");
	if(!OutFile_Mesh)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		getchar();
		return -99999;
	}
	ostringstream fname,zonename,dataNE;
	zonename<<"ZONE T = \"Mesh\"\n";
	dataNE<<"Nodes="<<Nodes<<", Elements="<<Cells<<", ZONETYPE=FEQuadrilateral\n";
	string tecformat[5]={"VARIABLES = \"X\",\"Y\"\n",
						zonename.str().c_str(),dataNE.str().c_str(),"DATAPACKING=BLOCK\n",
						"VarLocation=([1-2]=NODAL)\n"};
	OutFile_Mesh << tecformat[0]<<tecformat[1]<<tecformat[2]<<tecformat[3]<<tecformat[4];
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_Mesh << NodeArray[i].xN <<"   ";
		if((i+1)%16 == 0)
			OutFile_Mesh << endl;
	}
	OutFile_Mesh << endl;
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_Mesh << NodeArray[i].yN <<"   ";
		if((i+1)%16 == 0)
			OutFile_Mesh << endl;
	}
	OutFile_Mesh << endl;
	for(int i = 0;i != Cells;++i)
	{
		OutFile_Mesh << MeshIndex(CellArray[i].cellNodes[0] , NodeArray)<< " "
					 << MeshIndex(CellArray[i].cellNodes[1] , NodeArray)<< " "
					 << MeshIndex(CellArray[i].cellNodes[2] , NodeArray)<< " " 
					 << MeshIndex(CellArray[i].cellNodes[3] , NodeArray)<< endl;
	}
	cout <<"Mesh Output Done" <<endl;
	OutFile_Mesh.close();
	cout <<"Face Output verifing..."<<endl;
	ofstream OutFile_Face;
	OutFile_Face.open("../MeshOutput/13_" + s + ".dat");
	if(!OutFile_Face)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		getchar();
	}
	OutFile_Face << std::hex;
	for(int i = 0;i != InteriorFaceNum;++i)
	{
		OutFile_Face << MeshIndex(FaceArray[i].faceNodes[0] , NodeArray)<<" "
					 << MeshIndex(FaceArray[i].faceNodes[1] , NodeArray)<<" "
					 << MeshIndex(FaceArray[i].owner ,  CellArray)<<" "
					 << MeshIndex(FaceArray[i].neigh , CellArray)<<endl;
	}
	OutFile_Face.close();
//
	OutFile_Face.open("../MeshOutput/18_" + s +".dat");
	if(!OutFile_Face)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		getchar();
	}
	OutFile_Face << std::hex;
	for(int i = 0;i != Faces;++i)
	{
		if(8 == FaceArray[i].bc_type)
		OutFile_Face << MeshIndex(FaceArray[i].shadowF , FaceArray)<<" "
					 << i + 1 <<endl;
	}
	OutFile_Face.close();
	cout <<"Face Output Done" <<endl;
	cout <<"Node Output verifing..."<<endl;
	ofstream OutFile_Node;
	OutFile_Node.open("../MeshOutput/10_" + s + ".dat");
	if(!OutFile_Node)
	{
		cout << __FILE__ <<"  " << __func__ <<"  " << __LINE__ 
			 <<"  "<<"file open failed" << endl; 
		getchar();
	}
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_Node <<std::setiosflags(ios::scientific) <<std::setprecision(16)
					 <<NodeArray[i].xN<<"  "<<NodeArray[i].yN << endl;
	}
	OutFile_Node.close();
	cout <<"Node Output Done"<<endl;
	return 0;
}
void faceCellsCheck(Face_2D const &face)
{
	if(face.Vx == 0)
	{
		if(!doubleEqual(face.owner->yc , face.faceCells[0]->yc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[0] yc doesn't match"<<endl;
			printErrorLine(nl);
			exit(-1);
		}
		if(!doubleEqual(face.owner->yc , face.faceCells[2]->yc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[2] yc doesn't match"<<endl;
			printErrorLine(nl);
			exit(-1);
		}
		if(!doubleEqual(face.neigh->yc , face.faceCells[1]->yc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[1] yc doesn't match"<<endl;
			printErrorLine(nl);
			exit(-1);
		}
		if(!doubleEqual(face.neigh->yc , face.faceCells[3]->yc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[3] yc doesn't match"<<endl;
			printErrorLine(nl);
			exit(-1);
		}
	}
	else if(face.Vy == 0)
	{
		if(!doubleEqual(face.owner->xc , face.faceCells[0]->xc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[0] xc doesn't match"<<endl;
			printErrorLine(nl);
			exit(-1);
		}
		if(!doubleEqual(face.owner->xc , face.faceCells[2]->xc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[2] xc doesn't match"<<endl;
			printErrorLine(nl);
			exit(-1);
		}
		if(!doubleEqual(face.neigh->xc , face.faceCells[1]->xc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[1] xc doesn't match"<<endl;
			printErrorLine(nl);
			exit(-1);
		}
		if(!doubleEqual(face.neigh->xc , face.faceCells[3]->xc))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : faceCells[3] xc doesn't match"<<endl;
			cout <<face.neigh->xc<<"    "<<face.faceCells[3]->xc<<"    "
				 <<(face.neigh->xc-face.faceCells[3]->xc)<<endl;
			printErrorLine(nl);
			exit(-1);
		}
	}
	else
	{
		printErrorLine();
		_PRINT_ERROR_MSG_FLIP
		cout <<"Vx != 0 && Vy != 0"<<nl;
		printErrorLine(nl);
		exit(-1);
	}	
}
int MeshCheck()
{
	_PRINT_SPLITLINE_ARK
	cout <<"Face Checking..." <<endl;
	for(int i = 0;i != Faces;++i)
	{
		/*if(FaceArray[i].bc_type != 12 &&  FaceArray[i].bc_type != 8)
		cout << "FaceArray : " << i <<"   bc_type : " << FaceArray[i].bc_type <<std::hex
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[0] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[1] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].owner , CellArray)
			 <<"  "<< MeshIndex(FaceArray[i].neigh , CellArray)<< std::dec <<endl;
		else
		cout << "FaceArray : " << i <<"   bc_type : " << FaceArray[i].bc_type <<std::hex
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[0] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].NodeX_F[1] , NodeX)
			 <<"  "<< MeshIndex(FaceArray[i].owner , CellArray)
			 <<"  "<< MeshIndex(FaceArray[i].neigh->ShadowC , CellArray)<< std::dec <<endl;*/
		#ifdef _CARTESIAN_MESH_FLIP
		if(FaceArray[i].Vx * FaceArray[i].Vy != 0.0)
		{
			cout << "FaceArray : " << i <<" Vx : " << FaceArray[i].Vx
				 <<" Vx : " << FaceArray[i].Vy <<" Area: "
				 <<std::setiosflags(ios::scientific) <<std::setprecision(16)
				 <<FaceArray[i].Area
				 <<std::resetiosflags(ios::scientific) <<std::setprecision(6)
				 << endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
		}
		if(!doubleEqual(FaceArray[i].Area,MinL))
		{
			cout << "FaceArray : " << i <<" Vx : " << FaceArray[i].Vx
				 <<" Vx : " << FaceArray[i].Vy <<" Area: "
				 <<std::setiosflags(ios::scientific) <<std::setprecision(16)
				 <<FaceArray[i].Area
				 <<std::resetiosflags(ios::scientific) <<std::setprecision(6)
				 << endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
		}
		zoneCheck(FaceArray[i]);
		#endif
		if(2 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"Cell = -1" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"construct shadow faces for interior faces"<< endl;
				getchar();
			}
		}
		else if(3 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"owner == nullptr || neigh == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for wall boundaries"<< endl;
				getchar();
			}
		}
		else if(8 == FaceArray[i].bc_type || 12 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"Cell = -1" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"Failed to construct shadowF faces for 8 & 12"<< endl;
				getchar();
			}
		}
		else if(4 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"owner == nullptr || neigh == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Pressure Inlet boundaries"<< endl;
				getchar();
			}
		}
		else if(5 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"owner == nullptr || neigh == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Pressure Outlet boundaries"<< endl;
				getchar();
			}
		}
		else if(7 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"owner == nullptr || neigh == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Symmetry boundaries"<< endl;
				getchar();
			}
		}
		else if(9 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"owner == nullptr || neigh == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Pressure Farfield boundaries"<< endl;
				getchar();
			}
		}
		else if(10 == FaceArray[i].bc_type)
		{
			if(FaceArray[i].owner == nullptr || FaceArray[i].neigh == nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout << "FaceArray : " << i <<" bc_type : " << FaceArray[i].bc_type
					 <<"owner == nullptr || neigh == nullptr" << endl;
				getchar();
			}
			if(FaceArray[i].shadowF != nullptr)
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"constructed shadow faces for Velocity Inlet boundaries"<< endl;
				getchar();
			}
		}
		else
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "FaceArray : " << i <<"unknown bc_type : " << FaceArray[i].bc_type << endl;
			getchar();
		}
	}
	LoopPS(Faces)
	{
		Face_2D &face = FaceArray[n];
		double VxCell = face.neigh->xc - face.owner->xc;
		double VyCell = face.neigh->yc - face.owner->yc;
		if((face.Vx*VxCell + face.Vy*VyCell) <= 0)
		{
			_PRINT_SPLITLINE_ARK
			_PRINT_ERROR_MSG_FLIP
			cout << "Fatal Error : normal vector doesn't conform with owner-neigh relation"
				 << endl;
			getchar();
			_PRINT_SPLITLINE_ARK
			exit(-1);
		}
	}
	cout <<"Face Check Done" <<endl;
	printSplitLine(nl);
	//
	//!Cell Checking
	//
	printSplitLine();
	cout <<"Cell Checking..." <<endl;
	for(int i = 0;i < Cells;++i)
	{
		for(int iF = 0;iF < CellArray[i].celltype;++iF)
		{
			if(nullptr == CellArray[i].Cell_C[iF])
			{
				_PRINT_ERROR_MSG_FLIP
				cout<<"Face Addr : "<<CellArray[i].Face_C[iF]<<'\n'
					<<"Face owner : "<<CellArray[i].Face_C[iF]->owner<<'\n'
					<<"Face neigh : "<<CellArray[i].Face_C[iF]->neigh<<'\n'
					<<"Cell Addr : "<<&CellArray[i]<<'\n';
				exit(-1);
			}
			if(CellArray[i].Face_C[iF]->owner == &CellArray[i])
			{
				if(CellArray[i].Face_C[iF]->neigh != CellArray[i].Cell_C[iF])
				{
					_PRINT_ERROR_MSG_FLIP
					cout <<"Construct neighbour cells failed : "<<"Face_C[iF] != Cell_C[iF]"<<endl
						 <<"CellArray : " << i <<endl;
					cout<<"Face_C[iF]->neigh : "<<CellArray[i].Face_C[iF]->neigh<<endl
						<<"Cell_C[iF]                      : "<<CellArray[i].Cell_C[iF]<<endl;
					getchar();
				}
			}
			else if(CellArray[i].Face_C[iF]->neigh == &CellArray[i])
			{
				if(CellArray[i].Face_C[iF]->owner != CellArray[i].Cell_C[iF])
				{
					_PRINT_ERROR_MSG_FLIP
					cout <<"Construct neighbour cells failed : "<<"Face_C[iF] != Cell_C[iF]"<<endl
						 <<"CellArray : " << i <<endl;
					cout<<"Face_C[iF]->owner : "<<CellArray[i].Face_C[iF]->owner<<endl
						<<"Cell_C[iF]                      : "<<CellArray[i].Cell_C[iF]<<endl;
					getchar();
				}
			}
			else
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Construct neighbour cells failed"<<endl<<"CellArray : " << i <<endl;
				getchar();
			}
		}
		if(nullptr != CellArray[i].ShadowC)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"Inner Cell has a shadow Cell : " << i <<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if(3 == CellArray[i].celltype)
		{
			for(int k = 0;k != NumPerCell;++k)
			if(nullptr == CellArray[i].cellNodes[k])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Node Index : " << k <<" = nullptr"<<endl;
				getchar();
			}

			if((CellArray[i].cellNodes[2]->xN != CellArray[i].cellNodes[3]->xN) ||
			   (CellArray[i].cellNodes[2]->yN != CellArray[i].cellNodes[3]->yN))
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Addrs of last two Nodes don't equal"<<endl;
				getchar();
			}

			for(int k = 0;k != NumPerCell - 1;++k)
			if(nullptr == CellArray[i].Cell_C[k])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Cell Index : " << k <<" = nullptr"<<endl;
				getchar();
			}

			if(nullptr != CellArray[i].Cell_C[NumPerCell - 1])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Triangular Cell "<<"CellArray : " << i 
					 <<"  Cell Index : " << NumPerCell - 1 <<" != nullptr"<<endl;
				getchar();
			}
		}
		else if(4 == CellArray[i].celltype)
		{
			for(int k = 0;k != NumPerCell;++k)
			if(nullptr == CellArray[i].cellNodes[k])
			{
				_PRINT_ERROR_MSG_FLIP
				cout <<"Quadrilateral Cell "<<"CellArray : " << i 
					 <<"  Node Index : " << k <<" == nullptr"<<endl;
				getchar();
			}
			for(int k = 0;k != NumPerCell;++k)
			if(nullptr == CellArray[i].Cell_C[k])
			{
				cout <<"Quadrilateral Cell "<<"CellArray : " << i 
					 <<"  Cell Index : " << k <<" == nullptr"<<endl;
				getchar();
			}
		}
		else
		{
			cout << "CellArray : " << i <<" unknown celltype : "
			     << CellArray[i].celltype << endl;
			getchar();
		}
		#ifdef _CARTESIAN_MESH_FLIP
		if(!doubleEqual(MinL*MinL,CellArray[i].volume))
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "CellArray : " << i <<" unmatch cell volume : "
				 <<std::setiosflags(ios::scientific) <<std::setprecision(16)
			     << MinL*MinL <<"  "<< CellArray[i].volume <<"    "<<MinL*MinL - CellArray[i].volume
			     <<endl;
			getchar();
		}
		#endif
	}
	for(int n = 0;n < PeriodicFaceNum;++n)
	{
		const string &s = "PeriodicShadowCA";
		if(nullptr == PeriodicShadowCA[n].ShadowC)
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<n<<" ShadowC != nullptr"<<endl;
			getchar();
		}
		for(int k = 1;k < NumPerCell;++k)		
		if(nullptr != PeriodicShadowCA[n].Cell_C[k])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<n<<" Cell_C[k>0] != nullptr"<<endl;
			getchar();
		}
		for(int k = 1;k < NumPerCell;++k)
		if(nullptr != PeriodicShadowCA[n].Face_C[k])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << s <<" : i = "<<n<<" Face_C[k] == nullptr"<<endl;
			getchar();
		}
	}
	ShadowCellCheck(WallFaceNum,WallShadowCA,"WallShadowCA");
	ShadowCellCheck(P_InletFaceNum,P_InletShadowCA,"P_InletShadowCA");
	ShadowCellCheck(P_OutletFaceNum,P_OutletShadowCA,"P_OutletShadowCA");
	ShadowCellCheck(SymmetryFaceNum,SymmetryShadowCA,"SymmetryShadowCA");
	ShadowCellCheck(P_FarfieldFaceNum,P_FarfieldShadowCA,"P_FarfieldShadowCA");
	ShadowCellCheck(V_InletFaceNum,V_InletShadowCA,"V_InletShadowCA");
	cout <<"Cell Check Done" <<endl;
	_PRINT_SPLITLINE_ARK
	cout <<nl;
	#ifdef _CARTESIAN_MESH_FLIP
	printSplitLine();
	cout <<"Cartesian Cell Checking..." <<endl;
	//!------------Check if faceFaces were in right order or not---------------------
	#if !defined _Wall_3_BCs_FLIP
	LoopPS(Faces)
	{
		Face_2D const &face = FaceArray[n];
		if
		(
			doubleEqual(face.faceFaces[1]->xf-face.xf,dx)
			&& 
		 	doubleEqual(face.faceFaces[1]->xf-face.xf+Lx,dx)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[1]->xf - face.xf : "<<fs
				 <<face.faceFaces[1]->xf - face.xf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if
		(
			doubleEqual(face.xf-face.faceFaces[5]->xf,dx)
			&& 
		 	doubleEqual(face.xf-face.faceFaces[5]->xf+Lx,dx)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"face.xf-face.faceFaces[5]->xf : "<<fs
				 <<face.xf-face.faceFaces[5]->xf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if
		(
			doubleEqual(face.faceFaces[3]->yf-face.yf,dy)
			&& 
		 	doubleEqual(face.faceFaces[3]->yf-face.yf+Ly,dy)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"face.yf-face.faceFaces[3]->yf : "<<fs
				 <<face.yf-face.faceFaces[3]->yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if
		(
			doubleEqual(face.yf-face.faceFaces[7]->yf,dy)
			&& 
		 	doubleEqual(face.yf-face.faceFaces[7]->yf+Ly,dy)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"face.yf-face.faceFaces[7]->yf : "<<fs
				 <<face.yf-face.faceFaces[7]->yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if(!doubleEqual(face.faceFaces[1]->yf,face.yf))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[1]->yf != face.yf"<<endl;
			cout <<"faceFaces[1]->yf : "<<face.faceFaces[1]->yf<<fs
				 <<"face.yf :"<<face.yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if(!doubleEqual(face.faceFaces[3]->xf,face.xf))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[3]->xf != face.xf"<<endl;
			cout <<"faceFaces[3]->xf : "<<face.faceFaces[3]->xf<<fs
				 <<"face.xf :"<<face.xf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if(!doubleEqual(face.faceFaces[5]->yf,face.yf))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[5]->yf != face.yf"<<endl;
			cout <<"faceFaces[5]->yf : "<<face.faceFaces[5]->yf<<fs
				 <<"face.yf :"<<face.yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if(!doubleEqual(face.faceFaces[7]->xf,face.xf))
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[7]->xf != face.xf"<<endl;
			cout <<"faceFaces[7]->xf : "<<face.faceFaces[7]->xf<<fs
				 <<"face.xf :"<<face.xf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if
		(
			!doubleEqual(face.faceFaces[2]->xf,face.faceFaces[1]->xf)
			||
			!doubleEqual(face.faceFaces[2]->yf,face.faceFaces[3]->yf)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[2] doesn't match"<<endl;
			cout <<"faceFaces[2]->xf : "<<face.faceFaces[2]->xf<<fs
				 <<"faceFaces[1]->xf : "<<face.faceFaces[1]->xf<<fs
				 <<"faceFaces[2]->yf : "<<face.faceFaces[2]->yf<<fs
				 <<"faceFaces[3]->yf : "<<face.faceFaces[3]->yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if
		(
			!doubleEqual(face.faceFaces[4]->xf,face.faceFaces[5]->xf)
			||
			!doubleEqual(face.faceFaces[4]->yf,face.faceFaces[3]->yf)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[4] doesn't match"<<endl;
			cout <<"face.xf : "<<face.xf<<fs<<"face.yf : "<<face.yf<<endl;
			cout <<"faceFaces[4]->xf : "<<face.faceFaces[4]->xf<<fs
				 <<"faceFaces[5]->xf : "<<face.faceFaces[5]->xf<<fs
				 <<"faceFaces[4]->yf : "<<face.faceFaces[4]->yf<<fs
				 <<"faceFaces[3]->yf : "<<face.faceFaces[3]->yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if
		(
			!doubleEqual(face.faceFaces[6]->xf,face.faceFaces[5]->xf)
			||
			!doubleEqual(face.faceFaces[6]->yf,face.faceFaces[7]->yf)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[6] doesn't match"<<endl;
			cout <<"faceFaces[6]->xf : "<<face.faceFaces[6]->xf<<fs
				 <<"faceFaces[5]->xf : "<<face.faceFaces[5]->xf<<fs
				 <<"faceFaces[6]->yf : "<<face.faceFaces[6]->yf<<fs
				 <<"faceFaces[7]->yf : "<<face.faceFaces[7]->yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
		if
		(
			!doubleEqual(face.faceFaces[8]->xf,face.faceFaces[1]->xf)
			||
			!doubleEqual(face.faceFaces[8]->yf,face.faceFaces[7]->yf)
		)
		{
			printErrorLine();
			_PRINT_ERROR_MSG_FLIP
			cout <<"faceFaces[8] doesn't match"<<endl;
			cout <<"faceFaces[8]->xf : "<<face.faceFaces[8]->xf<<fs
				 <<"faceFaces[1]->xf : "<<face.faceFaces[1]->xf<<fs
				 <<"faceFaces[8]->yf : "<<face.faceFaces[8]->yf<<fs
				 <<"faceFaces[7]->yf : "<<face.faceFaces[7]->yf<<endl;
			printErrorLine('\n');
			exit(-1);
		}
	}
	#endif
	LoopPS(Faces)
	{
		faceCellsCheck(FaceArray[n]);
	}
	for(int i = 1;i < Nxp1;++i)
	{
		if(0 == CarCellArray[i][Nyp1] || 0 == CarCellArray[i][0])
		{
			_PRINT_SPLITLINE_ARK
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : CarCellArray for i : "<<i<<fs<<"j : 0 or Nyp1"<<" is 0"<<nl;
			_PRINT_SPLITLINE_ARK
			exit(-1);
		}
	}
	for(int j = 1;j < Nyp1;++j)
	{
		if(0 == CarCellArray[Nxp1][j] || 0 == CarCellArray[0][j])
		{
			_PRINT_SPLITLINE_ARK
			_PRINT_ERROR_MSG_FLIP
			cout <<"Fatal Error : CarCellArray for i : 0 or Nxp1"<<fs<<"j : "<<j<<" is 0"<<nl;
			_PRINT_SPLITLINE_ARK
			exit(-1);
		}
	}
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		if(
			CarCellArray[i][j]->Cell_C[0] != CarCellArray[i+1][j]
		  )
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_C[0] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout <<CarCellArray[i][j]->Cell_C[0] <<"    "<< CarCellArray[i+1][j]<<endl;
			getchar();
		}
		if(CarCellArray[i][j]->Cell_C[1] != CarCellArray[i][j+1])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_C[1] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout <<CarCellArray[i][j]->Cell_C[1] << fs << CarCellArray[i][j+1]<<endl;
			getchar();
		}
		if(CarCellArray[i][j]->Cell_C[2] != CarCellArray[i-1][j])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_C[2] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout <<CarCellArray[i][j]->Cell_C[2] << CarCellArray[i-1][j]<<nl;
			getchar();
		}
		if(CarCellArray[i][j]->Cell_C[3] != CarCellArray[i][j-1])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_C[3] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout <<CarCellArray[i][j]->Cell_C[3]<<"    "<<CarCellArray[i][j-1]<<endl;
			getchar();
		}
		if(CarCellArray[i][j]->Cell_Diag[0] != CarCellArray[i+1][j+1])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_Diag[0] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout <<CarCellArray[i][j]->Cell_Diag[0] <<fs<< CarCellArray[i+1][j+1]<<nl;
			getchar();
		}
		if(CarCellArray[i][j]->Cell_Diag[1] != CarCellArray[i-1][j+1])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_Diag[1] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout << CarCellArray[i][j]->Cell_Diag[1] <<fs<< CarCellArray[i-1][j+1]<<nl;
			getchar();
		}
		if(CarCellArray[i][j]->Cell_Diag[2] != CarCellArray[i-1][j-1])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_Diag[2] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout << CarCellArray[i][j]->Cell_Diag[2]<< fs<< CarCellArray[i-1][j-1]<<nl;
			getchar();
		}
		if(CarCellArray[i][j]->Cell_Diag[3] != CarCellArray[i+1][j-1])
		{
			_PRINT_ERROR_MSG_FLIP
			cout << "Cell_Diag[3] doesn't match"<<endl;
			cout <<"i : "<<i<<" j : "<<j<<endl;
			cout <<CarCellArray[i][j]->Cell_Diag[3] <<fs<< CarCellArray[i+1][j-1]<<nl;
			getchar();
		}
	}
	if(LeftRight)
	{
		if(PeriodicShadowC_NW->ShadowC != CarCellArray[Nx][Nyp1])
		{
			cout <<"NW shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
		if(PeriodicShadowC_NE->ShadowC != CarCellArray[1][Nyp1])
		{
			cout <<"NE shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
		if(PeriodicShadowC_SW->ShadowC != CarCellArray[Nx][0])
		{
			cout <<"SW shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
		if(PeriodicShadowC_SE->ShadowC != CarCellArray[1][0])
		{
			cout <<"SE shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
	}
	else if(TopBottom)
	{
		if(PeriodicShadowC_NW->ShadowC != CarCellArray[0][1])
		{
			cout <<"NW shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
		if(PeriodicShadowC_NE->ShadowC != CarCellArray[Nxp1][1])
		{
			cout <<"NE shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
		if(PeriodicShadowC_SW->ShadowC != CarCellArray[0][Ny])
		{
			cout <<"SW shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
		if(PeriodicShadowC_SE->ShadowC != CarCellArray[Nxp1][Ny])
		{
			cout <<"SE shadow cell doesn't match";
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
	}
	// else
	// {
	// 	cout <<"warning : no periodic Boundary";
	// 	getchar();
	// }
	if(nullptr == PeriodicShadowC_SW->ShadowC->ShadowC)
	{
		cout <<"SW->ShadowC->ShadowC doesn't exit";
		_PRINT_ERROR_MSG_FLIP
		getchar();
		exit(0);
	}
	if(nullptr == PeriodicShadowC_SE->ShadowC->ShadowC)
	{
		cout <<"SE->ShadowC->ShadowC doesn't exit";
		_PRINT_ERROR_MSG_FLIP
		getchar();
		exit(0);
	}
	if(nullptr == PeriodicShadowC_NW->ShadowC->ShadowC)
	{
		cout <<"NW->ShadowC->ShadowC doesn't exit";
		_PRINT_ERROR_MSG_FLIP
		getchar();
		exit(0);
	}
	if(nullptr == PeriodicShadowC_NE->ShadowC->ShadowC)
	{
		cout <<"NE->ShadowC->ShadowC doesn't exit";
		_PRINT_ERROR_MSG_FLIP
		getchar();
		exit(0);
	}
	cout <<"Cartesian Cell Check Done" <<endl;
	_PRINT_SPLITLINE_ARK
	cout<<nl;
	#endif
	return 0;
}