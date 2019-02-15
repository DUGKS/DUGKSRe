#include "DUGKSDeclaration.h"
#include "nonidealEoS.h"
//

extern void update_Psi_x(Cell_2D *cellptr);

extern void update_Psi_y(Cell_2D *cellptr);

extern void update_Psi_x(Face_2D *faceptr);

extern void update_Psi_y(Face_2D *faceptr);
//
using PseudoPotentialSC::K;

void Update_PseudoPsi(Cell_2D &cell)
{
  //cell.p = K * reduced_mKM(cell.Rho);
  //cell.p = K * VanderWaals(cell.Rho);
  cell.MsQ().p = K * CarnahanStarling(cell.MsQ().Rho);
//
  cell.MsQ().calcPsi();
}
void Update_PseudoPsi(Face_2D &face)
{
  //cell.p = K * reduced_mKM(cell.Rho);
  //cell.p = K * VanderWaals(cell.Rho);
  face.MsQh().p = K * CarnahanStarling(face.MsQh().Rho);
//
  face.MsQh().calcPsi();
}
void Update_PseudoForce(Cell_2D &cell)
{
	update_Psi_x(&cell);
  update_Psi_y(&cell);
//
  cell.MsQ().Fx *= -cell.MsQ().scG*cell.MsQ().Psi;
  cell.MsQ().Fy *= -cell.MsQ().scG*cell.MsQ().Psi;
}
void Update_PseudoForce(Face_2D &face)
{
  update_Psi_x(&face);
  update_Psi_y(&face);
  //
  face.MsQh().Fx *= -face.MsQh().scG*face.MsQh().Psi;
  face.MsQh().Fy *= -face.MsQh().scG*face.MsQh().Psi;
}