#include "DUGKSDeclaration.h"
#include "nonidealEoS.h"
//

extern void update_Psi_x(Cell_2D *cellptr);

extern void update_Psi_y(Cell_2D *cellptr);
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
void Update_PseudoForce(Cell_2D &cell)
{
	update_Psi_x(&cell);
  update_Psi_y(&cell);
//
  cell.MsQ().Fx *= -cell.MsQ().scG*cell.MsQ().Psi;
  cell.MsQ().Fy *= -cell.MsQ().scG*cell.MsQ().Psi;
}