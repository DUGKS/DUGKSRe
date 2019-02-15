double const 

T_r = PseudoPotentialSC::Tr;

double const 

Tcr = 1.0*T_r;
//-----------------------------reduced Carnahan-Starling EoS--------------------------
double const 

arCS = 3.852462257,

brCS = 0.1304438842,

crCS = 2.785855166;

//-----------------------------modified Kaplun-Meshalkin EoS--------------------------
double const

cmKM = 2.78,

bmKM = 3.0 - cmKM,

amKM = 1.0/bmKM,

dmKM = (12.0*cmKM - 6*cmKM*cmKM + cmKM*cmKM*cmKM - 8.0)/(cmKM*bmKM);



double reduced_CS(double rho) //Carnahan-Starling
{
	double brCSrho = rho*brCS;
	return  rho*crCS*Tcr*(1.0 + brCSrho*(1 + brCSrho*(1.0-brCSrho)))
			/
			((1.0-brCSrho)*(1.0-brCSrho)*(1.0-brCSrho))
			- arCS*rho*rho;
}
double reduced_mKM(double rho) //modified Kaplun-Meshalkin EoS
{
	return cmKM*rho*Tcr*(1.0 + dmKM/(1.0/rho - bmKM)) - amKM*rho*rho;
}

//------------------------------VdW-------------------------------
double const

TcrVdW = 4.0/7,

aVdW = 9.0/49,

bVdW = 2.0/21,

T_VdW = T_r*TcrVdW;

double VanderWaals(double rho)
{
	return rho*R0*T_VdW/(1.0-bVdW*rho) - aVdW*rho*rho;
}
//-----------------------------Carnahan-Starling-------------------
double const

aCS = 1.0,

bCS = 4.0,

Tc_CS = 0.0943,

T_CS = T_r*Tc_CS;

double CarnahanStarling(double rho) //rho*(1+rho*(1-rho))
{
	return rho*R0*T_CS*(1.0+rho+rho*rho-rho*rho*rho)/((1.0-rho)*(1.0-rho)*(1.0-rho))
			-aCS*rho*rho;
}