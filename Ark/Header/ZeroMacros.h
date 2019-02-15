#ifndef _ZERO_MACRO_ARK
#define _ZERO_MACRO_ARK

#ifndef _PRINT_ERROR_MSG_FLIP
#define _PRINT_ERROR_MSG_FLIP  cout<<"File : "<<__FILE__<<"  Line : "\
<<__LINE__<<"  fun : "<<__func__<<'\n';
#endif

#ifndef _PRINT_SPLITLINE_ARK
#define _PRINT_SPLITLINE_ARK	cout<<"----------------------------------"<<'\n';
#endif

#define LoopPS(MESHNUM) for(int n=0;n<MESHNUM;++n)

#define LoopPSB(nBoundary) for(int nB = 0;nB < nBoundary;++nB)

#define LoopVS(VNUM) for(int k = 0;k < VNUM;++k)

#define nl '\n'

#define fs "    "

#define Info cout<<setiosflags(ios::scientific)<<setprecision(12)

#define OUTCASE "../FlowField/"

double const

infinitesimal = 1.0E-14,

PI = 3.141592653589793;

#endif