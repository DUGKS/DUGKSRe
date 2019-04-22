#ifndef _ZERO_FLIP_H
#define _ZERO_FLIP_H
#include "ZeroMacros.h"
//
//--------------------------DEBUG MACRO---------------------
// #ifndef _ZERO_NDEBUG_FLIP
// #define _ZERO_NDEBUG_FLIP
// #endif

// #ifndef _CARTESIAN_LS_DEBUG_FLIP
// #define _CARTESIAN_LS_DEBUG_FLIP
// #endif

// #ifndef _SINGLE_STEP_DEBUG_FLIP
// #define _SINGLE_STEP_DEBUG_FLIP
// #endif
//----------------------------------------------------------
//----------------------------------------------------------
#ifndef _CARTESIAN_MESH_FLIP
#define _CARTESIAN_MESH_FLIP
#endif

#ifndef _MESHTYPE_ARK
#define _MESHTYPE_ARK "Car"
#endif

// #ifndef _FLUX_SCHEME_UW_ARK
// #define _FLUX_SCHEME_UW_ARK "UW"	
// #define _FLUX_SCHEME_ARK _FLUX_SCHEME_UW_ARK
// #endif

#ifndef _FLUX_SCHEME_CD_ARK
#define _FLUX_SCHEME_CD_ARK "CD"
#define _FLUX_SCHEME_ARK _FLUX_SCHEME_CD_ARK
#endif

//BB = bounce back,NEE = non-equilibrium extrapolation,DS = diffusive scattering

#ifndef _BC_ARK
#define _BC_ARK "periodic"
#endif

#ifndef _EOS_MODEL_ARK
#define _EOS_MODEL_ARK "null"
#endif

#ifndef _MESHFILE_NAME_ARK
#define _MESHFILE_NAME_ARK "512_512_Car_Periodic_Square_LBM"
#endif
//----------------Boundary Condition Macro------------------

// #ifndef _P_INLET_4_BCS_FLIP
// #define _P_INLET_4_BCS_FLIP	4
// #endif

// #ifndef _P_OUTLET_5_BCS_FLIP
// #define _P_OUTLET_5_BCS_FLIP	5
// #endif

#ifndef _PERIODIC_12_8_BCs_FLIP
#define _PERIODIC_12_8_BCs_FLIP 12
#endif

// #ifndef _Wall_3_BCs_FLIP
// #define _Wall_3_BCs_FLIP	3
// 	// #ifndef _Wall_3_BCs_DS
// 	// #define _Wall_3_BCs_DS
// 	// #define _Wall_3_BCs_MODEL DS
// 	// #endif
// //
// 	#ifndef _Wall_3_BCs_NEE
// 	#define _Wall_3_BCs_NEE
// 	#define _Wall_3_BCs_MODEL "NEE"
// 	#endif
// //
// 	// #ifndef _Wall_3_BCs_BB
// 	// #define _Wall_3_BCs_BB
// 	// #define _Wall_3_BCs_MODEL "BB"
// 	// #endif
// #endif
//--------------------------------------multiphase model-------------------

// #ifndef _ARK_PSEUDOPSI_FLIP
// #define _ARK_PSEUDOPSI_FLIP
// #endif

#ifndef _ARK_ALLENCAHN_FLIP
#define _ARK_ALLENCAHN_FLIP
#endif

#ifndef _GRAD_SCHEME_9Points
#define _GRAD_SCHEME_9Points 
#define _ARK_LAPLACIAN_SCHEME "9Points" 
#endif

// #ifndef _GRAD_SCHEME_5Points
// #define _GRAD_SCHEME_5Points 
// #define _ARK_LAPLACIAN_SCHEME "5Points" 
// #endif
//
//-------------------------------Force model------------------------------

#ifndef _ARK_FORCE_FLIP
#define _ARK_FORCE_FLIP
#define _ARK_FORCE_MODEL "He-Shan-Doolean"
// #define _ARK_FORCE_H_Liang "none"
	#ifndef _ARK_FORCE_H_Liang2014
	#define _ARK_FORCE_H_Liang2014 
	#define _ARK_FORCE_H_Liang "H.Liang PRE 2014"
	#endif
// //
// 	// #ifndef _ARK_FORCE_H_Liang2018
// 	// #define _ARK_FORCE_H_Liang2018
// 	// #define _ARK_FORCE_H_Liang  "H.Liang PRE 2018"
// 	// #endif
#endif

// #ifndef _ARK_STRANGSPLIT_FLIP
// #define _ARK_STRANGSPLIT_FLIP
// #define _ARK_FORCE_MODEL "StrangSplitting"
// #endif

//-------------------------------Momentum Energy-------------------------------
// #ifndef _ARK_MOMENTUM_FLIP
// #define _ARK_MOMENTUM_FLIP
// #endif

#ifndef _ARK_ISOTHERMAL_FLIP
#define _ARK_ISOTHERMAL_FLIP
#endif
//
// #ifndef _ARK_LIMITER_FLIP
// #define _ARK_LIMITER_FLIP
// #endif
//----------------------------------------------------------------------------

// #ifndef _ARK_ENDTIME_FLIP
// #define _ARK_ENDTIME_FLIP
// #endif

// #ifndef _OUTPUT_L2NORM_ERROR_FLIP
// #define _OUTPUT_L2NORM_ERROR_FLIP
// #endif

// #ifndef _ARK_NOHUP_FLIP	//Flip on for server
// #define _ARK_NOHUP_FLIP
// #endif


//
//
//
//
//
//
//
//
//
#endif
