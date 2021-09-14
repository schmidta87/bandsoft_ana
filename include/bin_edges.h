#ifndef __BIN_EDGES_H__
#define __BIN_EDGES_H__


// Bins for tagged/inclusive analysis
const int 	bins_Q2 		= 1;
const double 	Q2Bins[bins_Q2 + 1] 	= {2.0, 6.0};	// irregular bins

const int 	bins_Pt			= 2;
const double 	Pt_step 		= 0.1;
const double 	Pt_min 			= 0.0;
const double	Pt_max			= Pt_min + Pt_step*bins_Pt;

const int 	bins_Xb 		= 6;
const double 	Xb_step 		= 0.025; 
const double 	Xb_min[bins_Q2]		= {0.15};//{0.15, 0.2, 0.25};

const int 	bins_As 		= 3;
const double 	As_step 		= 0.1;
const double 	As_min			= 1.3;
const double 	As_max 			= As_min + As_step*bins_As; 

const int 	bins_Xp 		= 3;
const double 	Xp_step 		= 0.1; 
const double 	Xp_min			= 0.2;
const double 	Xp_max			= 0.5;
const double 	Xp_Ref[bins_Q2]		= {0.4};//{0.4,0.5,0.6};

const int 	bins_Virt 		= 13;
const double 	Virt_step 		= 0.02; 
const double 	Virt_min			= -0.45;
const double 	Virt_max			= -0.19;
const double 	Virt_Ref[bins_Xp]		= {-0.25,-0.25,-0.25};//{-0.22,-0.26,-0.36};

// Bins for checks on inclusive acceptance

// Pe = electron momentum
const int 	bins_Pe       	= 12;
const double 	Pe_step    	= 0.5;
const double 	Pe_min     	= 2.;
const double 	Pe_max     	= Pe_min + Pe_step*bins_Pe;

//The = electron theta (in degrees)
const int 	bins_The      	= 30;
const double 	The_step   	= 1.;
const double 	The_min    	= 5.;
const double 	The_max    	= The_min + The_step*bins_The;

//Phe = electron phi (in degrees)
const int 	bins_Phe      	= 72;
const double	Phe_min    	= -180.;
const double 	Phe_max    	= +180.;

const int 	bins_W2       	= 9;
const double 	W2_step    	= 1.;
const double 	W2_min     	= 4.;
const double 	W2_max     	= W2_min + W2_step*bins_W2;


/*
const double Q2Bins[3] = {2,2.7,6};

const int bins_Q2 	= 2;
const int lowQ2_bins 	= 21;
const int highQ2_bins 	= 20;

const int bins_Xb	= 11;
const double Xb_step	= 0.05;
const double Xb_min	= 0.125;
const double Xb_max	= Xb_min + Xb_step*bins_Xb;

const double lowQ2_AsPtBins[21][4]	= 
	{ 
	  // lower slice 
	  {1.23, 1.25, 0.04, 0.08},
	  {1.22, 1.25, 0.08, 0.10},
	  {1.21, 1.25, 0.10, 0.13},
  	  
	  {1.25, 1.30, 0.01, 0.04},
	  {1.25, 1.30, 0.04, 0.06},
	  {1.25, 1.30, 0.06, 0.08},
	  {1.25, 1.30, 0.08, 0.10},
	  {1.25, 1.30, 0.10, 0.12},
	  {1.25, 1.30, 0.12, 0.15},

	  {1.30, 1.35, 0.02, 0.06},
	  {1.30, 1.35, 0.06, 0.09},
	  {1.30, 1.35, 0.09, 0.11},
	  {1.30, 1.35, 0.11, 0.14},

	  {1.35, 1.40, 0.03, 0.07},
	  {1.35, 1.40, 0.07, 0.10},
	  {1.35, 1.40, 0.10, 0.13},
	  {1.35, 1.40, 0.13, 0.16},

	  {1.40, 1.45, 0.04, 0.11},
	  {1.40, 1.45, 0.11, 0.17},
	
	  {1.45, 1.50, 0.05, 0.17},

	  {1.50, 1.60, 0.06, 0.17}

	};
	
const double highQ2_AsPtBins[20][4]	=
	{ 
	  // lower slice 
	  {1.23, 1.25, 0.04, 0.08},
	  {1.22, 1.25, 0.08, 0.10},
	  {1.21, 1.25, 0.10, 0.13},
  	  
	  {1.25, 1.30, 0.01, 0.04},
	  {1.25, 1.30, 0.04, 0.06},
	  {1.25, 1.30, 0.06, 0.08},
	  {1.25, 1.30, 0.08, 0.10},
	  {1.25, 1.30, 0.10, 0.12},
	  {1.25, 1.30, 0.12, 0.15},

	  {1.30, 1.35, 0.02, 0.06},
	  {1.30, 1.35, 0.06, 0.09},
	  {1.30, 1.35, 0.09, 0.11},
	  {1.30, 1.35, 0.11, 0.14},

	  {1.35, 1.40, 0.03, 0.07},
	  {1.35, 1.40, 0.07, 0.10},
	  {1.35, 1.40, 0.10, 0.13},
	  {1.35, 1.40, 0.13, 0.16},

	  {1.40, 1.45, 0.04, 0.11},
	  {1.40, 1.45, 0.11, 0.17},
	
	  {1.45, 1.50, 0.05, 0.17},

	};
*/

#endif
