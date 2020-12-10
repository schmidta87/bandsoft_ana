#ifndef __KINEMATIC_CUTS_H__
#define __KINEMATIC_CUTS_H__

// Electron cuts
const int 	ECUT_PID = 11;
const int 	ECUT_charge = -1;
const double 	ECUT_EoP_min = 0.17;
const double 	ECUT_EoP_max = 0.3;
const double 	ECUT_Epcal_min = 0.07;
const double 	ECUT_V_min = 15;
const double 	ECUT_W_min = 15;
const double 	ECUT_vtx_min = -8;
const double 	ECUT_vtx_max =  3;
const double 	ECUT_pE_min = 3;
const double 	ECUT_pE_max = 10.6;
const double 	ECUT_Q2_min = 2;
const double 	ECUT_Q2_max = 10;
const double 	ECUT_W2_min = 2*2;
// Neutron PID cuts
const int	NCUT_goodneutron = 1;
const int	NCUT_leadindex = -1;
const int	NCUT_status = 0;
const double	SimAdcToMeVee = 1E4;
const double 	DataAdcToMeVee = 2300;
const double 	NCUT_Edep = 5;
const double 	NCUT_THETANQ_min = acos(-0.8);
const double	NCUT_THETANQ_max = acos(-1);
// Neutron signal region cuts
const double	NCUT_Tofabove0	= 0;
const double 	NCUT_Pn_min = 0.3;
const double 	NCUT_Pn_max = 0.7;
const double 	NCUT_Wp_min = 1.8;
const double 	NCUT_Wp_max = 4.5;
const double 	NCUT_As_min = 1.2;
const double	NCUT_As_max = 1.8;

#endif
