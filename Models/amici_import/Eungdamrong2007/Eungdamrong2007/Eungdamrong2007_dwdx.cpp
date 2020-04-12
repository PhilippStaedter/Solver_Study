#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Eungdamrong2007(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -0.048000000000000001*RasGTP_depal_translocate_Kr;
    dwdx[1] = 0.048000000000000001*Ras_GTP_palm1_Kf;
    dwdx[2] = -0.048000000000000001*RasGM_basal_GAP_Vmax*RasGTP_Golgi_GM/pow(RasGM_basal_GAP_Km + RasGTP_Golgi_GM, 2) + 0.048000000000000001*RasGM_basal_GAP_Vmax/(RasGM_basal_GAP_Km + RasGTP_Golgi_GM);
    dwdx[3] = 7.9734219269102881e-5*Ca_PLCe_cyt*ras_act_PLCe_Kf;
    dwdx[4] = 0.00099667774086378597*EGFR_PM*EGFR_binding_Kf;
    dwdx[5] = 2.6931270074281662e-6*Ca*Ca_bind_CAPRI_Kf/KMOLE;
    dwdx[6] = 1.2582642575688975e-7*pow(Ca, 2)*ER_erMembrane*flux1_vP/(KMOLE*(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2)));
    dwdx[7] = 0.59999999999999998*((PIP2_PM < PIP2_basal_PIP2_synthesis) ? (
   0.58099999999999996*kBasalSynPIP2_PIP2_synthesis*(exp((-PIP2_PM + PIP2_basal_PIP2_synthesis)/PIP2_basal_PIP2_synthesis) - 1)
)
: (
   0
)) + 0.59999999999999998*((t > tauPIP2syn_PIP2_synthesis) ? (
   kStimSynPIP2_PIP2_synthesis*exp(-(t - tauPIP2syn_PIP2_synthesis)/PIP2syndecay_PIP2_synthesis)
)
: (
   0
));
    dwdx[8] = 0.59999999999999998*PI_PM*((PIP_PM < PIP_basal_PIP_synthesis) ? (
   -0.58099999999999996*kBasalSynPIP_PIP_synthesis*exp((-PIP_PM + PIP_basal_PIP_synthesis)/PIP_basal_PIP_synthesis)/PIP_basal_PIP_synthesis
)
: (
   0
));
    dwdx[9] = 0.59999999999999998*((PIP_PM < PIP_basal_PIP_synthesis) ? (
   0.58099999999999996*kBasalSynPIP_PIP_synthesis*(exp((-PIP_PM + PIP_basal_PIP_synthesis)/PIP_basal_PIP_synthesis) - 1)
)
: (
   0
)) + 0.59999999999999998*((t > tauPIPsyn_PIP_synthesis) ? (
   kStimSynPIP_PIP_synthesis*exp(-(t - tauPIPsyn_PIP_synthesis)/PIPsyndecay_PIP_synthesis)
)
: (
   0
));
    dwdx[10] = -0.12*Activated_EGFR_PM*Shc_PM/pow(Shc_PM + Shc_phosphorylation_Km, 2) + 0.12*Activated_EGFR_PM/(Shc_PM + Shc_phosphorylation_Km);
    dwdx[11] = -0.59999999999999998*reaction0_Kr;
    dwdx[12] = -0.59999999999999998*CAPRI_translocation_Kr;
    dwdx[13] = 6.0*RasGTP_PM/(CAPRI_GAP_Km + RasGTP_PM);
    dwdx[14] = -45.600000000000001*dact_Ca_binds_IP3R;
    dwdx[15] = 1.0450699813695148e-10*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 3)*pow(Rinh + RinhCa, 3)) - 1.5676049720542721e-10*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 2)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
    dwdx[16] = 0.00099667774086378597*SosGrb2_cyt*Sos_activation_Kf;
    dwdx[17] = 0.59999999999999998*reaction0_Kf;
    dwdx[18] = 0.00099667774086378597*EGFR_binding_Kf*EGF_EC;
    dwdx[19] = 0.59999999999999998*PLCg_dephos_Kf;
    dwdx[20] = 0.59999999999999998*PIP2_PM*PIP2_hydrolysis_k_PIP2hyd;
    dwdx[21] = 0.00099667774086378597*rasGTP_pal_translocation_Kf;
    dwdx[22] = -7.9734219269102881e-5*Ras_GTP_palm1_Kr;
    dwdx[23] = 0.0016212624584717584*RasPal_basal_GAP_Kf/KMOLE;
    dwdx[24] = -0.59999999999999998*PLCg_dephos_Kr;
    dwdx[25] = -0.17999999999999999*Activated_EGFR_PM*PLC_PM/pow(EGF_act_PLCgamma_Km + PLC_PM, 2) + 0.17999999999999999*Activated_EGFR_PM/(EGF_act_PLCgamma_Km + PLC_PM);
    dwdx[26] = 0.59999999999999998*PIP_PM*((PIP2_PM < PIP2_basal_PIP2_synthesis) ? (
   -0.58099999999999996*kBasalSynPIP2_PIP2_synthesis*exp((-PIP2_PM + PIP2_basal_PIP2_synthesis)/PIP2_basal_PIP2_synthesis)/PIP2_basal_PIP2_synthesis
)
: (
   0
));
    dwdx[27] = 0.59999999999999998*PIP2_hydrolysis_k_PIP2hyd*PLC_act_PM;
    dwdx[28] = 0.12*Shc_PM/(Shc_PM + Shc_phosphorylation_Km);
    dwdx[29] = 0.59999999999999998*EGF_internalization_Kf;
    dwdx[30] = -0.59999999999999998*EGFR_binding_Kr;
    dwdx[31] = 0.17999999999999999*PLC_PM/(EGF_act_PLCgamma_Km + PLC_PM);
    dwdx[32] = 2.6931270074281662e-6*CAPRI_cyt*Ca_bind_CAPRI_Kf/KMOLE;
    dwdx[33] = 2.6931270074281662e-6*RasGRP_cyt*ca_bind_rasGRP_Kf/KMOLE;
    dwdx[34] = 2.6931270074281662e-6*buffer_cyt*calcium_buffer_Kf/KMOLE;
    dwdx[35] = 2.6931270074281662e-6*PLCe_cyt*ca_act_PLCe_Kf/KMOLE;
    dwdx[36] = 7.5747508305647741e-5*Kon_reaction2*Rinh;
    dwdx[37] = 7.5747508305647741e-5*Ca_binds_IP3R_Kf*Ract;
    dwdx[38] = -6.9439865871728446e-13*pow(Ca, 3)*ER_erMembrane*flux1_vP*serca/(KMOLE*pow(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2), 2)) + 2.5165285151377951e-7*Ca*ER_erMembrane*flux1_vP*serca/(KMOLE*(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2)));
    dwdx[39] = 8.6799832339660571e-14*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
    dwdx[40] = 7.5747508305647741e-5*ER_erMembrane*flux2_vL/KMOLE;
    dwdx[41] = 7.5747508305647741e-5*Ca*Ca_binds_IP3R_Kf;
    dwdx[42] = 1.0450699813695148e-10*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 3)*pow(Rinh + RinhCa, 3));
    dwdx[43] = 7.5747508305647741e-5*Ca*Kon_reaction2;
    dwdx[44] = 1.5676049720542721e-10*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 4)) - 1.5676049720542721e-10*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 2)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
    dwdx[45] = -0.045600000000000002*Kon_reaction2*dinh_reaction2;
    dwdx[46] = 1.5676049720542721e-10*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 4));
    dwdx[47] = 0.0016212624584717584*IP3_degradation_kIP3deg/KMOLE;
    dwdx[48] = 2.6039949701898172e-13*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 4)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3)) - 1.5676049720542721e-10*ER_erMembrane*pow(IP3, 2)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
    dwdx[49] = -0.0024000000000000002*RasGDP_Golgi_GM*RasGRP_DAG_GM/pow(RasGDP_Golgi_GM + RasGRP_DAG_GEF_Km, 2) + 0.0024000000000000002*RasGRP_DAG_GM/(RasGDP_Golgi_GM + RasGRP_DAG_GEF_Km);
    dwdx[50] = 0.048000000000000001*RasGDP_pal_Kf;
    dwdx[51] = -0.048000000000000001*RasGDP_depal_translocate_Kr;
    dwdx[52] = -0.00048000000000000001*Ca_RasGRP_GM_GM*RasGDP_Golgi_GM/pow(CaRasGRP_act_RasGM_Km + RasGDP_Golgi_GM, 2) + 0.00048000000000000001*Ca_RasGRP_GM_GM/(CaRasGRP_act_RasGM_Km + RasGDP_Golgi_GM);
    dwdx[53] = 0.00048000000000000001*RasGDP_Golgi_GM/(CaRasGRP_act_RasGM_Km + RasGDP_Golgi_GM);
    dwdx[54] = -0.048000000000000001*CaRasGRP_translocation_Kr;
    dwdx[55] = 7.9734219269102881e-5*RasGRP_cyt*reaction5_Kf;
    dwdx[56] = 0.048000000000000001*reaction7_Kf;
    dwdx[57] = 0.0024000000000000002*RasGDP_Golgi_GM/(RasGDP_Golgi_GM + RasGRP_DAG_GEF_Km);
    dwdx[58] = -0.048000000000000001*reaction5_Kr;
    dwdx[59] = -0.0016212624584717584*Ca_bind_CAPRI_Kr/KMOLE;
    dwdx[60] = 0.00099667774086378597*CAPRI_translocation_Kf;
    dwdx[61] = 7.9734219269102881e-5*RasGTP_depal_translocate_Kf;
    dwdx[62] = 0.0016212624584717584*basal_cyt_depal_GEF_Kf/KMOLE;
    dwdx[63] = -0.00099667774086378597*Ras_PM_depal1_Kr;
    dwdx[64] = 7.9734219269102881e-5*RasGDP_depal_translocate_Kf;
    dwdx[65] = -0.0016212624584717584*basal_cyt_depal_GEF_Kr/KMOLE;
    dwdx[66] = -0.00099667774086378597*RasGDP_depal2_Kr;
    dwdx[67] = -7.9734219269102881e-5*RasGDP_pal_Kr;
    dwdx[68] = 0.00099667774086378597*RasGDPpal_translocation_Kf;
    dwdx[69] = -0.0016212624584717584*RasPal_basal_GAP_Kr/KMOLE;
    dwdx[70] = -0.0016212624584717584*ca_act_PLCe_Kr/KMOLE;
    dwdx[71] = 7.9734219269102881e-5*RasGTP_Golgi_GM*ras_act_PLCe_Kf;
    dwdx[72] = 0.048000000000000001*PIP2_GM_GM*caPLCe_gen_DAG_kact;
    dwdx[73] = -0.048000000000000001*ras_act_PLCe_Kr;
    dwdx[74] = 0.048000000000000001*Ras_CaPLCe_GM*caPLCe_gen_DAG_kact;
    dwdx[75] = 1.2582642575688975e-7*pow(Ca, 2)*flux1_vP*serca/(KMOLE*(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2)));
    dwdx[76] = -5.2253499068475741e-11*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
    dwdx[77] = -0.045600000000000002*flux2_vL*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/KMOLE;
    dwdx[78] = -8.6799832339660571e-14*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
    dwdx[79] = -7.5747508305647741e-5*ER_erMembrane*flux2_vL/KMOLE;
    dwdx[80] = 2.6931270074281662e-6*Grb2_cyt*sos_grb2_binding_Kf/KMOLE;
    dwdx[81] = 2.6931270074281662e-6*Sos_cyt*sos_grb2_binding_Kf/KMOLE;
    dwdx[82] = 2.6931270074281662e-6*Ca*ca_act_PLCe_Kf/KMOLE;
    dwdx[83] = 2.6931270074281662e-6*Ca*calcium_buffer_Kf/KMOLE;
    dwdx[84] = -0.0016212624584717584*calcium_buffer_Kr/KMOLE;
    dwdx[85] = -0.0016212624584717584*sos_grb2_binding_Kr/KMOLE;
    dwdx[86] = 0.00099667774086378597*Shc_star_PM*Sos_activation_Kf;
    dwdx[87] = -0.59999999999999998*Sos_activation_Kr;
    dwdx[88] = 0.012*RasGDP_PM/(RasGDP_PM + Sos_act_RasPM_Km);
    dwdx[89] = -0.59999999999999998*rasGTP_pal_translocation_Kr;
    dwdx[90] = 0.59999999999999998*basal_GAP_Kf;
    dwdx[91] = -6.0*CaCAPRI_PM_PM*RasGTP_PM/pow(CAPRI_GAP_Km + RasGTP_PM, 2) + 6.0*CaCAPRI_PM_PM/(CAPRI_GAP_Km + RasGTP_PM);
    dwdx[92] = 0.59999999999999998*Ras_PM_depal1_Kf;
    dwdx[93] = -0.59999999999999998*basal_GAP_Kr;
    dwdx[94] = -0.59999999999999998*RasGDPpal_translocation_Kr;
    dwdx[95] = 0.59999999999999998*RasGDP_depal2_Kf;
    dwdx[96] = -0.012*RasGDP_PM*SGS_PM/pow(RasGDP_PM + Sos_act_RasPM_Km, 2) + 0.012*SGS_PM/(RasGDP_PM + Sos_act_RasPM_Km);
    dwdx[97] = 2.6931270074281662e-6*Ca*ca_bind_rasGRP_Kf/KMOLE;
    dwdx[98] = 7.9734219269102881e-5*DAG_GM_GM*reaction5_Kf;
    dwdx[99] = -0.0016212624584717584*ca_bind_rasGRP_Kr/KMOLE;
    dwdx[100] = 7.9734219269102881e-5*CaRasGRP_translocation_Kf;
}