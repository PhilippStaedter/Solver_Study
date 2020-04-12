#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Eungdamrong2007(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 0.59999999999999998*PIP_PM*(((PIP2_PM < PIP2_basal_PIP2_synthesis) ? (
   0.58099999999999996*kBasalSynPIP2_PIP2_synthesis*(exp((-PIP2_PM + PIP2_basal_PIP2_synthesis)/PIP2_basal_PIP2_synthesis) - 1)
)
: (
   0
)) + ((t > tauPIP2syn_PIP2_synthesis) ? (
   kStimSynPIP2_PIP2_synthesis*exp(-(t - tauPIP2syn_PIP2_synthesis)/PIP2syndecay_PIP2_synthesis)
)
: (
   0
)));
    w[1] = 0.12*Activated_EGFR_PM*Shc_PM/(Shc_PM + Shc_phosphorylation_Km);
    w[2] = 0.97599999999999998*(2.7593514420370556e-6*CAPRI_cyt*Ca*Ca_bind_CAPRI_Kf - 0.0016611295681063099*CaCAPRI_cyt*Ca_bind_CAPRI_Kr)/KMOLE;
    w[3] = 0.97599999999999998*IP3_degradation_kIP3deg*(0.0016611295681063099*IP3 - IP3_degradation_IP3_basal)/KMOLE;
    w[4] = -0.048000000000000001*RasGTP_Golgi_GM*RasGTP_depal_translocate_Kr + 7.9734219269102881e-5*RasGTP_depal_cyt*RasGTP_depal_translocate_Kf;
    w[5] = 0.97599999999999998*(2.7593514420370556e-6*Ca*RasGRP_cyt*ca_bind_rasGRP_Kf - 0.0016611295681063099*CaRasGRP1_cyt*ca_bind_rasGRP_Kr)/KMOLE;
    w[6] = 0.0024000000000000002*RasGDP_Golgi_GM*RasGRP_DAG_GM/(RasGDP_Golgi_GM + RasGRP_DAG_GEF_Km);
    w[7] = -0.59999999999999998*RasGTP_PM*rasGTP_pal_translocation_Kr + 0.00099667774086378597*RasGTP_pal_cyt*rasGTP_pal_translocation_Kf;
    w[8] = -0.59999999999999998*PLC_PM*PLCg_dephos_Kr + 0.59999999999999998*PLC_act_PM*PLCg_dephos_Kf;
    w[9] = -0.59999999999999998*RasGDP_PM*basal_GAP_Kr + 0.59999999999999998*RasGTP_PM*basal_GAP_Kf;
    w[10] = 0.00099667774086378597*CAPRI_translocation_Kf*CaCAPRI_cyt - 0.59999999999999998*CAPRI_translocation_Kr*CaCAPRI_PM_PM;
    w[11] = 7.9734219269102881e-5*DAG_GM_GM*RasGRP_cyt*reaction5_Kf - 0.048000000000000001*RasGRP_DAG_GM*reaction5_Kr;
    w[12] = 0.048000000000000001*RasGDP_Golgi_GM*RasGDP_pal_Kf - 7.9734219269102881e-5*RasGDP_pal_Kr*RasGDP_pal_cyt;
    w[13] = 6.0*CaCAPRI_PM_PM*RasGTP_PM/(CAPRI_GAP_Km + RasGTP_PM);
    w[14] = -0.59999999999999998*RasGDP_PM*RasGDPpal_translocation_Kr + 0.00099667774086378597*RasGDP_pal_cyt*RasGDPpal_translocation_Kf;
    w[15] = 0.97599999999999998*(2.7593514420370556e-6*Grb2_cyt*Sos_cyt*sos_grb2_binding_Kf - 0.0016611295681063099*SosGrb2_cyt*sos_grb2_binding_Kr)/KMOLE;
    w[16] = -0.048000000000000001*RasGDP_Golgi_GM*RasGDP_depal_translocate_Kr + 7.9734219269102881e-5*RasGDP_depal_cyt*RasGDP_depal_translocate_Kf;
    w[17] = 0.048000000000000001*RasGTP_Golgi_GM*Ras_GTP_palm1_Kf - 7.9734219269102881e-5*RasGTP_pal_cyt*Ras_GTP_palm1_Kr;
    w[18] = 0.97599999999999998*(-0.0016611295681063099*RasGDP_pal_cyt*RasPal_basal_GAP_Kr + 0.0016611295681063099*RasGTP_pal_cyt*RasPal_basal_GAP_Kf)/KMOLE;
    w[19] = 0.97599999999999998*(-0.0016611295681063099*RasGDP_depal_cyt*basal_cyt_depal_GEF_Kr + 0.0016611295681063099*RasGTP_depal_cyt*basal_cyt_depal_GEF_Kf)/KMOLE;
    w[20] = 0.048000000000000001*PIP2_GM_GM*Ras_CaPLCe_GM*caPLCe_gen_DAG_kact;
    w[21] = 0.00048000000000000001*Ca_RasGRP_GM_GM*RasGDP_Golgi_GM/(CaRasGRP_act_RasGM_Km + RasGDP_Golgi_GM);
    w[22] = 0.59999999999999998*PIP2_PM*PIP2_hydrolysis_k_PIP2hyd*PLC_act_PM;
    w[23] = -0.59999999999999998*SGS_PM*Sos_activation_Kr + 0.00099667774086378597*Shc_star_PM*SosGrb2_cyt*Sos_activation_Kf;
    w[24] = 0.59999999999999998*PI_PM*(((PIP_PM < PIP_basal_PIP_synthesis) ? (
   0.58099999999999996*kBasalSynPIP_PIP_synthesis*(exp((-PIP_PM + PIP_basal_PIP_synthesis)/PIP_basal_PIP_synthesis) - 1)
)
: (
   0
)) + ((t > tauPIPsyn_PIP_synthesis) ? (
   kStimSynPIP_PIP_synthesis*exp(-(t - tauPIPsyn_PIP_synthesis)/PIPsyndecay_PIP_synthesis)
)
: (
   0
)));
    w[25] = 0.59999999999999998*Activated_EGFR_PM*EGF_internalization_Kf;
    w[26] = 0.97599999999999998*(2.7593514420370556e-6*Ca*buffer_cyt*calcium_buffer_Kf - 0.0016611295681063099*ca_buffer_cyt*calcium_buffer_Kr)/KMOLE;
    w[27] = 0.048000000000000001*RasGM_basal_GAP_Vmax*RasGTP_Golgi_GM/(RasGM_basal_GAP_Km + RasGTP_Golgi_GM);
    w[28] = -0.59999999999999998*Shc_PM*reaction0_Kr + 0.59999999999999998*Shc_star_PM*reaction0_Kf;
    w[29] = -0.59999999999999998*Activated_EGFR_PM*EGFR_binding_Kr + 0.00099667774086378597*EGFR_PM*EGFR_binding_Kf*EGF_EC;
    w[30] = 0.97599999999999998*(2.7593514420370556e-6*Ca*PLCe_cyt*ca_act_PLCe_Kf - 0.0016611295681063099*Ca_PLCe_cyt*ca_act_PLCe_Kr)/KMOLE;
    w[31] = 7.9734219269102881e-5*Ca_PLCe_cyt*RasGTP_Golgi_GM*ras_act_PLCe_Kf - 0.048000000000000001*Ras_CaPLCe_GM*ras_act_PLCe_Kr;
    w[32] = 0.59999999999999998*RasGDP_PM*RasGDP_depal2_Kf - 0.00099667774086378597*RasGDP_depal2_Kr*RasGDP_depal_cyt;
    w[33] = 7.9734219269102881e-5*CaRasGRP1_cyt*CaRasGRP_translocation_Kf - 0.048000000000000001*CaRasGRP_translocation_Kr*Ca_RasGRP_GM_GM;
    w[34] = 7.5747508305647741e-5*Ca*Kon_reaction2*Rinh - 0.045600000000000002*Kon_reaction2*RinhCa*dinh_reaction2;
    w[35] = 0.17999999999999999*Activated_EGFR_PM*PLC_PM/(EGF_act_PLCgamma_Km + PLC_PM);
    w[36] = 7.5747508305647741e-5*Ca*Ca_binds_IP3R_Kf*Ract - 45.600000000000001*RactCa*dact_Ca_binds_IP3R;
    w[37] = 0.048000000000000001*DAG_GM_GM*reaction7_Kf;
    w[38] = 0.012*RasGDP_PM*SGS_PM/(RasGDP_PM + Sos_act_RasPM_Km);
    w[39] = 1.2582642575688975e-7*pow(Ca, 2)*ER_erMembrane*flux1_vP*serca/(KMOLE*(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2)));
    w[40] = -5.2253499068475741e-11*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
    w[41] = -0.045600000000000002*ER_erMembrane*flux2_vL*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/KMOLE;
    w[42] = 0.59999999999999998*RasGTP_PM*Ras_PM_depal1_Kf - 0.00099667774086378597*RasGTP_depal_cyt*Ras_PM_depal1_Kr;
}