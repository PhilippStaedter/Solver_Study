#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Eungdamrong2007(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[2] = -(2.6931270074281662e-6*CAPRI_cyt*Ca*Ca_bind_CAPRI_Kf - 0.0016212624584717584*CaCAPRI_cyt*Ca_bind_CAPRI_Kr)/pow(KMOLE, 2);
            dwdp[3] = -0.97599999999999998*IP3_degradation_kIP3deg*(0.0016611295681063099*IP3 - IP3_degradation_IP3_basal)/pow(KMOLE, 2);
            dwdp[5] = -(2.6931270074281662e-6*Ca*RasGRP_cyt*ca_bind_rasGRP_Kf - 0.0016212624584717584*CaRasGRP1_cyt*ca_bind_rasGRP_Kr)/pow(KMOLE, 2);
            dwdp[15] = -(2.6931270074281662e-6*Grb2_cyt*Sos_cyt*sos_grb2_binding_Kf - 0.0016212624584717584*SosGrb2_cyt*sos_grb2_binding_Kr)/pow(KMOLE, 2);
            dwdp[18] = -(-0.0016212624584717584*RasGDP_pal_cyt*RasPal_basal_GAP_Kr + 0.0016212624584717584*RasGTP_pal_cyt*RasPal_basal_GAP_Kf)/pow(KMOLE, 2);
            dwdp[19] = -(-0.0016212624584717584*RasGDP_depal_cyt*basal_cyt_depal_GEF_Kr + 0.0016212624584717584*RasGTP_depal_cyt*basal_cyt_depal_GEF_Kf)/pow(KMOLE, 2);
            dwdp[26] = -(2.6931270074281662e-6*Ca*buffer_cyt*calcium_buffer_Kf - 0.0016212624584717584*ca_buffer_cyt*calcium_buffer_Kr)/pow(KMOLE, 2);
            dwdp[30] = -(2.6931270074281662e-6*Ca*PLCe_cyt*ca_act_PLCe_Kf - 0.0016212624584717584*Ca_PLCe_cyt*ca_act_PLCe_Kr)/pow(KMOLE, 2);
            dwdp[39] = -1.2582642575688975e-7*pow(Ca, 2)*ER_erMembrane*flux1_vP*serca/(pow(KMOLE, 2)*(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2)));
            dwdp[40] = 5.2253499068475741e-11*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(pow(KMOLE, 2)*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
            dwdp[41] = 0.045600000000000002*ER_erMembrane*flux2_vL*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/pow(KMOLE, 2);
            break;
        case 1:
            dwdp[0] = 0.59999999999999998*PIP_PM*((t > tauPIP2syn_PIP2_synthesis) ? (
   exp(-(t - tauPIP2syn_PIP2_synthesis)/PIP2syndecay_PIP2_synthesis)
)
: (
   0
));
            break;
        case 2:
            dwdp[0] = 0.59999999999999998*PIP_PM*((t > tauPIP2syn_PIP2_synthesis) ? (
   kStimSynPIP2_PIP2_synthesis*exp(-(t - tauPIP2syn_PIP2_synthesis)/PIP2syndecay_PIP2_synthesis)/PIP2syndecay_PIP2_synthesis
)
: (
   0
));
            break;
        case 3:
            dwdp[0] = 0.59999999999999998*PIP_PM*((t > tauPIP2syn_PIP2_synthesis) ? (
   -kStimSynPIP2_PIP2_synthesis*(-t + tauPIP2syn_PIP2_synthesis)*exp(-(t - tauPIP2syn_PIP2_synthesis)/PIP2syndecay_PIP2_synthesis)/pow(PIP2syndecay_PIP2_synthesis, 2)
)
: (
   0
));
            break;
        case 4:
            dwdp[0] = 0.59999999999999998*PIP_PM*((PIP2_PM < PIP2_basal_PIP2_synthesis) ? (
   0.58099999999999996*kBasalSynPIP2_PIP2_synthesis*(1.0/PIP2_basal_PIP2_synthesis - (-PIP2_PM + PIP2_basal_PIP2_synthesis)/pow(PIP2_basal_PIP2_synthesis, 2))*exp((-PIP2_PM + PIP2_basal_PIP2_synthesis)/PIP2_basal_PIP2_synthesis)
)
: (
   0
));
            break;
        case 5:
            dwdp[0] = 0.59999999999999998*PIP_PM*((PIP2_PM < PIP2_basal_PIP2_synthesis) ? (
   0.58099999999999996*exp((-PIP2_PM + PIP2_basal_PIP2_synthesis)/PIP2_basal_PIP2_synthesis) - 0.58099999999999996
)
: (
   0
));
            break;
        case 6:
            dwdp[24] = 0.59999999999999998*PI_PM*((PIP_PM < PIP_basal_PIP_synthesis) ? (
   0.58099999999999996*kBasalSynPIP_PIP_synthesis*(1.0/PIP_basal_PIP_synthesis - (-PIP_PM + PIP_basal_PIP_synthesis)/pow(PIP_basal_PIP_synthesis, 2))*exp((-PIP_PM + PIP_basal_PIP_synthesis)/PIP_basal_PIP_synthesis)
)
: (
   0
));
            break;
        case 7:
            dwdp[24] = 0.59999999999999998*PI_PM*((PIP_PM < PIP_basal_PIP_synthesis) ? (
   0.58099999999999996*exp((-PIP_PM + PIP_basal_PIP_synthesis)/PIP_basal_PIP_synthesis) - 0.58099999999999996
)
: (
   0
));
            break;
        case 8:
            dwdp[24] = 0.59999999999999998*PI_PM*((t > tauPIPsyn_PIP_synthesis) ? (
   exp(-(t - tauPIPsyn_PIP_synthesis)/PIPsyndecay_PIP_synthesis)
)
: (
   0
));
            break;
        case 9:
            dwdp[24] = 0.59999999999999998*PI_PM*((t > tauPIPsyn_PIP_synthesis) ? (
   kStimSynPIP_PIP_synthesis*exp(-(t - tauPIPsyn_PIP_synthesis)/PIPsyndecay_PIP_synthesis)/PIPsyndecay_PIP_synthesis
)
: (
   0
));
            break;
        case 10:
            dwdp[24] = 0.59999999999999998*PI_PM*((t > tauPIPsyn_PIP_synthesis) ? (
   -kStimSynPIP_PIP_synthesis*(-t + tauPIPsyn_PIP_synthesis)*exp(-(t - tauPIPsyn_PIP_synthesis)/PIPsyndecay_PIP_synthesis)/pow(PIPsyndecay_PIP_synthesis, 2)
)
: (
   0
));
            break;
        case 11:
            dwdp[34] = 7.5747508305647741e-5*Ca*Rinh - 0.045600000000000002*RinhCa*dinh_reaction2;
            break;
        case 12:
            dwdp[34] = -0.045600000000000002*Kon_reaction2*RinhCa;
            break;
        case 13:
            dwdp[36] = -45.600000000000001*RactCa;
            break;
        case 15:
            dwdp[1] = -0.12*Activated_EGFR_PM*Shc_PM/pow(Shc_PM + Shc_phosphorylation_Km, 2);
            break;
        case 17:
            dwdp[2] = -0.0016212624584717584*CaCAPRI_cyt/KMOLE;
            break;
        case 18:
            dwdp[2] = 2.6931270074281662e-6*CAPRI_cyt*Ca/KMOLE;
            break;
        case 19:
            dwdp[3] = -0.97599999999999998*IP3_degradation_kIP3deg/KMOLE;
            break;
        case 20:
            dwdp[3] = 0.97599999999999998*(0.0016611295681063099*IP3 - IP3_degradation_IP3_basal)/KMOLE;
            break;
        case 21:
            dwdp[4] = -0.048000000000000001*RasGTP_Golgi_GM;
            break;
        case 22:
            dwdp[4] = 7.9734219269102881e-5*RasGTP_depal_cyt;
            break;
        case 24:
            dwdp[5] = -0.0016212624584717584*CaRasGRP1_cyt/KMOLE;
            break;
        case 25:
            dwdp[5] = 2.6931270074281662e-6*Ca*RasGRP_cyt/KMOLE;
            break;
        case 26:
            dwdp[6] = -0.0024000000000000002*RasGDP_Golgi_GM*RasGRP_DAG_GM/pow(RasGDP_Golgi_GM + RasGRP_DAG_GEF_Km, 2);
            break;
        case 28:
            dwdp[7] = -0.59999999999999998*RasGTP_PM;
            break;
        case 29:
            dwdp[7] = 0.00099667774086378597*RasGTP_pal_cyt;
            break;
        case 31:
            dwdp[8] = -0.59999999999999998*PLC_PM;
            break;
        case 32:
            dwdp[8] = 0.59999999999999998*PLC_act_PM;
            break;
        case 34:
            dwdp[9] = -0.59999999999999998*RasGDP_PM;
            break;
        case 35:
            dwdp[9] = 0.59999999999999998*RasGTP_PM;
            break;
        case 37:
            dwdp[10] = -0.59999999999999998*CaCAPRI_PM_PM;
            break;
        case 38:
            dwdp[10] = 0.00099667774086378597*CaCAPRI_cyt;
            break;
        case 40:
            dwdp[11] = -0.048000000000000001*RasGRP_DAG_GM;
            break;
        case 41:
            dwdp[11] = 7.9734219269102881e-5*DAG_GM_GM*RasGRP_cyt;
            break;
        case 43:
            dwdp[12] = -7.9734219269102881e-5*RasGDP_pal_cyt;
            break;
        case 44:
            dwdp[12] = 0.048000000000000001*RasGDP_Golgi_GM;
            break;
        case 46:
            dwdp[13] = -6.0*CaCAPRI_PM_PM*RasGTP_PM/pow(CAPRI_GAP_Km + RasGTP_PM, 2);
            break;
        case 48:
            dwdp[14] = -0.59999999999999998*RasGDP_PM;
            break;
        case 49:
            dwdp[14] = 0.00099667774086378597*RasGDP_pal_cyt;
            break;
        case 51:
            dwdp[15] = -0.0016212624584717584*SosGrb2_cyt/KMOLE;
            break;
        case 52:
            dwdp[15] = 2.6931270074281662e-6*Grb2_cyt*Sos_cyt/KMOLE;
            break;
        case 53:
            dwdp[16] = -0.048000000000000001*RasGDP_Golgi_GM;
            break;
        case 54:
            dwdp[16] = 7.9734219269102881e-5*RasGDP_depal_cyt;
            break;
        case 56:
            dwdp[17] = -7.9734219269102881e-5*RasGTP_pal_cyt;
            break;
        case 57:
            dwdp[17] = 0.048000000000000001*RasGTP_Golgi_GM;
            break;
        case 59:
            dwdp[18] = -0.0016212624584717584*RasGDP_pal_cyt/KMOLE;
            break;
        case 60:
            dwdp[18] = 0.0016212624584717584*RasGTP_pal_cyt/KMOLE;
            break;
        case 61:
            dwdp[19] = -0.0016212624584717584*RasGDP_depal_cyt/KMOLE;
            break;
        case 62:
            dwdp[19] = 0.0016212624584717584*RasGTP_depal_cyt/KMOLE;
            break;
        case 63:
            dwdp[20] = 0.048000000000000001*PIP2_GM_GM*Ras_CaPLCe_GM;
            break;
        case 65:
            dwdp[21] = -0.00048000000000000001*Ca_RasGRP_GM_GM*RasGDP_Golgi_GM/pow(CaRasGRP_act_RasGM_Km + RasGDP_Golgi_GM, 2);
            break;
        case 67:
            dwdp[22] = 0.59999999999999998*PIP2_PM*PLC_act_PM;
            break;
        case 69:
            dwdp[23] = -0.59999999999999998*SGS_PM;
            break;
        case 70:
            dwdp[23] = 0.00099667774086378597*Shc_star_PM*SosGrb2_cyt;
            break;
        case 74:
            dwdp[25] = 0.59999999999999998*Activated_EGFR_PM;
            break;
        case 76:
            dwdp[26] = -0.0016212624584717584*ca_buffer_cyt/KMOLE;
            break;
        case 77:
            dwdp[26] = 2.6931270074281662e-6*Ca*buffer_cyt/KMOLE;
            break;
        case 78:
            dwdp[27] = 0.048000000000000001*RasGTP_Golgi_GM/(RasGM_basal_GAP_Km + RasGTP_Golgi_GM);
            break;
        case 79:
            dwdp[27] = -0.048000000000000001*RasGM_basal_GAP_Vmax*RasGTP_Golgi_GM/pow(RasGM_basal_GAP_Km + RasGTP_Golgi_GM, 2);
            break;
        case 81:
            dwdp[28] = -0.59999999999999998*Shc_PM;
            break;
        case 82:
            dwdp[28] = 0.59999999999999998*Shc_star_PM;
            break;
        case 84:
            dwdp[29] = -0.59999999999999998*Activated_EGFR_PM;
            break;
        case 85:
            dwdp[29] = 0.00099667774086378597*EGFR_PM*EGF_EC;
            break;
        case 87:
            dwdp[30] = -0.0016212624584717584*Ca_PLCe_cyt/KMOLE;
            break;
        case 88:
            dwdp[30] = 2.6931270074281662e-6*Ca*PLCe_cyt/KMOLE;
            break;
        case 89:
            dwdp[31] = -0.048000000000000001*Ras_CaPLCe_GM;
            break;
        case 90:
            dwdp[31] = 7.9734219269102881e-5*Ca_PLCe_cyt*RasGTP_Golgi_GM;
            break;
        case 92:
            dwdp[32] = -0.00099667774086378597*RasGDP_depal_cyt;
            break;
        case 93:
            dwdp[32] = 0.59999999999999998*RasGDP_PM;
            break;
        case 95:
            dwdp[33] = -0.048000000000000001*Ca_RasGRP_GM_GM;
            break;
        case 96:
            dwdp[33] = 7.9734219269102881e-5*CaRasGRP1_cyt;
            break;
        case 99:
            dwdp[35] = -0.17999999999999999*Activated_EGFR_PM*PLC_PM/pow(EGF_act_PLCgamma_Km + PLC_PM, 2);
            break;
        case 101:
            dwdp[36] = 7.5747508305647741e-5*Ca*Ract;
            break;
        case 104:
            dwdp[37] = 0.048000000000000001*DAG_GM_GM;
            break;
        case 106:
            dwdp[38] = -0.012*RasGDP_PM*SGS_PM/pow(RasGDP_PM + Sos_act_RasPM_Km, 2);
            break;
        case 108:
            dwdp[39] = -2.5165285151377951e-7*pow(Ca, 2)*ER_erMembrane*flux1_kP*flux1_vP*serca/(KMOLE*pow(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2), 2));
            break;
        case 109:
            dwdp[39] = 1.2582642575688975e-7*pow(Ca, 2)*ER_erMembrane*serca/(KMOLE*(2.7593514420370556e-6*pow(Ca, 2) + pow(flux1_kP, 2)));
            break;
        case 111:
            dwdp[40] = -5.2253499068475741e-11*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 3)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
            break;
        case 112:
            dwdp[40] = 1.5676049720542721e-10*ER_erMembrane*pow(IP3, 3)*pow(RactCa, 3)*pow(Rinh, 3)*flux0_singleChanFlux*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/(KMOLE*pow(0.0016611295681063099*IP3 + flux0_dI, 4)*pow(Ract + RactCa, 2)*pow(Rinh + RinhCa, 3));
            break;
        case 114:
            dwdp[41] = -0.045600000000000002*ER_erMembrane*(-0.0016611295681063099*Ca + 0.0016611295681063099*Ca_ER)/KMOLE;
            break;
        case 116:
            dwdp[42] = -0.00099667774086378597*RasGTP_depal_cyt;
            break;
        case 117:
            dwdp[42] = 0.59999999999999998*RasGTP_PM;
            break;
    }
}