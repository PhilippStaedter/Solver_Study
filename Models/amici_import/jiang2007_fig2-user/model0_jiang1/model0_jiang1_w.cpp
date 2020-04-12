#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_jiang1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*GAP*flow;
    w[1] = 1.0*flow*(-GLC + (1.0/1000.0)*GLCflow_Glc_F);
    w[2] = 1.0*LAC*flow;
    w[3] = -1.0*pow(ADP_cyt, 2)*hidden_1_k9b + 1.0*AMP*ATP_cyt*hidden_1_k9f;
    w[4] = 1.0*ATP_cyt*GLC*v1_V1/((ATP_cyt + v1_K1ATP)*(GLC + v1_K1GLC));
    w[5] = 1.0*Acetyl_CoA*OXA*v10_V*v10_v10_CS/(Acetyl_CoA*OXA + Acetyl_CoA*v10_Kb + OXA*v10_Ka + v10_Kia*v10_Kib);
    w[6] = 1.0*v11_v11_ACO*(Cit*v11_KcF*v11_Kp - IsoCit*v11_KcR*v11_Ks)/(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks);
    w[7] = 1.0*v12_KcF*v12_v12_IDHa*(ADP*IsoCit*v12_b + pow(IsoCit, 2))/(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f);
    w[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC/(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB);
    w[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)/(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB);
    w[10] = 1.0*v16_KcF*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)/(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2);
    w[11] = 1.0*v17_v17_FM*(Fum*v17_KcF*v17_Kp - Mal*v17_KcR*v17_Ks)/(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks);
    w[12] = 1.0*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))/(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1);
    w[13] = 1.0*ATP_cyt*pow(F6P, 2)*v2_V2/((ATP_cyt + v2_K2ATP)*(pow(F6P, 2) + v2_K2*(1 + pow(ATP_cyt, 2)*v2_k2/pow(AMP, 2))));
    w[14] = 1.0*v20_KcF*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)/(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq);
    w[15] = 1.0*v21_KcF*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)/(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2);
    w[16] = 1.0*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))/(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1);
    w[17] = 1.0*v24_KcF*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)/(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq);
    w[18] = 1.0*Cytc3p*QH2*v25_KcF*v25_v25_Complex_III/(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB);
    w[19] = 1.0*Cytc2p*v26_KcF*v26_v26_Complex_IV/(Cytc2p + v26_Ks);
    w[20] = 1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*v27_Kia*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
    w[21] = 1.0*ADP*v28_V*v28_v28_Complex_V/(pow(ADP, 2)/v28_Ki + ADP + v28_Km);
    w[22] = 1.0*v29_v29_ACO*(Cit_cyt*v29_KcF*v29_Kp - IsoCitcyt*v29_KcR*v29_Ks)/(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks);
    w[23] = 1.0*FBP*v3_k3f - 1.0*pow(GAP, 2)*v3_k3b;
    w[24] = 1.0*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))/(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1);
    w[25] = 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2);
    w[26] = 1.0*v32_KcF*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq);
    w[27] = 1.0*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))/(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1);
    w[28] = 1.0*v34_KcF*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)/(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq);
    w[29] = 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2);
    w[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip);
    w[31] = 1.0*G3P*v37_V*v37_v37_GUT2P/(G3P + v37_K);
    w[32] = 1.0*NADH_cyt*v38_V*v38_v38_GUT2P/(NADH_cyt + v38_K);
    w[33] = 1.0*Mal_cyt*NADP_cyt*v39_Kcat*v39_v39_MDH/((Mal_cyt + v39_Kmal)*(NADP_cyt + v39_Knadp));
    w[34] = 1.0*GAP*NAD*v4_V4/((GAP + v4_K4GAP)*(NAD + v4_K4NAD));
    w[35] = 1.0*ADP_cyt*v40_V*v40_v40_AAC/(ADP_cyt + v40_K);
    w[36] = 1.0*v41_v41_IDHc*(-CO2*NADPH_cyt*OG_cyt/(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123) + IsoCitcyt*NADP_cyt/(IsoCitcyt*NADP_cyt*v41_phi0 + IsoCitcyt*v41_phi2 + NADP_cyt*v41_phi1 + v41_phi12));
    w[37] = 1.0*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))/(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1);
    w[38] = 1.0*ATP*v43_V*v43_v43_AAC/(ATP + v43_K);
    w[39] = 1.0*Mal*v44_Kcat*v44_v44_MDH/(Mal + v44_Km);
    w[40] = 1.0*ADP_cyt*DPG*v5_k5f - 1.0*ATP_cyt*PEP*v5_k5b;
    w[41] = 1.0*ADP_cyt*PEP*v6_V6/((ADP_cyt + v6_K6ADP)*(PEP + v6_K6PEP));
    w[42] = -1.0*LAC*NAD*v7_k8b + 1.0*NADH_cyt*PYR_cyt*v7_k8f;
    w[43] = 1.0*PYR_cyt*v8_V*v8_v8_PYC/(PYR_cyt + v8_K);
    w[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC/(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB);
}