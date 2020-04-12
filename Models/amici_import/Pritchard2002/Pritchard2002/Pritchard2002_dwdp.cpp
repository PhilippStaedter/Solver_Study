#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Pritchard2002(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -GLCi*GLCo*HXT_Vmax_1*(-GLCi + GLCo)/(pow(HXT_Kglc_1, 3)*pow(GLCi*GLCo*HXT_Ki_1/pow(HXT_Kglc_1, 2) + 1 + (GLCi + GLCo)/HXT_Kglc_1, 2));
            break;
        case 1:
            dwdp[0] = HXT_Vmax_1*(-GLCi + GLCo)*(2*GLCi*GLCo*HXT_Ki_1/pow(HXT_Kglc_1, 3) + (GLCi + GLCo)/pow(HXT_Kglc_1, 2))/(HXT_Kglc_1*pow(GLCi*GLCo*HXT_Ki_1/pow(HXT_Kglc_1, 2) + 1 + (GLCi + GLCo)/HXT_Kglc_1, 2)) - HXT_Vmax_1*(-GLCi + GLCo)/(pow(HXT_Kglc_1, 2)*(GLCi*GLCo*HXT_Ki_1/pow(HXT_Kglc_1, 2) + 1 + (GLCi + GLCo)/HXT_Kglc_1));
            break;
        case 2:
            dwdp[0] = (-GLCi + GLCo)/(HXT_Kglc_1*(GLCi*GLCo*HXT_Ki_1/pow(HXT_Kglc_1, 2) + 1 + (GLCi + GLCo)/HXT_Kglc_1));
            break;
        case 3:
            dwdp[1] = 1.0*ADP*HK_Vmax_2*(-ADP*G6P/(HK_Katp_2*HK_Keq_2*HK_Kglc_2) + ATP*GLCi/(HK_Katp_2*HK_Kglc_2))/(pow(HK_Kadp_2, 2)*pow(ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1, 2)*(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1));
            break;
        case 4:
            dwdp[1] = 1.0*G6P*HK_Vmax_2*(-ADP*G6P/(HK_Katp_2*HK_Keq_2*HK_Kglc_2) + ATP*GLCi/(HK_Katp_2*HK_Kglc_2))/(pow(HK_Kg6p_2, 2)*(ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1)*pow(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1, 2));
            break;
        case 5:
            dwdp[1] = 1.0*ADP*G6P*HK_Vmax_2/(HK_Katp_2*pow(HK_Keq_2, 2)*HK_Kglc_2*(ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1)*(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1));
            break;
        case 6:
            dwdp[1] = 1.0*ATP*HK_Vmax_2*(-ADP*G6P/(HK_Katp_2*HK_Keq_2*HK_Kglc_2) + ATP*GLCi/(HK_Katp_2*HK_Kglc_2))/(pow(HK_Katp_2, 2)*pow(ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1, 2)*(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1)) + 1.0*HK_Vmax_2*(ADP*G6P/(pow(HK_Katp_2, 2)*HK_Keq_2*HK_Kglc_2) - ATP*GLCi/(pow(HK_Katp_2, 2)*HK_Kglc_2))/((ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1)*(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1));
            break;
        case 7:
            dwdp[1] = 1.0*GLCi*HK_Vmax_2*(-ADP*G6P/(HK_Katp_2*HK_Keq_2*HK_Kglc_2) + ATP*GLCi/(HK_Katp_2*HK_Kglc_2))/(pow(HK_Kglc_2, 2)*(ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1)*pow(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1, 2)) + 1.0*HK_Vmax_2*(ADP*G6P/(HK_Katp_2*HK_Keq_2*pow(HK_Kglc_2, 2)) - ATP*GLCi/(HK_Katp_2*pow(HK_Kglc_2, 2)))/((ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1)*(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1));
            break;
        case 8:
            dwdp[1] = 1.0*(-ADP*G6P/(HK_Katp_2*HK_Keq_2*HK_Kglc_2) + ATP*GLCi/(HK_Katp_2*HK_Kglc_2))/((ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1)*(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1));
            break;
        case 9:
            dwdp[2] = 1.0*F6P*PGI_Vmax_3*(-F6P/(PGI_Keq_3*PGI_Kg6p_3) + G6P/PGI_Kg6p_3)/(pow(PGI_Kf6p_3, 2)*pow(F6P/PGI_Kf6p_3 + G6P/PGI_Kg6p_3 + 1, 2));
            break;
        case 10:
            dwdp[2] = 1.0*F6P*PGI_Vmax_3/(pow(PGI_Keq_3, 2)*PGI_Kg6p_3*(F6P/PGI_Kf6p_3 + G6P/PGI_Kg6p_3 + 1));
            break;
        case 11:
            dwdp[2] = 1.0*G6P*PGI_Vmax_3*(-F6P/(PGI_Keq_3*PGI_Kg6p_3) + G6P/PGI_Kg6p_3)/(pow(PGI_Kg6p_3, 2)*pow(F6P/PGI_Kf6p_3 + G6P/PGI_Kg6p_3 + 1, 2)) + 1.0*PGI_Vmax_3*(F6P/(PGI_Keq_3*pow(PGI_Kg6p_3, 2)) - G6P/pow(PGI_Kg6p_3, 2))/(F6P/PGI_Kf6p_3 + G6P/PGI_Kg6p_3 + 1);
            break;
        case 12:
            dwdp[2] = 1.0*(-F6P/(PGI_Keq_3*PGI_Kg6p_3) + G6P/PGI_Kg6p_3)/(F6P/PGI_Kf6p_3 + G6P/PGI_Kg6p_3 + 1);
            break;
        case 13:
            dwdp[3] = -2.0*pow(ATP, 2)*F6P*PFK_L0_4*PFK_Vmax_4*PFK_gR_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*(ATP*PFK_Catp_4/PFK_Katp_4 + 1)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(pow(PFK_Katp_4, 2)*PFK_Kf6p_4*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2));
            break;
        case 14:
            dwdp[3] = 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(2*F16bP*PFK_Cf16_4*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1)/(pow(PFK_Kf16_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) - 2*F16bP*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(PFK_Kf16_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 3)))*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2));
            break;
        case 15:
            dwdp[3] = -2.0*ATP*F16bP*F6P*PFK_L0_4*PFK_Vmax_4*PFK_gR_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1)*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf16_4*PFK_Kf6p_4*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2));
            break;
        case 16:
            dwdp[3] = 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(2*F26bP*PFK_Cf26_4*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1)/(pow(PFK_Kf26_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) - 2*F26bP*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(PFK_Kf26_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 3)))*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2));
            break;
        case 17:
            dwdp[3] = -2.0*ATP*F26bP*F6P*PFK_L0_4*PFK_Vmax_4*PFK_gR_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1)*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf26_4*PFK_Kf6p_4*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2));
            break;
        case 18:
            dwdp[3] = 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(2*AMP*PFK_Camp_4*PFK_L0_4*(AMP*PFK_Camp_4/PFK_Kamp_4 + 1)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(PFK_Kamp_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) - 2*AMP*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(PFK_Kamp_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 3)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)))*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2));
            break;
        case 19:
            dwdp[3] = -2.0*AMP*ATP*F6P*PFK_L0_4*PFK_Vmax_4*PFK_gR_4*(AMP*PFK_Camp_4/PFK_Kamp_4 + 1)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Kamp_4*PFK_Katp_4*PFK_Kf6p_4*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2));
            break;
        case 20:
            dwdp[3] = 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(2*ATP*PFK_Ciatp_4*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(PFK_Kiatp_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) - 2*ATP*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(PFK_Kiatp_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 3)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)))*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2));
            break;
        case 21:
            dwdp[3] = -2.0*pow(ATP, 2)*F6P*PFK_L0_4*PFK_Vmax_4*PFK_gR_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*PFK_Kiatp_4*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2));
            break;
        case 22:
            dwdp[3] = -1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2));
            break;
        case 23:
            dwdp[3] = 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(-ATP*F6P*PFK_gR_4/(pow(PFK_Katp_4, 2)*PFK_Kf6p_4) - ATP/pow(PFK_Katp_4, 2))/(PFK_Katp_4*PFK_Kf6p_4*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2))) + 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(2*ATP*PFK_Catp_4*PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*(ATP*PFK_Catp_4/PFK_Katp_4 + 1)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(PFK_Katp_4, 2)*pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) - (-2*ATP*F6P*PFK_gR_4/(pow(PFK_Katp_4, 2)*PFK_Kf6p_4) - 2*ATP/pow(PFK_Katp_4, 2))*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1))*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)) - 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(pow(PFK_Katp_4, 2)*PFK_Kf6p_4*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2)));
            break;
        case 24:
            dwdp[3] = -1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(-2*ATP*F6P*PFK_gR_4/(PFK_Katp_4*pow(PFK_Kf6p_4, 2)) - 2*F6P/pow(PFK_Kf6p_4, 2))*pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2)/(PFK_Katp_4*PFK_Kf6p_4*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)) + 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(-ATP*F6P*PFK_gR_4/(PFK_Katp_4*pow(PFK_Kf6p_4, 2)) - F6P/pow(PFK_Kf6p_4, 2))/(PFK_Katp_4*PFK_Kf6p_4*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2))) - 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*pow(PFK_Kf6p_4, 2)*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2)));
            break;
        case 25:
            dwdp[3] = 1.0*pow(ATP, 2)*pow(F6P, 2)*PFK_Vmax_4*PFK_gR_4/(pow(PFK_Katp_4, 2)*pow(PFK_Kf6p_4, 2)*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2))) - 2.0*pow(ATP, 2)*pow(F6P, 2)*PFK_Vmax_4*PFK_gR_4*pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2)/(pow(PFK_Katp_4, 2)*pow(PFK_Kf6p_4, 2)*pow(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2), 2)) + 1.0*ATP*F6P*PFK_Vmax_4*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2)));
            break;
        case 26:
            dwdp[3] = 1.0*ATP*F6P*PFK_gR_4*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2)));
            break;
        case 27:
            dwdp[4] = 1.0*ALD_Vmax_5*F16bP*GAP*(F16bP/ALD_Kf16bp_5 - DHAP*GAP/(ALD_Keq_5*ALD_Kf16bp_5))/(ALD_Kf16bp_5*pow(ALD_Kigap_5, 2)*pow(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5), 2));
            break;
        case 28:
            dwdp[4] = 1.0*ALD_Vmax_5*(F16bP/ALD_Kf16bp_5 - DHAP*GAP/(ALD_Keq_5*ALD_Kf16bp_5))*(GAP/pow(ALD_Kgap_5, 2) + DHAP*GAP/(ALD_Kdhap_5*pow(ALD_Kgap_5, 2)))/pow(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5), 2);
            break;
        case 29:
            dwdp[4] = 1.0*ALD_Vmax_5*(DHAP/pow(ALD_Kdhap_5, 2) + DHAP*GAP/(pow(ALD_Kdhap_5, 2)*ALD_Kgap_5))*(F16bP/ALD_Kf16bp_5 - DHAP*GAP/(ALD_Keq_5*ALD_Kf16bp_5))/pow(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5), 2);
            break;
        case 30:
            dwdp[4] = 1.0*ALD_Vmax_5*DHAP*GAP/(pow(ALD_Keq_5, 2)*ALD_Kf16bp_5*(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5)));
            break;
        case 31:
            dwdp[4] = 1.0*ALD_Vmax_5*(-F16bP/pow(ALD_Kf16bp_5, 2) + DHAP*GAP/(ALD_Keq_5*pow(ALD_Kf16bp_5, 2)))/(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5)) + 1.0*ALD_Vmax_5*(F16bP/pow(ALD_Kf16bp_5, 2) + F16bP*GAP/(pow(ALD_Kf16bp_5, 2)*ALD_Kigap_5))*(F16bP/ALD_Kf16bp_5 - DHAP*GAP/(ALD_Keq_5*ALD_Kf16bp_5))/pow(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5), 2);
            break;
        case 32:
            dwdp[4] = 1.0*(F16bP/ALD_Kf16bp_5 - DHAP*GAP/(ALD_Keq_5*ALD_Kf16bp_5))/(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5));
            break;
        case 33:
            dwdp[5] = -1.0*GAP;
            break;
        case 34:
            dwdp[5] = 1.0*DHAP;
            break;
        case 35:
            dwdp[6] = 1.0*BPG*GAPDH_C_7*GAPDH_Vmaxr_7*NADH/(GAPDH_Kbpg_7*pow(GAPDH_Knadh_7, 2)*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1)) + 1.0*GAPDH_C_7*NADH*(-BPG*GAPDH_Vmaxr_7*NADH/(GAPDH_Kbpg_7*GAPDH_Knadh_7) + GAP*GAPDH_Vmaxf_7*NAD/(GAPDH_Kgap_7*GAPDH_Knad_7))/(pow(GAPDH_Knadh_7, 2)*pow(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7, 2)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
            break;
        case 36:
            dwdp[6] = 1.0*BPG*GAPDH_C_7*(-BPG*GAPDH_Vmaxr_7*NADH/(GAPDH_Kbpg_7*GAPDH_Knadh_7) + GAP*GAPDH_Vmaxf_7*NAD/(GAPDH_Kgap_7*GAPDH_Knad_7))/(pow(GAPDH_Kbpg_7, 2)*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*pow(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1, 2)) + 1.0*BPG*GAPDH_C_7*GAPDH_Vmaxr_7*NADH/(pow(GAPDH_Kbpg_7, 2)*GAPDH_Knadh_7*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
            break;
        case 37:
            dwdp[6] = -1.0*BPG*GAPDH_C_7*NADH/(GAPDH_Kbpg_7*GAPDH_Knadh_7*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
            break;
        case 38:
            dwdp[6] = -1.0*GAP*GAPDH_C_7*GAPDH_Vmaxf_7*NAD/(GAPDH_Kgap_7*pow(GAPDH_Knad_7, 2)*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1)) + 1.0*GAPDH_C_7*NAD*(-BPG*GAPDH_Vmaxr_7*NADH/(GAPDH_Kbpg_7*GAPDH_Knadh_7) + GAP*GAPDH_Vmaxf_7*NAD/(GAPDH_Kgap_7*GAPDH_Knad_7))/(pow(GAPDH_Knad_7, 2)*pow(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7, 2)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
            break;
        case 39:
            dwdp[6] = 1.0*GAP*GAPDH_C_7*(-BPG*GAPDH_Vmaxr_7*NADH/(GAPDH_Kbpg_7*GAPDH_Knadh_7) + GAP*GAPDH_Vmaxf_7*NAD/(GAPDH_Kgap_7*GAPDH_Knad_7))/(pow(GAPDH_Kgap_7, 2)*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*pow(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1, 2)) - 1.0*GAP*GAPDH_C_7*GAPDH_Vmaxf_7*NAD/(pow(GAPDH_Kgap_7, 2)*GAPDH_Knad_7*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
            break;
        case 40:
            dwdp[6] = 1.0*GAP*GAPDH_C_7*NAD/(GAPDH_Kgap_7*GAPDH_Knad_7*(1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
            break;
        case 41:
            dwdp[6] = 1.0*(-BPG*GAPDH_Vmaxr_7*NADH/(GAPDH_Kbpg_7*GAPDH_Knadh_7) + GAP*GAPDH_Vmaxf_7*NAD/(GAPDH_Kgap_7*GAPDH_Knad_7))/((1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
            break;
        case 42:
            dwdp[7] = 1.0*ADP*PGK_Vmax_8*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(pow(PGK_Kadp_8, 2)*PGK_Katp_8*PGK_Kp3g_8*pow(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1, 2)*(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1));
            break;
        case 43:
            dwdp[7] = 1.0*BPG*PGK_Vmax_8*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(PGK_Katp_8*pow(PGK_Kbpg_8, 2)*PGK_Kp3g_8*(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1)*pow(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1, 2));
            break;
        case 44:
            dwdp[7] = 1.0*ATP*PGK_Vmax_8*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(pow(PGK_Katp_8, 3)*PGK_Kp3g_8*pow(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1, 2)*(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1)) - 1.0*PGK_Vmax_8*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(pow(PGK_Katp_8, 2)*PGK_Kp3g_8*(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1)*(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1));
            break;
        case 45:
            dwdp[7] = 1.0*P3G*PGK_Vmax_8*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(PGK_Katp_8*pow(PGK_Kp3g_8, 3)*(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1)*pow(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1, 2)) - 1.0*PGK_Vmax_8*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(PGK_Katp_8*pow(PGK_Kp3g_8, 2)*(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1)*(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1));
            break;
        case 46:
            dwdp[7] = 1.0*ADP*BPG*PGK_Vmax_8/(PGK_Katp_8*PGK_Kp3g_8*(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1)*(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1));
            break;
        case 47:
            dwdp[7] = 1.0*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(PGK_Katp_8*PGK_Kp3g_8*(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1)*(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1));
            break;
        case 48:
            dwdp[8] = 1.0*P2G*PGM_Vmax_9*(-P2G/(PGM_Keq_9*PGM_Kp3g_9) + P3G/PGM_Kp3g_9)/(pow(PGM_Kp2g_9, 2)*pow(P2G/PGM_Kp2g_9 + P3G/PGM_Kp3g_9 + 1, 2));
            break;
        case 49:
            dwdp[8] = 1.0*P2G*PGM_Vmax_9/(pow(PGM_Keq_9, 2)*PGM_Kp3g_9*(P2G/PGM_Kp2g_9 + P3G/PGM_Kp3g_9 + 1));
            break;
        case 50:
            dwdp[8] = 1.0*P3G*PGM_Vmax_9*(-P2G/(PGM_Keq_9*PGM_Kp3g_9) + P3G/PGM_Kp3g_9)/(pow(PGM_Kp3g_9, 2)*pow(P2G/PGM_Kp2g_9 + P3G/PGM_Kp3g_9 + 1, 2)) + 1.0*PGM_Vmax_9*(P2G/(PGM_Keq_9*pow(PGM_Kp3g_9, 2)) - P3G/pow(PGM_Kp3g_9, 2))/(P2G/PGM_Kp2g_9 + P3G/PGM_Kp3g_9 + 1);
            break;
        case 51:
            dwdp[8] = 1.0*(-P2G/(PGM_Keq_9*PGM_Kp3g_9) + P3G/PGM_Kp3g_9)/(P2G/PGM_Kp2g_9 + P3G/PGM_Kp3g_9 + 1);
            break;
        case 52:
            dwdp[9] = 1.0*ENO_Vmax_10*PEP*(P2G/ENO_Kp2g_10 - PEP/(ENO_Keq_10*ENO_Kp2g_10))/(pow(ENO_Kpep_10, 2)*pow(1 + PEP/ENO_Kpep_10 + P2G/ENO_Kp2g_10, 2));
            break;
        case 53:
            dwdp[9] = 1.0*ENO_Vmax_10*PEP/(pow(ENO_Keq_10, 2)*ENO_Kp2g_10*(1 + PEP/ENO_Kpep_10 + P2G/ENO_Kp2g_10));
            break;
        case 54:
            dwdp[9] = 1.0*ENO_Vmax_10*(-P2G/pow(ENO_Kp2g_10, 2) + PEP/(ENO_Keq_10*pow(ENO_Kp2g_10, 2)))/(1 + PEP/ENO_Kpep_10 + P2G/ENO_Kp2g_10) + 1.0*ENO_Vmax_10*P2G*(P2G/ENO_Kp2g_10 - PEP/(ENO_Keq_10*ENO_Kp2g_10))/(pow(ENO_Kp2g_10, 2)*pow(1 + PEP/ENO_Kpep_10 + P2G/ENO_Kp2g_10, 2));
            break;
        case 55:
            dwdp[9] = 1.0*(P2G/ENO_Kp2g_10 - PEP/(ENO_Keq_10*ENO_Kp2g_10))/(1 + PEP/ENO_Kpep_10 + P2G/ENO_Kp2g_10);
            break;
        case 56:
            dwdp[10] = 1.0*ATP*PYK_Vmax_11*(ADP*PEP/(PYK_Kadp_11*PYK_Kpep_11) - ATP*PYR/(PYK_Kadp_11*PYK_Keq_11*PYK_Kpep_11))/(pow(PYK_Katp_11, 2)*pow(ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1, 2)*(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11));
            break;
        case 57:
            dwdp[10] = 1.0*PYK_Vmax_11*PYR*(ADP*PEP/(PYK_Kadp_11*PYK_Kpep_11) - ATP*PYR/(PYK_Kadp_11*PYK_Keq_11*PYK_Kpep_11))/(pow(PYK_Kpyr_11, 2)*(ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1)*pow(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11, 2));
            break;
        case 58:
            dwdp[10] = 1.0*ATP*PYK_Vmax_11*PYR/(PYK_Kadp_11*pow(PYK_Keq_11, 2)*PYK_Kpep_11*(ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1)*(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11));
            break;
        case 59:
            dwdp[10] = 1.0*ADP*PYK_Vmax_11*(ADP*PEP/(PYK_Kadp_11*PYK_Kpep_11) - ATP*PYR/(PYK_Kadp_11*PYK_Keq_11*PYK_Kpep_11))/(pow(PYK_Kadp_11, 2)*pow(ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1, 2)*(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11)) + 1.0*PYK_Vmax_11*(-ADP*PEP/(pow(PYK_Kadp_11, 2)*PYK_Kpep_11) + ATP*PYR/(pow(PYK_Kadp_11, 2)*PYK_Keq_11*PYK_Kpep_11))/((ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1)*(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11));
            break;
        case 60:
            dwdp[10] = 1.0*PEP*PYK_Vmax_11*(ADP*PEP/(PYK_Kadp_11*PYK_Kpep_11) - ATP*PYR/(PYK_Kadp_11*PYK_Keq_11*PYK_Kpep_11))/(pow(PYK_Kpep_11, 2)*(ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1)*pow(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11, 2)) + 1.0*PYK_Vmax_11*(-ADP*PEP/(PYK_Kadp_11*pow(PYK_Kpep_11, 2)) + ATP*PYR/(PYK_Kadp_11*PYK_Keq_11*pow(PYK_Kpep_11, 2)))/((ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1)*(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11));
            break;
        case 61:
            dwdp[10] = 1.0*(ADP*PEP/(PYK_Kadp_11*PYK_Kpep_11) - ATP*PYR/(PYK_Kadp_11*PYK_Keq_11*PYK_Kpep_11))/((ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1)*(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11));
            break;
        case 62:
            dwdp[11] = -1.0*PDC_Vmax_12*pow(PYR/PDC_Kpyr_12, 2*PDC_nH_12)*log(PYR/PDC_Kpyr_12)/pow(pow(PYR/PDC_Kpyr_12, PDC_nH_12) + 1, 2) + 1.0*PDC_Vmax_12*pow(PYR/PDC_Kpyr_12, PDC_nH_12)*log(PYR/PDC_Kpyr_12)/(pow(PYR/PDC_Kpyr_12, PDC_nH_12) + 1);
            break;
        case 63:
            dwdp[11] = 1.0*PDC_Vmax_12*PDC_nH_12*pow(PYR/PDC_Kpyr_12, 2*PDC_nH_12)/(PDC_Kpyr_12*pow(pow(PYR/PDC_Kpyr_12, PDC_nH_12) + 1, 2)) - 1.0*PDC_Vmax_12*PDC_nH_12*pow(PYR/PDC_Kpyr_12, PDC_nH_12)/(PDC_Kpyr_12*(pow(PYR/PDC_Kpyr_12, PDC_nH_12) + 1));
            break;
        case 64:
            dwdp[11] = 1.0*pow(PYR/PDC_Kpyr_12, PDC_nH_12)/(pow(PYR/PDC_Kpyr_12, PDC_nH_12) + 1);
            break;
        case 65:
            dwdp[12] = 1.0*ADH_Vmax_13*AcAld*EtOH*NADH*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))/(ADH_Kacald_13*pow(ADH_Kietoh_13, 2)*ADH_Kinadh_13*pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2));
            break;
        case 66:
            dwdp[12] = 1.0*ADH_Vmax_13*AcAld*EtOH*NAD*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))/(ADH_Ketoh_13*pow(ADH_Kiacald_13, 2)*ADH_Kinad_13*pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2));
            break;
        case 67:
            dwdp[12] = 1.0*ADH_Vmax_13*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))*(ADH_Knadh_13*AcAld/(pow(ADH_Kacald_13, 2)*ADH_Kinadh_13) + AcAld*NADH/(pow(ADH_Kacald_13, 2)*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(pow(ADH_Kacald_13, 2)*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(pow(ADH_Kacald_13, 2)*ADH_Kietoh_13*ADH_Kinadh_13))/pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2);
            break;
        case 68:
            dwdp[12] = 1.0*ADH_Vmax_13*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))*(NADH/pow(ADH_Kinadh_13, 2) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*pow(ADH_Kinadh_13, 2)) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*pow(ADH_Kinadh_13, 2)) + AcAld*NADH/(ADH_Kacald_13*pow(ADH_Kinadh_13, 2)) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*pow(ADH_Kinadh_13, 2)) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*pow(ADH_Kinadh_13, 2)))/pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2);
            break;
        case 69:
            dwdp[12] = 1.0*ADH_Vmax_13*(-AcAld/(ADH_Kacald_13*ADH_Kinadh_13) - AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13))*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))/pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2);
            break;
        case 70:
            dwdp[12] = 1.0*ADH_Vmax_13*(-EtOH/(ADH_Ketoh_13*ADH_Kinad_13) - EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13))*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))/pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2);
            break;
        case 71:
            dwdp[12] = 1.0*ADH_Vmax_13*AcAld*NADH/(pow(ADH_Keq_13, 2)*ADH_Ketoh_13*ADH_Kinad_13*(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13)));
            break;
        case 72:
            dwdp[12] = 1.0*ADH_Vmax_13*(-EtOH*NAD/(ADH_Ketoh_13*pow(ADH_Kinad_13, 2)) + AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*pow(ADH_Kinad_13, 2)))/(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13)) + 1.0*ADH_Vmax_13*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))*(NAD/pow(ADH_Kinad_13, 2) + ADH_Knad_13*EtOH/(ADH_Ketoh_13*pow(ADH_Kinad_13, 2)) + EtOH*NAD/(ADH_Ketoh_13*pow(ADH_Kinad_13, 2)) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*pow(ADH_Kinad_13, 2)*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*pow(ADH_Kinad_13, 2)) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*pow(ADH_Kinad_13, 2)*ADH_Kinadh_13))/pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2);
            break;
        case 73:
            dwdp[12] = 1.0*ADH_Vmax_13*(-EtOH*NAD/(pow(ADH_Ketoh_13, 2)*ADH_Kinad_13) + AcAld*NADH/(ADH_Keq_13*pow(ADH_Ketoh_13, 2)*ADH_Kinad_13))/(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13)) + 1.0*ADH_Vmax_13*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))*(ADH_Knad_13*EtOH/(pow(ADH_Ketoh_13, 2)*ADH_Kinad_13) + EtOH*NAD/(pow(ADH_Ketoh_13, 2)*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(pow(ADH_Ketoh_13, 2)*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(pow(ADH_Ketoh_13, 2)*ADH_Kiacald_13*ADH_Kinad_13))/pow(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13), 2);
            break;
        case 74:
            dwdp[12] = 1.0*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))/(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13));
            break;
        case 75:
            dwdp[13] = 1.0*ATP;
            break;
        case 76:
            dwdp[14] = -1.0*AMP*ATP;
            break;
        case 77:
            dwdp[14] = 1.0*pow(ADP, 2);
            break;
        case 78:
            dwdp[15] = 1.0*G3PDH_Vmax_16*NAD*(DHAP*NADH/(G3PDH_Kdhap_16*G3PDH_Knadh_16) - Glycerol*NAD/(G3PDH_Kdhap_16*G3PDH_Keq_16*G3PDH_Knadh_16))/(pow(G3PDH_Knad_16, 2)*pow(1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16, 2)*(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16));
            break;
        case 79:
            dwdp[15] = 1.0*G3PDH_Vmax_16*Glycerol*(DHAP*NADH/(G3PDH_Kdhap_16*G3PDH_Knadh_16) - Glycerol*NAD/(G3PDH_Kdhap_16*G3PDH_Keq_16*G3PDH_Knadh_16))/(pow(G3PDH_Kglycerol_16, 2)*(1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16)*pow(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16, 2));
            break;
        case 80:
            dwdp[15] = 1.0*G3PDH_Vmax_16*Glycerol*NAD/(G3PDH_Kdhap_16*pow(G3PDH_Keq_16, 2)*G3PDH_Knadh_16*(1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16)*(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16));
            break;
        case 81:
            dwdp[15] = 1.0*G3PDH_Vmax_16*(-DHAP*NADH/(G3PDH_Kdhap_16*pow(G3PDH_Knadh_16, 2)) + Glycerol*NAD/(G3PDH_Kdhap_16*G3PDH_Keq_16*pow(G3PDH_Knadh_16, 2)))/((1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16)*(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16)) + 1.0*G3PDH_Vmax_16*NADH*(DHAP*NADH/(G3PDH_Kdhap_16*G3PDH_Knadh_16) - Glycerol*NAD/(G3PDH_Kdhap_16*G3PDH_Keq_16*G3PDH_Knadh_16))/(pow(G3PDH_Knadh_16, 2)*pow(1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16, 2)*(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16));
            break;
        case 82:
            dwdp[15] = 1.0*DHAP*G3PDH_Vmax_16*(DHAP*NADH/(G3PDH_Kdhap_16*G3PDH_Knadh_16) - Glycerol*NAD/(G3PDH_Kdhap_16*G3PDH_Keq_16*G3PDH_Knadh_16))/(pow(G3PDH_Kdhap_16, 2)*(1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16)*pow(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16, 2)) + 1.0*G3PDH_Vmax_16*(-DHAP*NADH/(pow(G3PDH_Kdhap_16, 2)*G3PDH_Knadh_16) + Glycerol*NAD/(pow(G3PDH_Kdhap_16, 2)*G3PDH_Keq_16*G3PDH_Knadh_16))/((1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16)*(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16));
            break;
        case 83:
            dwdp[15] = 1.0*(DHAP*NADH/(G3PDH_Kdhap_16*G3PDH_Knadh_16) - Glycerol*NAD/(G3PDH_Kdhap_16*G3PDH_Keq_16*G3PDH_Knadh_16))/((1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16)*(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16));
            break;
        case 84:
            dwdp[16] = 1.0;
            break;
        case 85:
            dwdp[17] = 1.0;
            break;
        case 86:
            dwdp[18] = 1.0*AcAld;
            break;
    }
}