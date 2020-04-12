#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_bucher1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[8] = pow(ASL_c, 2)*CYP3A4_ASLoOH_Vmax/(pow(CYP3A4_ASLoOH_Km1, 3)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) - ASL_c*CYP3A4_ASLoOH_Vmax/(pow(CYP3A4_ASLoOH_Km1, 2)*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            dwdp[9] = pow(ASL_c, 2)*CYP3A4_ASLpOH_Vmax/(pow(CYP3A4_ASLoOH_Km1, 2)*CYP3A4_ASLpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[10] = ASL_c*AS_c*CYP3A4_ASoOH_Vmax/(pow(CYP3A4_ASLoOH_Km1, 2)*CYP3A4_ASoOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[11] = ASL_c*AS_c*CYP3A4_ASpOH_Vmax/(pow(CYP3A4_ASLoOH_Km1, 2)*CYP3A4_ASpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            break;
        case 1:
            dwdp[8] = ASL_c/(CYP3A4_ASLoOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            break;
        case 2:
            dwdp[8] = pow(ASL_c, 2)*CYP3A4_ASLoOH_Vmax/(CYP3A4_ASLoOH_Km1*pow(CYP3A4_ASLpOH_Km1, 2)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[9] = pow(ASL_c, 2)*CYP3A4_ASLpOH_Vmax/(pow(CYP3A4_ASLpOH_Km1, 3)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) - ASL_c*CYP3A4_ASLpOH_Vmax/(pow(CYP3A4_ASLpOH_Km1, 2)*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            dwdp[10] = ASL_c*AS_c*CYP3A4_ASoOH_Vmax/(pow(CYP3A4_ASLpOH_Km1, 2)*CYP3A4_ASoOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[11] = ASL_c*AS_c*CYP3A4_ASpOH_Vmax/(pow(CYP3A4_ASLpOH_Km1, 2)*CYP3A4_ASpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            break;
        case 3:
            dwdp[9] = ASL_c/(CYP3A4_ASLpOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            break;
        case 4:
            dwdp[8] = ASL_c*AS_c*CYP3A4_ASLoOH_Vmax/(CYP3A4_ASLoOH_Km1*pow(CYP3A4_ASoOH_Km1, 2)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[9] = ASL_c*AS_c*CYP3A4_ASLpOH_Vmax/(CYP3A4_ASLpOH_Km1*pow(CYP3A4_ASoOH_Km1, 2)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[10] = pow(AS_c, 2)*CYP3A4_ASoOH_Vmax/(pow(CYP3A4_ASoOH_Km1, 3)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) - AS_c*CYP3A4_ASoOH_Vmax/(pow(CYP3A4_ASoOH_Km1, 2)*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            dwdp[11] = pow(AS_c, 2)*CYP3A4_ASpOH_Vmax/(pow(CYP3A4_ASoOH_Km1, 2)*CYP3A4_ASpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            break;
        case 5:
            dwdp[10] = AS_c/(CYP3A4_ASoOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            break;
        case 6:
            dwdp[8] = ASL_c*AS_c*CYP3A4_ASLoOH_Vmax/(CYP3A4_ASLoOH_Km1*pow(CYP3A4_ASpOH_Km1, 2)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[9] = ASL_c*AS_c*CYP3A4_ASLpOH_Vmax/(CYP3A4_ASLpOH_Km1*pow(CYP3A4_ASpOH_Km1, 2)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[10] = pow(AS_c, 2)*CYP3A4_ASoOH_Vmax/(CYP3A4_ASoOH_Km1*pow(CYP3A4_ASpOH_Km1, 2)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
            dwdp[11] = pow(AS_c, 2)*CYP3A4_ASpOH_Vmax/(pow(CYP3A4_ASpOH_Km1, 3)*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) - AS_c*CYP3A4_ASpOH_Vmax/(pow(CYP3A4_ASpOH_Km1, 2)*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            break;
        case 7:
            dwdp[11] = AS_c/(CYP3A4_ASpOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
            break;
        case 8:
            dwdp[13] = ASL_c;
            break;
        case 9:
            dwdp[14] = ASLoOH_c;
            break;
        case 10:
            dwdp[15] = ASLpOH_c;
            break;
        case 11:
            dwdp[12] = AS_c;
            break;
        case 12:
            dwdp[16] = ASoOH_c;
            break;
        case 13:
            dwdp[17] = ASpOH_c;
            break;
        case 14:
            dwdp[19] = ASL_m;
            break;
        case 15:
            dwdp[20] = ASLoOH_m;
            break;
        case 16:
            dwdp[21] = ASLpOH_m;
            break;
        case 17:
            dwdp[18] = AS_m;
            break;
        case 18:
            dwdp[22] = ASoOH_m;
            break;
        case 19:
            dwdp[23] = ASpOH_m;
            break;
        case 20:
            dwdp[0] = -ASL_b + ASL_c*(-fu_ASL + 1)/fu_ASL;
            dwdp[1] = -ASLoOH_b + ASLoOH_c*(-fu_ASL + 1)/fu_ASL;
            dwdp[2] = -ASLpOH_b + ASLpOH_c*(-fu_ASL + 1)/fu_ASL;
            dwdp[3] = -AS_b + AS_c*(-fu_AS + 1)/fu_AS;
            dwdp[4] = -ASoOH_b + ASoOH_c*(-fu_AS + 1)/fu_AS;
            dwdp[5] = -ASpOH_b + ASpOH_c*(-fu_AS + 1)/fu_AS;
            break;
        case 21:
            dwdp[28] = pow(AS_c, 3)*UGT1A3_AS_Vmax/(pow(UGT1A3_AS_KI1, 2)*pow(pow(AS_c, 2)/UGT1A3_AS_KI1 + AS_c + UGT1A3_AS_Km1, 2));
            break;
        case 22:
            dwdp[28] = -AS_c*UGT1A3_AS_Vmax/pow(pow(AS_c, 2)/UGT1A3_AS_KI1 + AS_c + UGT1A3_AS_Km1, 2);
            break;
        case 23:
            dwdp[28] = AS_c/(pow(AS_c, 2)/UGT1A3_AS_KI1 + AS_c + UGT1A3_AS_Km1);
            break;
        case 24:
            dwdp[3] = Prot_k1*(-AS_c/fu_AS - AS_c*(-fu_AS + 1)/pow(fu_AS, 2));
            dwdp[4] = Prot_k1*(-ASoOH_c/fu_AS - ASoOH_c*(-fu_AS + 1)/pow(fu_AS, 2));
            dwdp[5] = Prot_k1*(-ASpOH_c/fu_AS - ASpOH_c*(-fu_AS + 1)/pow(fu_AS, 2));
            break;
        case 25:
            dwdp[0] = Prot_k1*(-ASL_c/fu_ASL - ASL_c*(-fu_ASL + 1)/pow(fu_ASL, 2));
            dwdp[1] = Prot_k1*(-ASLoOH_c/fu_ASL - ASLoOH_c*(-fu_ASL + 1)/pow(fu_ASL, 2));
            dwdp[2] = Prot_k1*(-ASLpOH_c/fu_ASL - ASLpOH_c*(-fu_ASL + 1)/pow(fu_ASL, 2));
            break;
        case 26:
            dwdp[6] = ASLoOH_c;
            dwdp[7] = ASLpOH_c;
            dwdp[24] = ASL_c;
            break;
        case 27:
            dwdp[25] = ASL_m;
            dwdp[26] = ASLoOH_m;
            dwdp[27] = ASLpOH_m;
            break;
        case 28:
            dwdp[24] = ASL_c;
            break;
        case 29:
            dwdp[6] = ASLoOH_c;
            dwdp[7] = ASLpOH_c;
            break;
    }
}