#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_bucher1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -Prot_k1;
    dwdx[1] = Prot_k1*(-fu_ASL + 1)/fu_ASL;
    dwdx[2] = ASL_c*CYP3A4_ASLoOH_Vmax*(-1/CYP3A4_ASLpOH_Km1 - 1/CYP3A4_ASLoOH_Km1)/(CYP3A4_ASLoOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) + CYP3A4_ASLoOH_Vmax/(CYP3A4_ASLoOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    dwdx[3] = ASL_c*CYP3A4_ASLpOH_Vmax*(-1/CYP3A4_ASLpOH_Km1 - 1/CYP3A4_ASLoOH_Km1)/(CYP3A4_ASLpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) + CYP3A4_ASLpOH_Vmax/(CYP3A4_ASLpOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    dwdx[4] = AS_c*CYP3A4_ASoOH_Vmax*(-1/CYP3A4_ASLpOH_Km1 - 1/CYP3A4_ASLoOH_Km1)/(CYP3A4_ASoOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
    dwdx[5] = AS_c*CYP3A4_ASpOH_Vmax*(-1/CYP3A4_ASLpOH_Km1 - 1/CYP3A4_ASLoOH_Km1)/(CYP3A4_ASpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
    dwdx[6] = Export_ASL_k;
    dwdx[7] = k_CR_ASL_c + k_PON_ASL_c;
    dwdx[8] = Import_ASL_k;
    dwdx[9] = k_CR_ASL_m;
    dwdx[10] = -Prot_k1;
    dwdx[11] = Prot_k1*(-fu_ASL + 1)/fu_ASL;
    dwdx[12] = k_CR_ASL_c + k_PON_OH_c;
    dwdx[13] = Export_ASLoOH_k;
    dwdx[14] = Import_ASLoOH_k;
    dwdx[15] = k_CR_ASL_m;
    dwdx[16] = -Prot_k1;
    dwdx[17] = Prot_k1*(-fu_ASL + 1)/fu_ASL;
    dwdx[18] = k_CR_ASL_c + k_PON_OH_c;
    dwdx[19] = Export_ASLpOH_k;
    dwdx[20] = Import_ASLpOH_k;
    dwdx[21] = k_CR_ASL_m;
    dwdx[22] = -Prot_k1;
    dwdx[23] = Prot_k1*(-fu_AS + 1)/fu_AS;
    dwdx[24] = ASL_c*CYP3A4_ASLoOH_Vmax*(-1/CYP3A4_ASpOH_Km1 - 1/CYP3A4_ASoOH_Km1)/(CYP3A4_ASLoOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
    dwdx[25] = ASL_c*CYP3A4_ASLpOH_Vmax*(-1/CYP3A4_ASpOH_Km1 - 1/CYP3A4_ASoOH_Km1)/(CYP3A4_ASLpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2));
    dwdx[26] = AS_c*CYP3A4_ASoOH_Vmax*(-1/CYP3A4_ASpOH_Km1 - 1/CYP3A4_ASoOH_Km1)/(CYP3A4_ASoOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) + CYP3A4_ASoOH_Vmax/(CYP3A4_ASoOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    dwdx[27] = AS_c*CYP3A4_ASpOH_Vmax*(-1/CYP3A4_ASpOH_Km1 - 1/CYP3A4_ASoOH_Km1)/(CYP3A4_ASpOH_Km1*pow(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1, 2)) + CYP3A4_ASpOH_Vmax/(CYP3A4_ASpOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    dwdx[28] = Export_AS_k;
    dwdx[29] = AS_c*UGT1A3_AS_Vmax*(-2*AS_c/UGT1A3_AS_KI1 - 1)/pow(pow(AS_c, 2)/UGT1A3_AS_KI1 + AS_c + UGT1A3_AS_Km1, 2) + UGT1A3_AS_Vmax/(pow(AS_c, 2)/UGT1A3_AS_KI1 + AS_c + UGT1A3_AS_Km1);
    dwdx[30] = Import_AS_k;
    dwdx[31] = -Prot_k1;
    dwdx[32] = Prot_k1*(-fu_AS + 1)/fu_AS;
    dwdx[33] = Export_ASoOH_k;
    dwdx[34] = Import_ASoOH_k;
    dwdx[35] = -Prot_k1;
    dwdx[36] = Prot_k1*(-fu_AS + 1)/fu_AS;
    dwdx[37] = Export_ASpOH_k;
    dwdx[38] = Import_ASpOH_k;
}