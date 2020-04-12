#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_bucher1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = Prot_k1*(-ASL_b + ASL_c*(-fu_ASL + 1)/fu_ASL);
    w[1] = Prot_k1*(-ASLoOH_b + ASLoOH_c*(-fu_ASL + 1)/fu_ASL);
    w[2] = Prot_k1*(-ASLpOH_b + ASLpOH_c*(-fu_ASL + 1)/fu_ASL);
    w[3] = Prot_k1*(-AS_b + AS_c*(-fu_AS + 1)/fu_AS);
    w[4] = Prot_k1*(-ASoOH_b + ASoOH_c*(-fu_AS + 1)/fu_AS);
    w[5] = Prot_k1*(-ASpOH_b + ASpOH_c*(-fu_AS + 1)/fu_AS);
    w[6] = ASLoOH_c*(k_CR_ASL_c + k_PON_OH_c);
    w[7] = ASLpOH_c*(k_CR_ASL_c + k_PON_OH_c);
    w[8] = ASL_c*CYP3A4_ASLoOH_Vmax/(CYP3A4_ASLoOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    w[9] = ASL_c*CYP3A4_ASLpOH_Vmax/(CYP3A4_ASLpOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    w[10] = AS_c*CYP3A4_ASoOH_Vmax/(CYP3A4_ASoOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    w[11] = AS_c*CYP3A4_ASpOH_Vmax/(CYP3A4_ASpOH_Km1*(ASL_c/CYP3A4_ASLpOH_Km1 + ASL_c/CYP3A4_ASLoOH_Km1 + AS_c/CYP3A4_ASpOH_Km1 + AS_c/CYP3A4_ASoOH_Km1 + 1));
    w[12] = AS_c*Export_AS_k;
    w[13] = ASL_c*Export_ASL_k;
    w[14] = ASLoOH_c*Export_ASLoOH_k;
    w[15] = ASLpOH_c*Export_ASLpOH_k;
    w[16] = ASoOH_c*Export_ASoOH_k;
    w[17] = ASpOH_c*Export_ASpOH_k;
    w[18] = AS_m*Import_AS_k;
    w[19] = ASL_m*Import_ASL_k;
    w[20] = ASLoOH_m*Import_ASLoOH_k;
    w[21] = ASLpOH_m*Import_ASLpOH_k;
    w[22] = ASoOH_m*Import_ASoOH_k;
    w[23] = ASpOH_m*Import_ASpOH_k;
    w[24] = ASL_c*(k_CR_ASL_c + k_PON_ASL_c);
    w[25] = ASL_m*k_CR_ASL_m;
    w[26] = ASLoOH_m*k_CR_ASL_m;
    w[27] = ASLpOH_m*k_CR_ASL_m;
    w[28] = AS_c*UGT1A3_AS_Vmax/(pow(AS_c, 2)/UGT1A3_AS_KI1 + AS_c + UGT1A3_AS_Km1);
}