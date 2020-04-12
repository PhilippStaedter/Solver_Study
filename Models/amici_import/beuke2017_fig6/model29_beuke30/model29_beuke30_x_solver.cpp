#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model29_beuke30(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = IkBa_cytoplasm;
    x_solver[1] = IkBa_mRNA_cytoplasm;
    x_solver[2] = IkBa_mRNA_nucleus;
    x_solver[3] = IkBa_pre_mRNA_1;
    x_solver[4] = IkBa_pre_mRNA_2;
    x_solver[5] = JNK;
    x_solver[6] = JNK_P;
    x_solver[7] = LPS;
    x_solver[8] = MSK1;
    x_solver[9] = MSK1_P;
    x_solver[10] = TNFR1_EL;
    x_solver[11] = TNFR1_cytoplasm;
    x_solver[12] = TNFa_TNFR1_EL;
    x_solver[13] = TNFa_TNFR1_cytoplasm;
    x_solver[14] = TNFa_mRNA_LSEC;
    x_solver[15] = TNFa_mRNA_MC;
    x_solver[16] = TNFa_pre_mRNA_LSEC;
    x_solver[17] = TNFa_pre_mRNA_MC;
    x_solver[18] = p38;
    x_solver[19] = p38_P;
    x_solver[20] = p65_2P;
    x_solver[21] = p65_IkBa_nucleus;
    x_solver[22] = p65_cytoplasm;
    x_solver[23] = p65_mRNA;
    x_solver[24] = p65_nucleus;
    x_solver[25] = s11;
    x_solver[26] = s14;
    x_solver[27] = s15;
    x_solver[28] = s16;
    x_solver[29] = s19;
    x_solver[30] = s3;
    x_solver[31] = s5;
    x_solver[32] = s8;
}