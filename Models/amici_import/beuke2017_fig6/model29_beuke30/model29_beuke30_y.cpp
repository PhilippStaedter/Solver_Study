#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model29_beuke30(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = IkBa_cytoplasm;
    y[1] = IkBa_mRNA_cytoplasm;
    y[2] = IkBa_mRNA_nucleus;
    y[3] = IkBa_pre_mRNA_1;
    y[4] = IkBa_pre_mRNA_2;
    y[5] = JNK;
    y[6] = JNK_P;
    y[7] = LPS;
    y[8] = MSK1;
    y[9] = MSK1_P;
    y[10] = TNFR1_EL;
    y[11] = TNFR1_cytoplasm;
    y[12] = TNFa_TNFR1_EL;
    y[13] = TNFa_TNFR1_cytoplasm;
    y[14] = TNFa_mRNA_LSEC;
    y[15] = TNFa_mRNA_MC;
    y[16] = TNFa_pre_mRNA_LSEC;
    y[17] = TNFa_pre_mRNA_MC;
    y[18] = p38;
    y[19] = p38_P;
    y[20] = p65_2P;
    y[21] = p65_IkBa_nucleus;
    y[22] = p65_cytoplasm;
    y[23] = p65_mRNA;
    y[24] = p65_nucleus;
    y[25] = s11;
    y[26] = s14;
    y[27] = s15;
    y[28] = s16;
    y[29] = s19;
    y[30] = s3;
    y[31] = s5;
    y[32] = s8;
}