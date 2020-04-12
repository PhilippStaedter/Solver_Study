#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model29_beuke30(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = IkBa_cytoplasm;
    x_rdata[1] = IkBa_mRNA_cytoplasm;
    x_rdata[2] = IkBa_mRNA_nucleus;
    x_rdata[3] = IkBa_pre_mRNA_1;
    x_rdata[4] = IkBa_pre_mRNA_2;
    x_rdata[5] = JNK;
    x_rdata[6] = JNK_P;
    x_rdata[7] = LPS;
    x_rdata[8] = MSK1;
    x_rdata[9] = MSK1_P;
    x_rdata[10] = TNFR1_EL;
    x_rdata[11] = TNFR1_cytoplasm;
    x_rdata[12] = TNFa_TNFR1_EL;
    x_rdata[13] = TNFa_TNFR1_cytoplasm;
    x_rdata[14] = TNFa_mRNA_LSEC;
    x_rdata[15] = TNFa_mRNA_MC;
    x_rdata[16] = TNFa_pre_mRNA_LSEC;
    x_rdata[17] = TNFa_pre_mRNA_MC;
    x_rdata[18] = p38;
    x_rdata[19] = p38_P;
    x_rdata[20] = p65_2P;
    x_rdata[21] = p65_IkBa_nucleus;
    x_rdata[22] = p65_cytoplasm;
    x_rdata[23] = p65_mRNA;
    x_rdata[24] = p65_nucleus;
    x_rdata[25] = s11;
    x_rdata[26] = s14;
    x_rdata[27] = s15;
    x_rdata[28] = s16;
    x_rdata[29] = s19;
    x_rdata[30] = s3;
    x_rdata[31] = s5;
    x_rdata[32] = s8;
}