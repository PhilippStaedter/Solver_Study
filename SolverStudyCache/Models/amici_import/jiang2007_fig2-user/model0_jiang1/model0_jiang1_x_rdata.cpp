#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_jiang1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ADP;
    x_rdata[1] = ADP_cyt;
    x_rdata[2] = AMP;
    x_rdata[3] = ATP;
    x_rdata[4] = ATP_cyt;
    x_rdata[5] = Acetyl_CoA;
    x_rdata[6] = Acetyl_CoA_cyt;
    x_rdata[7] = Ala;
    x_rdata[8] = Asp;
    x_rdata[9] = Asp_cyt;
    x_rdata[10] = CO2;
    x_rdata[11] = Cit;
    x_rdata[12] = Cit_cyt;
    x_rdata[13] = CoA;
    x_rdata[14] = CoA_cyt;
    x_rdata[15] = Cytc2p;
    x_rdata[16] = Cytc3p;
    x_rdata[17] = DHAP;
    x_rdata[18] = DPG;
    x_rdata[19] = ETFox;
    x_rdata[20] = ETFred;
    x_rdata[21] = F6P;
    x_rdata[22] = FAD;
    x_rdata[23] = FADH2;
    x_rdata[24] = FBP;
    x_rdata[25] = Fum;
    x_rdata[26] = G3P;
    x_rdata[27] = GAP;
    x_rdata[28] = GDP;
    x_rdata[29] = GLC;
    x_rdata[30] = GTP;
    x_rdata[31] = Glu;
    x_rdata[32] = Glu_cyt;
    x_rdata[33] = H2O;
    x_rdata[34] = IsoCit;
    x_rdata[35] = IsoCitcyt;
    x_rdata[36] = LAC;
    x_rdata[37] = Mal;
    x_rdata[38] = Mal_cyt;
    x_rdata[39] = NAD;
    x_rdata[40] = NADH;
    x_rdata[41] = NADH_cyt;
    x_rdata[42] = NADPH;
    x_rdata[43] = NADPH_cyt;
    x_rdata[44] = NADP_cyt;
    x_rdata[45] = NADP_p;
    x_rdata[46] = NAD_p;
    x_rdata[47] = OG;
    x_rdata[48] = OG_cyt;
    x_rdata[49] = OXA;
    x_rdata[50] = OXA_cyt;
    x_rdata[51] = PEP;
    x_rdata[52] = PYR_cyt;
    x_rdata[53] = Pi;
    x_rdata[54] = Pyr;
    x_rdata[55] = Q;
    x_rdata[56] = QH2;
    x_rdata[57] = SCoA;
    x_rdata[58] = Suc;
}