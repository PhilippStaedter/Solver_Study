#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_jiang1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ADP;
    y[1] = ADP_cyt;
    y[2] = AMP;
    y[3] = ATP;
    y[4] = ATP_cyt;
    y[5] = Acetyl_CoA;
    y[6] = Acetyl_CoA_cyt;
    y[7] = Ala;
    y[8] = Asp;
    y[9] = Asp_cyt;
    y[10] = CO2;
    y[11] = Cit;
    y[12] = Cit_cyt;
    y[13] = CoA;
    y[14] = CoA_cyt;
    y[15] = Cytc2p;
    y[16] = Cytc3p;
    y[17] = DHAP;
    y[18] = DPG;
    y[19] = ETFox;
    y[20] = ETFred;
    y[21] = F6P;
    y[22] = FAD;
    y[23] = FADH2;
    y[24] = FBP;
    y[25] = Fum;
    y[26] = G3P;
    y[27] = GAP;
    y[28] = GDP;
    y[29] = GLC;
    y[30] = GTP;
    y[31] = Glu;
    y[32] = Glu_cyt;
    y[33] = H2O;
    y[34] = IsoCit;
    y[35] = IsoCitcyt;
    y[36] = LAC;
    y[37] = Mal;
    y[38] = Mal_cyt;
    y[39] = NAD;
    y[40] = NADH;
    y[41] = NADH_cyt;
    y[42] = NADPH;
    y[43] = NADPH_cyt;
    y[44] = NADP_cyt;
    y[45] = NADP_p;
    y[46] = NAD_p;
    y[47] = OG;
    y[48] = OG_cyt;
    y[49] = OXA;
    y[50] = OXA_cyt;
    y[51] = PEP;
    y[52] = PYR_cyt;
    y[53] = Pi;
    y[54] = Pyr;
    y[55] = Q;
    y[56] = QH2;
    y[57] = SCoA;
    y[58] = Suc;
}