#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_jiang1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ADP;
    x_solver[1] = ADP_cyt;
    x_solver[2] = AMP;
    x_solver[3] = ATP;
    x_solver[4] = ATP_cyt;
    x_solver[5] = Acetyl_CoA;
    x_solver[6] = Acetyl_CoA_cyt;
    x_solver[7] = Ala;
    x_solver[8] = Asp;
    x_solver[9] = Asp_cyt;
    x_solver[10] = CO2;
    x_solver[11] = Cit;
    x_solver[12] = Cit_cyt;
    x_solver[13] = CoA;
    x_solver[14] = CoA_cyt;
    x_solver[15] = Cytc2p;
    x_solver[16] = Cytc3p;
    x_solver[17] = DHAP;
    x_solver[18] = DPG;
    x_solver[19] = ETFox;
    x_solver[20] = ETFred;
    x_solver[21] = F6P;
    x_solver[22] = FAD;
    x_solver[23] = FADH2;
    x_solver[24] = FBP;
    x_solver[25] = Fum;
    x_solver[26] = G3P;
    x_solver[27] = GAP;
    x_solver[28] = GDP;
    x_solver[29] = GLC;
    x_solver[30] = GTP;
    x_solver[31] = Glu;
    x_solver[32] = Glu_cyt;
    x_solver[33] = H2O;
    x_solver[34] = IsoCit;
    x_solver[35] = IsoCitcyt;
    x_solver[36] = LAC;
    x_solver[37] = Mal;
    x_solver[38] = Mal_cyt;
    x_solver[39] = NAD;
    x_solver[40] = NADH;
    x_solver[41] = NADH_cyt;
    x_solver[42] = NADPH;
    x_solver[43] = NADPH_cyt;
    x_solver[44] = NADP_cyt;
    x_solver[45] = NADP_p;
    x_solver[46] = NAD_p;
    x_solver[47] = OG;
    x_solver[48] = OG_cyt;
    x_solver[49] = OXA;
    x_solver[50] = OXA_cyt;
    x_solver[51] = PEP;
    x_solver[52] = PYR_cyt;
    x_solver[53] = Pi;
    x_solver[54] = Pyr;
    x_solver[55] = Q;
    x_solver[56] = QH2;
    x_solver[57] = SCoA;
    x_solver[58] = Suc;
}