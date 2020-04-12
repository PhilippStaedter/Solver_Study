#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_hald(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ACA;
    y[1] = ACA0;
    y[2] = ADP;
    y[3] = AMP;
    y[4] = ATP;
    y[5] = DHAP;
    y[6] = DHAPCN;
    y[7] = DPG;
    y[8] = EtOH;
    y[9] = EtOH0;
    y[10] = F6P;
    y[11] = FBP;
    y[12] = G6P;
    y[13] = GAP;
    y[14] = Glc;
    y[15] = Glc0;
    y[16] = Glyc;
    y[17] = Glyc0;
    y[18] = HCN;
    y[19] = HCN0;
    y[20] = NAD;
    y[21] = NADH;
    y[22] = OAc;
    y[23] = OAc0;
    y[24] = PEP;
    y[25] = Pyr;
    y[26] = PyrCN;
    y[27] = X;
    y[28] = drain;
    y[29] = glycogen;
    y[30] = lacto;
}