#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_hatakeyama1_Fig5G(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Akt;
    y[1] = AktPIP;
    y[2] = AktPIP3;
    y[3] = AktPIPP;
    y[4] = E;
    y[5] = ERK;
    y[6] = ERKP;
    y[7] = ERKPP;
    y[8] = GS;
    y[9] = HRG;
    y[10] = MEK;
    y[11] = MEKP;
    y[12] = MEKPP;
    y[13] = MKP3;
    y[14] = PI3K;
    y[15] = PI3Kstar;
    y[16] = PIP3;
    y[17] = PP2A;
    y[18] = P_I;
    y[19] = R;
    y[20] = RHRG;
    y[21] = RHRG2;
    y[22] = RP;
    y[23] = RPI3K;
    y[24] = RPI3Kstar;
    y[25] = RShGS;
    y[26] = RShP;
    y[27] = RShc;
    y[28] = Raf;
    y[29] = Rafstar;
    y[30] = RasGDP;
    y[31] = RasGTP;
    y[32] = ShGS;
    y[33] = ShP;
    y[34] = Shc;
    y[35] = internalization;
}