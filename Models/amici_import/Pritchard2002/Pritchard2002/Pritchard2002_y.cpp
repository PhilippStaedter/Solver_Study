#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Pritchard2002(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = GLCo;
    y[1] = GLCi;
    y[2] = ATP;
    y[3] = G6P;
    y[4] = ADP;
    y[5] = F6P;
    y[6] = F16bP;
    y[7] = AMP;
    y[8] = F26bP;
    y[9] = DHAP;
    y[10] = GAP;
    y[11] = NAD;
    y[12] = BPG;
    y[13] = NADH;
    y[14] = P3G;
    y[15] = P2G;
    y[16] = PEP;
    y[17] = PYR;
    y[18] = AcAld;
    y[19] = CO2;
    y[20] = EtOH;
    y[21] = Glycerol;
    y[22] = Glycogen;
    y[23] = Trehalose;
    y[24] = Succinate;
}