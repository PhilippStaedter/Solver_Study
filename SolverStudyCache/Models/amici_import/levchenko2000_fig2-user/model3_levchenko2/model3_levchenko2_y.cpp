#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model3_levchenko2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = C1;
    y[1] = C2;
    y[2] = C3;
    y[3] = C4;
    y[4] = C5;
    y[5] = C6;
    y[6] = C7;
    y[7] = C8;
    y[8] = C9;
    y[9] = MAPK;
    y[10] = MAPKMEKpp;
    y[11] = MAPKPH;
    y[12] = MAPKp;
    y[13] = MAPKpMAPKPH;
    y[14] = MAPKpMEKpp;
    y[15] = MAPKpp;
    y[16] = MAPKppMAPKPH;
    y[17] = MEK;
    y[18] = MEKPH;
    y[19] = MEKRAFp;
    y[20] = MEKp;
    y[21] = MEKpMEKPH;
    y[22] = MEKpRAFp;
    y[23] = MEKpp;
    y[24] = MEKppMEKPH;
    y[25] = RAF;
    y[26] = RAFK;
    y[27] = RAFPH;
    y[28] = RAFRAFK;
    y[29] = RAFp;
    y[30] = RAFpRAFPH;
}