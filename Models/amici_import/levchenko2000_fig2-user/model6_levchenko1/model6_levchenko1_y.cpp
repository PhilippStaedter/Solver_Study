#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model6_levchenko1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = MAPK;
    y[1] = MAPKMEKpp;
    y[2] = MAPKPH;
    y[3] = MAPKp;
    y[4] = MAPKpMAPKPH;
    y[5] = MAPKpMEKpp;
    y[6] = MAPKpp;
    y[7] = MAPKppMAPKPH;
    y[8] = MEK;
    y[9] = MEKPH;
    y[10] = MEKRAFp;
    y[11] = MEKp;
    y[12] = MEKpMEKPH;
    y[13] = MEKpRAFp;
    y[14] = MEKpp;
    y[15] = MEKppMEKPH;
    y[16] = RAF;
    y[17] = RAFK;
    y[18] = RAFPH;
    y[19] = RAFRAFK;
    y[20] = RAFp;
    y[21] = RAFpRAFPH;
}