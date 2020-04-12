#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_model1_piedrafita1(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = -1.0;
    dxdotdw[1] = 1.0;
    dxdotdw[2] = 1.0;
    dxdotdw[3] = -1.0;
    dxdotdw[4] = 1.0;
    dxdotdw[5] = -1.0;
    dxdotdw[6] = -1.0;
    dxdotdw[7] = 1.0;
    dxdotdw[8] = 1.0;
    dxdotdw[9] = 1.0;
    dxdotdw[10] = -1.0;
    dxdotdw[11] = -1.0;
    dxdotdw[12] = -1.0;
    dxdotdw[13] = -1.0;
    dxdotdw[14] = 1.0;
    dxdotdw[15] = -1.0;
    dxdotdw[16] = 1.0;
    dxdotdw[17] = 1.0;
    dxdotdw[18] = 1.0;
    dxdotdw[19] = -1.0;
    dxdotdw[20] = -1.0;
    dxdotdw[21] = -1.0;
    dxdotdw[22] = 1.0;
}