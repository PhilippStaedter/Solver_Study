#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_model0_gorlich1(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = 55555555555.555557;
    dxdotdw[1] = -83333333333.333328;
    dxdotdw[2] = 83333333333.333328;
    dxdotdw[3] = -83333333333.333328;
    dxdotdw[4] = -83333333333.333328;
    dxdotdw[5] = 83333333333.333328;
    dxdotdw[6] = 55555555555.555557;
    dxdotdw[7] = -83333333333.333328;
    dxdotdw[8] = -83333333333.333328;
    dxdotdw[9] = 83333333333.333328;
    dxdotdw[10] = -83333333333.333328;
    dxdotdw[11] = 55555555555.555557;
    dxdotdw[12] = 55555555555.555557;
    dxdotdw[13] = -55555555555.555557;
    dxdotdw[14] = 55555555555.555557;
    dxdotdw[15] = -55555555555.555557;
    dxdotdw[16] = -55555555555.555557;
    dxdotdw[17] = 55555555555.555557;
    dxdotdw[18] = -55555555555.555557;
    dxdotdw[19] = 83333333333.333328;
    dxdotdw[20] = -83333333333.333328;
    dxdotdw[21] = 83333333333.333328;
}