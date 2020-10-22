#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_gorlich1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = GDP;
    y[1] = GTP;
    y[2] = RCC1;
    y[3] = RCC1_Ran;
    y[4] = RCC1_RanGDP;
    y[5] = RCC1_RanGTP;
    y[6] = RanBP1;
    y[7] = RanGAP;
    y[8] = RanGDP_cy;
    y[9] = RanGDP_nuc;
    y[10] = RanGTP_RanBP1;
    y[11] = RanGTP_cy;
    y[12] = RanGTP_nuc;
}