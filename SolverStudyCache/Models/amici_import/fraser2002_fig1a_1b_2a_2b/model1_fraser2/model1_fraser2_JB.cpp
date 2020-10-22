#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model1_fraser2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = -1.0*dwdx0 + 1.0*dwdx1;
    JB[3] = -1.0*dwdx2;
    JB[7] = -1.0*dwdx3 + 1.0*dwdx4;
    JB[10] = -1.0*dwdx5;
    JB[14] = 1.0*dwdx6 - 1.0*dwdx8;
    JB[17] = -1.0*dwdx7;
    JB[18] = 1.0*dwdx10;
    JB[21] = -1.0*dwdx11 - 1.0*dwdx12 + 1.0*dwdx9;
    JB[25] = 1.0*dwdx14;
    JB[28] = 1.0*dwdx13 - 1.0*dwdx15 - 1.0*dwdx16;
    JB[32] = 1.0*dwdx17;
    JB[35] = -1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20;
}