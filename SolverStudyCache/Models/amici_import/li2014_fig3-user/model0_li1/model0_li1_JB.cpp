#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_li1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 - 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3;
    JB[1] = -1.0*dwdx4 + 1.0*dwdx5;
    JB[3] = 1.0*dwdx7;
    JB[4] = 1.0*dwdx10 + 1.0*dwdx8 - 1.0*dwdx9;
    JB[5] = -1.0*dwdx6;
    JB[6] = 1.0*dwdx12;
    JB[7] = -1.0*dwdx13;
    JB[8] = 1.0*dwdx11;
}