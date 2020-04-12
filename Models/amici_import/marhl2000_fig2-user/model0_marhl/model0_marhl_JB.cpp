#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_marhl(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 0.25*dwdx0 + 0.25*dwdx1;
    JB[3] = -1.0*dwdx0 - 1.0*dwdx1;
    JB[6] = 0.25*dwdx2;
    JB[8] = -1.0*dwdx2;
    JB[12] = 1.0*dwdx3;
    JB[13] = -1.0*dwdx3;
    JB[14] = -1.0*dwdx3;
    JB[15] = 0.25*dwdx4 + 0.25*dwdx6 - 0.25*dwdx7;
    JB[16] = 0.25*dwdx8 - 0.25*dwdx9;
    JB[17] = -1.0*dwdx5;
    JB[18] = -1.0*dwdx4 + 1.0*dwdx5 - 1.0*dwdx6 + 1.0*dwdx7 - 1.0*dwdx8 + 1.0*dwdx9;
    JB[19] = 1.0*dwdx5;
    JB[22] = -1.0*dwdx10;
    JB[23] = 1.0*dwdx10;
    JB[24] = 1.0*dwdx10;
}