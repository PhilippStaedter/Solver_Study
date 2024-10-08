#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_esteban2013(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = -1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2;
    JB[1] = -1.0*dwdx1;
    JB[6] = 1.0*dwdx3;
    JB[9] = -1.0*dwdx4;
    JB[12] = -1.0*dwdx5 + 1.0*dwdx6 + 1.0*dwdx7;
    JB[13] = -1.0*dwdx6;
    JB[18] = 1.0*dwdx9;
    JB[19] = -1.0*dwdx8;
    JB[20] = -1.0*dwdx13 + 1.0*dwdx14;
    JB[21] = -1.0*dwdx14;
    JB[22] = -1.0*dwdx11 + 1.0*dwdx12;
    JB[23] = -1.0*dwdx12;
    JB[24] = 1.0*dwdx10;
}