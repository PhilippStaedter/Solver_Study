#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_ma(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx1;
    J[1] = 1.0*dwdx3;
    J[3] = -1.0*dwdx8;
    J[8] = -1.0*dwdx4;
    J[12] = 1.0*dwdx14;
    J[15] = 1.0*dwdx5;
    J[16] = -1.0*dwdx6;
    J[17] = -1.0*dwdx10;
    J[24] = -1.0*dwdx9;
    J[27] = 1.0*dwdx16;
    J[30] = -1.0*dwdx7;
    J[32] = -1.0*dwdx12;
    J[35] = 1.0*dwdx0;
    J[40] = -1.0*dwdx13;
    J[42] = 1.0*dwdx2;
    J[46] = -1.0*dwdx11;
    J[48] = -1.0*dwdx15;
}