#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_karin5(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx1;
    J[1] = 1.0*dwdx5;
    J[4] = 1.0*dwdx11;
    J[6] = 1.0*dwdx3;
    J[7] = 1.0*dwdx7;
    J[10] = 1.0*dwdx0;
    J[11] = 1.0*dwdx4;
    J[12] = 1.0*dwdx8;
    J[13] = 1.0*dwdx9;
    J[15] = 1.0*dwdx2;
    J[16] = 1.0*dwdx6;
    J[18] = 1.0*dwdx10;
    J[19] = 1.0*dwdx12;
    J[24] = 1.0*dwdx13;
}