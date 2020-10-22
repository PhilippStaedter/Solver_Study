#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_hornberg1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[10] = -1.0*dwdx1;
    J[11] = 1.0*dwdx3;
    J[19] = 1.0*dwdx1;
    J[20] = -1.0*dwdx3;
    J[28] = -1.0*dwdx2;
    J[30] = -1.0*dwdx4;
    J[31] = 1.0*dwdx5;
    J[37] = 1.0*dwdx2;
    J[39] = 1.0*dwdx4;
    J[40] = -1.0*dwdx5;
    J[49] = -1.0*dwdx6;
    J[50] = -1.0*dwdx7;
    J[51] = 1.0*dwdx8;
    J[58] = 1.0*dwdx6;
    J[59] = 1.0*dwdx7;
    J[60] = -1.0*dwdx8;
    J[63] = 1.0*dwdx0;
    J[69] = -1.0*dwdx9;
    J[70] = -1.0*dwdx10;
    J[71] = 1.0*dwdx11;
    J[72] = -1.0*dwdx0;
    J[78] = 1.0*dwdx9;
    J[79] = 1.0*dwdx10;
    J[80] = -1.0*dwdx11;
}