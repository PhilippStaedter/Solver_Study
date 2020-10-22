#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model1_fraser2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx1;
    J[3] = -1.0*dwdx10;
    J[7] = 1.0*dwdx3 - 1.0*dwdx4;
    J[10] = -1.0*dwdx14;
    J[14] = -1.0*dwdx6 + 1.0*dwdx8;
    J[17] = -1.0*dwdx17;
    J[18] = 1.0*dwdx2;
    J[21] = 1.0*dwdx11 + 1.0*dwdx12 - 1.0*dwdx9;
    J[25] = 1.0*dwdx5;
    J[28] = -1.0*dwdx13 + 1.0*dwdx15 + 1.0*dwdx16;
    J[32] = 1.0*dwdx7;
    J[35] = 1.0*dwdx18 + 1.0*dwdx19 - 1.0*dwdx20;
}