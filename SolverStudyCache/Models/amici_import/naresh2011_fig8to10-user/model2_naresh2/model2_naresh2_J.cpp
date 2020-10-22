#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model2_naresh2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1;
    J[2] = 1.0*dwdx9;
    J[3] = 1.0*dwdx12;
    J[10] = 1.0*dwdx2 + 1.0*dwdx4;
    J[12] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx5 + 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    J[13] = 1.0*dwdx14 + 1.0*dwdx16 - 1.0*dwdx17;
    J[14] = 1.0*dwdx18 + 1.0*dwdx21;
    J[15] = 1.0*dwdx3;
    J[17] = 1.0*dwdx11 + 1.0*dwdx6 + 1.0*dwdx8;
    J[18] = -1.0*dwdx12 - 1.0*dwdx13 + 1.0*dwdx15 + 1.0*dwdx17;
    J[19] = 1.0*dwdx19;
    J[20] = -1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4;
    J[22] = -1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7;
    J[23] = -1.0*dwdx14 - 1.0*dwdx15 - 1.0*dwdx16;
    J[24] = -1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21;
}