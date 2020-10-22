#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_delcontezerail1_Fig4D(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx2 - 1.0*dwdx3;
    J[1] = -1.0*dwdx5;
    J[2] = 1.0*dwdx7;
    J[4] = 1.0*dwdx1;
    J[5] = 1.0*dwdx4 - 1.0*dwdx6;
    J[7] = 1.0*dwdx10 + 1.0*dwdx9;
    J[8] = -1.0*dwdx0 + 1.0*dwdx2 + 1.0*dwdx3;
    J[9] = 1.0*dwdx5;
    J[10] = -1.0*dwdx7 - 1.0*dwdx8;
    J[12] = -1.0*dwdx1;
    J[13] = -1.0*dwdx4 + 1.0*dwdx6;
    J[15] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx9;
}