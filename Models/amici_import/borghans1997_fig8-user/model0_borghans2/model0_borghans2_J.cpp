#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_borghans2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx1 - 1.0*dwdx2;
    J[3] = -1.0*dwdx8;
    J[7] = 1.0*dwdx7;
    J[8] = -1.0*dwdx0;
    J[10] = -1.0*dwdx3 - 1.0*dwdx4;
    J[11] = 1.0*dwdx5 - 1.0*dwdx6;
    J[12] = 1.0*dwdx0;
    J[14] = 1.0*dwdx3 + 1.0*dwdx4;
    J[15] = -1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7;
}