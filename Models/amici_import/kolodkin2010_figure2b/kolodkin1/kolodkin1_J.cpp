#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_kolodkin1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 + 1.0*dwdx1;
    J[1] = -1.0*dwdx2;
    J[2] = -1.0*dwdx3;
    J[3] = 1.0*dwdx0;
    J[4] = 1.0*dwdx2;
    J[5] = 1.0*dwdx3;
    J[6] = -1.0*dwdx0;
    J[7] = -1.0*dwdx2;
    J[8] = -1.0*dwdx3;
}