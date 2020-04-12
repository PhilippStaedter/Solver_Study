#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_lebeda2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[1] = 1.0*dwdx1;
    J[5] = -1.0*dwdx1;
    J[11] = 1.0*dwdx2;
    J[12] = 1.0*dwdx0;
    J[15] = -1.0*dwdx2;
}