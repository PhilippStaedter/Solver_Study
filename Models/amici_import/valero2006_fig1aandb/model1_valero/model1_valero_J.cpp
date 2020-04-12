#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model1_valero(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[1] = 2.0*dwdx1;
    J[2] = 2.0*dwdx3;
    J[7] = -1.0*dwdx1;
    J[8] = 1.0*dwdx2 - 1.0*dwdx3;
    J[12] = 1.0*dwdx0;
    J[13] = -1.0*dwdx1;
    J[14] = -1.0*dwdx2 - 1.0*dwdx3;
    J[29] = -1.0*dwdx4;
    J[30] = 1.0*dwdx0;
    J[35] = -1.0*dwdx4;
}