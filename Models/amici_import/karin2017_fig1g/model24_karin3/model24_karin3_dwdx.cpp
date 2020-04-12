#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model24_karin3(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = z*(-983040/(pow(amici_y, 6)*pow(1 + 32768/pow(amici_y, 5), 2)) + 403368.0/(pow(amici_y, 6)*pow(1 + 16807/pow(amici_y, 5), 2)));
    dwdx[1] = -mu*z;
    dwdx[2] = -0.10000000000000001 - 6/(1 + 32768/pow(amici_y, 5)) + 4.7999999999999998/(1 + 16807/pow(amici_y, 5));
    dwdx[3] = -amici_y*mu;
}