#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model24_karin3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = z*(-0.10000000000000001 - 6/(1 + 32768/pow(amici_y, 5)) + 4.7999999999999998/(1 + 16807/pow(amici_y, 5)));
    w[1] = mu*(M - amici_y*z);
}