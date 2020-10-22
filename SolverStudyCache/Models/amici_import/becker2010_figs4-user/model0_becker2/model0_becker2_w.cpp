#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_becker2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Bmax_SAv*kt;
    w[1] = 1.0*EpoR*kt;
    w[2] = 1.0*EpoR*SAv*kon_SAv;
    w[3] = 1.0*SAv_EpoR*koff_SAv;
    w[4] = 1.0*SAv_EpoR*kt;
    w[5] = 1.0*SAv_EpoRi*kex_SAv;
    w[6] = 1.0*SAv_EpoRi*kdi;
    w[7] = 1.0*SAv_EpoRi*kde;
}