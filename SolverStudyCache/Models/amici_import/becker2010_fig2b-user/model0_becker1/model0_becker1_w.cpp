#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_becker1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Bmax*kt;
    w[1] = 1.0*EpoR*kt;
    w[2] = 1.0*Epo*EpoR*kon;
    w[3] = 1.0*Epo_EpoR*koff;
    w[4] = 1.0*Epo_EpoR*ke;
    w[5] = 1.0*Epo_EpoRi*kex;
    w[6] = 1.0*Epo_EpoRi*kdi;
    w[7] = 1.0*Epo_EpoRi*kde;
}