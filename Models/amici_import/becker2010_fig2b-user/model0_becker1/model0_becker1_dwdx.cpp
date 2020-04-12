#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_becker1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*EpoR*kon;
    dwdx[1] = 1.0*kt;
    dwdx[2] = 1.0*Epo*kon;
    dwdx[3] = 1.0*koff;
    dwdx[4] = 1.0*ke;
    dwdx[5] = 1.0*kex;
    dwdx[6] = 1.0*kdi;
    dwdx[7] = 1.0*kde;
}