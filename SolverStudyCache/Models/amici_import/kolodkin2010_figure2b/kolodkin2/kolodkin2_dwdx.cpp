#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_kolodkin2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = Kapf/Vnucleus;
    dwdx[1] = -Kapb/Vnucleus;
    dwdx[2] = NRn*k12;
    dwdx[3] = RE*k1;
    dwdx[4] = -k22;
    dwdx[5] = Ln_*k12;
    dwdx[6] = NRLn*k1;
    dwdx[7] = -k2;
}