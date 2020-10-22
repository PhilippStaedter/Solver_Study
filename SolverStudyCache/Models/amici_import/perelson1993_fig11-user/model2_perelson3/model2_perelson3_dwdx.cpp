#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_perelson3(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = V*k1;
    dwdx[1] = muT;
    dwdx[2] = r*(1 - Ttot/Tmax);
    dwdx[3] = k2;
    dwdx[4] = muT;
    dwdx[5] = N0*muB;
    dwdx[6] = muB;
    dwdx[7] = -T*r/Tmax;
    dwdx[8] = T*k1;
    dwdx[9] = muVV;
    dwdx[10] = -s*theta/pow(V + theta, 2);
}