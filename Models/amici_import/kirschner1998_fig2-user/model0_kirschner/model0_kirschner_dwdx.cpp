#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kirschner(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = mu;
    dwdx[1] = V*amici_k;
    dwdx[2] = V*c;
    dwdx[3] = -V*s2/pow(V + b1, 2) + s2/(V + b1);
    dwdx[4] = T*amici_k;
    dwdx[5] = -V*g/pow(V + b2, 2) + g/(V + b2);
    dwdx[6] = T*c;
}