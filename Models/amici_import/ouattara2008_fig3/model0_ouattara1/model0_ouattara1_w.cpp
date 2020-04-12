#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_ouattara1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = s;
    w[1] = T*delta;
    w[2] = T*V*beta;
    w[3] = TStar*mu;
    w[4] = TStar*amici_k;
    w[5] = V*c;
}