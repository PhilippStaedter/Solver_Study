#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_kirschner(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = s1;
    w[1] = V*s2/(V + b1);
    w[2] = T*mu;
    w[3] = T*V*amici_k;
    w[4] = V*g/(V + b2);
    w[5] = T*V*c;
}