#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_bachar1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = b*(I1 + I2 + S);
    w[1] = Beta*S*c*(I1 + I2*a)/(I1 + I2 + S);
    w[2] = S*d;
    w[3] = I1*d;
    w[4] = I1*k1;
    w[5] = Alpha*I2;
    w[6] = I2*d;
    w[7] = I2*k2;
}