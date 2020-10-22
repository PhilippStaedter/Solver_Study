#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_tripathi1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = d;
    dwdx[1] = -S*(Beta_1*I1 + Beta_2*I2)/pow(A + I1 + I2 + S, 2);
    dwdx[2] = Alpha;
    dwdx[3] = Beta_1*S/(A + I1 + I2 + S) - S*(Beta_1*I1 + Beta_2*I2)/pow(A + I1 + I2 + S, 2);
    dwdx[4] = Theta;
    dwdx[5] = Delta;
    dwdx[6] = d;
    dwdx[7] = Beta_2*S/(A + I1 + I2 + S) - S*(Beta_1*I1 + Beta_2*I2)/pow(A + I1 + I2 + S, 2);
    dwdx[8] = Delta;
    dwdx[9] = d;
    dwdx[10] = -S*(Beta_1*I1 + Beta_2*I2)/pow(A + I1 + I2 + S, 2) + (Beta_1*I1 + Beta_2*I2)/(A + I1 + I2 + S);
    dwdx[11] = d;
}