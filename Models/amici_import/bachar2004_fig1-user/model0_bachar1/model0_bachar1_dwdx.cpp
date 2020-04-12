#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_bachar1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = b;
    dwdx[1] = -Beta*S*c*(I1 + I2*a)/pow(I1 + I2 + S, 2) + Beta*S*c/(I1 + I2 + S);
    dwdx[2] = d;
    dwdx[3] = k1;
    dwdx[4] = b;
    dwdx[5] = Beta*S*a*c/(I1 + I2 + S) - Beta*S*c*(I1 + I2*a)/pow(I1 + I2 + S, 2);
    dwdx[6] = Alpha;
    dwdx[7] = d;
    dwdx[8] = k2;
    dwdx[9] = b;
    dwdx[10] = -Beta*S*c*(I1 + I2*a)/pow(I1 + I2 + S, 2) + Beta*c*(I1 + I2*a)/(I1 + I2 + S);
    dwdx[11] = d;
}