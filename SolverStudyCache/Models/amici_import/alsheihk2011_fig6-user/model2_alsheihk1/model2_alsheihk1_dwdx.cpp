#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_alsheihk1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = d;
    dwdx[1] = Mu;
    dwdx[2] = Beta_1*S;
    dwdx[3] = Theta;
    dwdx[4] = Mu;
    dwdx[5] = Delta;
    dwdx[6] = Beta_2*S;
    dwdx[7] = Mu;
    dwdx[8] = Delta;
    dwdx[9] = Beta_1*I1 + Beta_2*I2;
    dwdx[10] = Mu;
}