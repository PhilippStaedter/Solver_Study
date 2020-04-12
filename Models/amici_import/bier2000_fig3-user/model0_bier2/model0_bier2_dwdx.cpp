#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_bier2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = T1*k1;
    dwdx[1] = T2*k1;
    dwdx[2] = G1*k1;
    dwdx[3] = -T1*kp/pow(T1 + km, 2) + kp/(T1 + km);
    dwdx[4] = -epsilon;
    dwdx[5] = G2*k1;
    dwdx[6] = -T2*kp/pow(T2 + km, 2) + kp/(T2 + km);
    dwdx[7] = epsilon;
}