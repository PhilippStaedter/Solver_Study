#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_bier2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = G1*T1*k1;
    w[1] = G2*T2*k1;
    w[2] = V_in;
    w[3] = T1*kp/(T1 + km);
    w[4] = T2*kp/(T2 + km);
    w[5] = epsilon*(-T1 + T2);
}