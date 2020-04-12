#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_lee2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*E*P*k1 - 1.0*E_P_1*j1;
    w[1] = 1.0*E*P2*k7a - 1.0*E_P2*j7a;
    w[2] = 1.0*E_P2*k8a;
    w[3] = 1.0*E_P_1*k2;
    w[4] = 1.0*E*M*k3a - 1.0*E_M*j3a;
    w[5] = 1.0*E_M*k4a;
    w[6] = 1.0*E*P*k5 - 1.0*E_P_2*j5;
    w[7] = 1.0*E_P_2*k6;
}