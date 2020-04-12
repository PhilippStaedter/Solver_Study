#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_gardner1_Fig4B(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = reaction1_vi;
    w[1] = Z*reaction10_alpha*reaction10_d1;
    w[2] = Z*reaction11_alpha*reaction11_kd;
    w[3] = reaction12_vs;
    w[4] = Y*reaction13_d1;
    w[5] = C*X*reaction2_k1/(C + reaction2_K5);
    w[6] = C*reaction3_kd;
    w[7] = C*V1p*(-M + 1)/((C + K6)*(-M + reaction4_K1 + 1));
    w[8] = M*reaction5_V2/(M + reaction5_K2);
    w[9] = M*V3p*(-X + 1)/(-X + reaction6_K3 + 1);
    w[10] = X*reaction7_V4/(X + reaction7_K4);
    w[11] = C*Y*reaction8_a1;
    w[12] = Z*reaction9_a2;
}