#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_goldbeter6(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Q*reaction_0_a;
    w[1] = 1.0*P*reaction_1_V/(1.0*P + reaction_1_Km);
    w[2] = reaction_2_V1*(-1.0*Q + 1)/(-1.0*Q + reaction_2_K1 + 1);
    w[3] = 1.0*Q*R*reaction_3_V2/(1.0*Q + reaction_3_K2);
    w[4] = 1.0*P*reaction_4_V3*(-1.0*R + 1)/(-1.0*R + reaction_4_k3 + 1);
    w[5] = 1.0*R*reaction_5_V/(1.0*R + reaction_5_Km);
}