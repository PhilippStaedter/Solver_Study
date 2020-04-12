#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_goldbeter6(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*P*reaction_1_V/pow(1.0*P + reaction_1_Km, 2) + 1.0*reaction_1_V/(1.0*P + reaction_1_Km);
    dwdx[1] = 1.0*reaction_4_V3*(-1.0*R + 1)/(-1.0*R + reaction_4_k3 + 1);
    dwdx[2] = 1.0*reaction_0_a;
    dwdx[3] = 1.0*reaction_2_V1*(-1.0*Q + 1)/pow(-1.0*Q + reaction_2_K1 + 1, 2) - 1.0*reaction_2_V1/(-1.0*Q + reaction_2_K1 + 1);
    dwdx[4] = -1.0*Q*R*reaction_3_V2/pow(1.0*Q + reaction_3_K2, 2) + 1.0*R*reaction_3_V2/(1.0*Q + reaction_3_K2);
    dwdx[5] = 1.0*Q*reaction_3_V2/(1.0*Q + reaction_3_K2);
    dwdx[6] = 1.0*P*reaction_4_V3*(-1.0*R + 1)/pow(-1.0*R + reaction_4_k3 + 1, 2) - 1.0*P*reaction_4_V3/(-1.0*R + reaction_4_k3 + 1);
    dwdx[7] = -1.0*R*reaction_5_V/pow(1.0*R + reaction_5_Km, 2) + 1.0*reaction_5_V/(1.0*R + reaction_5_Km);
}