#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_panteleev3(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*Xa*reaction_4_k1;
    dwdx[1] = 1.0*X*reaction_1_k1;
    dwdx[2] = -1.0*Xa*reaction_3_k2;
    dwdx[3] = 1.0*Xa_TFPI*reaction_5_k1;
    dwdx[4] = -1.0*reaction_1_k2;
    dwdx[5] = 1.0*reaction_2_k1;
    dwdx[6] = 1.0*reaction_3_k1;
    dwdx[7] = 1.0*VIIa_TF*reaction_1_k1;
    dwdx[8] = -1.0*VIIa_TF*reaction_3_k2;
    dwdx[9] = 1.0*TFPI*reaction_4_k1;
    dwdx[10] = -1.0*reaction_4_k2;
    dwdx[11] = 1.0*VIIa_TF*reaction_5_k1;
    dwdx[12] = -1.0*reaction_5_k2;
}