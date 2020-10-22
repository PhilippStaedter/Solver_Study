#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model5_panteleev3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*VIIa_TF*X*reaction_1_k1 - 1.0*VIIa_TF_X*reaction_1_k2;
    w[1] = 1.0*VIIa_TF_X*reaction_2_k1;
    w[2] = -1.0*VIIa_TF*Xa*reaction_3_k2 + 1.0*VIIa_TF_Xa*reaction_3_k1;
    w[3] = 1.0*TFPI*Xa*reaction_4_k1 - 1.0*Xa_TFPI*reaction_4_k2;
    w[4] = 1.0*VIIa_TF*Xa_TFPI*reaction_5_k1 - 1.0*Xa_TFPI_VIIa_TF*reaction_5_k2;
}