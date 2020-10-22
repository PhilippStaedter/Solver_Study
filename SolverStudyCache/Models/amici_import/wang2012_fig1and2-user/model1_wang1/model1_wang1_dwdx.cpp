#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model1_wang1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = alpha;
    dwdx[1] = r;
    dwdx[2] = T*r/Tmax + r*(T + TStar)/Tmax;
    dwdx[3] = V*amici_k;
    dwdx[4] = V*amici_k*sigma1;
    dwdx[5] = T*r/Tmax;
    dwdx[6] = beta;
    dwdx[7] = NN*beta;
    dwdx[8] = T*amici_k;
    dwdx[9] = T*amici_k*sigma1;
    dwdx[10] = gamma;
}