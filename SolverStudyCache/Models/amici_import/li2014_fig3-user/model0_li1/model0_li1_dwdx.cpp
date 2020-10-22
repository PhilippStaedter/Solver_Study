#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_li1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = alpha;
    dwdx[1] = r1;
    dwdx[2] = T*r1/Tmax + r1*(T + TStar)/Tmax;
    dwdx[3] = V*amici_k;
    dwdx[4] = V*amici_k*sigma;
    dwdx[5] = TStar*r2/Tmax;
    dwdx[6] = NN*beta;
    dwdx[7] = T*r1/Tmax;
    dwdx[8] = beta;
    dwdx[9] = r2;
    dwdx[10] = TStar*r2/Tmax + r2*(T + TStar)/Tmax;
    dwdx[11] = epsilon;
    dwdx[12] = T*amici_k;
    dwdx[13] = T*amici_k*sigma;
}