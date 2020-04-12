#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_wang1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = s;
    w[1] = T*alpha;
    w[2] = T*r;
    w[3] = T*r*(T + TStar)/Tmax;
    w[4] = T*V*amici_k;
    w[5] = T*V*amici_k*sigma1;
    w[6] = TStar*beta;
    w[7] = NN*TStar*beta;
    w[8] = V*gamma;
}