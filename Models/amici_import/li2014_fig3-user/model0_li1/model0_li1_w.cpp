#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_li1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = s;
    w[1] = NN*TStar*beta;
    w[2] = V*epsilon;
    w[3] = T*alpha;
    w[4] = T*r1;
    w[5] = T*r1*(T + TStar)/Tmax;
    w[6] = T*V*amici_k;
    w[7] = T*V*amici_k*sigma;
    w[8] = TStar*beta;
    w[9] = TStar*r2;
    w[10] = TStar*r2*(T + TStar)/Tmax;
}