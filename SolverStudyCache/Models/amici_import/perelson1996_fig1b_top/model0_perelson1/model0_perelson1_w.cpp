#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_perelson1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = T0*Vin*amici_k;
    w[1] = Tstar*delta;
    w[2] = Vin*c;
    w[3] = Vni*c;
    w[4] = NN*Tstar*delta;
}