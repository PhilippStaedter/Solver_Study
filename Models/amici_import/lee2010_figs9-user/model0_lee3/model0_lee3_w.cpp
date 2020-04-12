#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_lee3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*II*r1_k1;
    w[1] = 1.0*M*r2_k1;
    w[2] = 1.0*II*r3_k1;
    w[3] = 1.0*P2*r4_k1;
}