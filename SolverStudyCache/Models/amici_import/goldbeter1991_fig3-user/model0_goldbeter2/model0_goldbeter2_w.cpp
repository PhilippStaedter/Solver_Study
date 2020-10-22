#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_goldbeter2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = cell*vi;
    w[1] = Cvar*cell*kd;
    w[2] = Cvar*X*cell*vd/(Cvar + Kd);
    w[3] = Cvar*VM1*cell*(-Mvar + 1)/((Cvar + Kc)*(K1 - Mvar + 1));
    w[4] = Mvar*V2*cell/(K2 + Mvar);
    w[5] = Mvar*VM3*cell*(-X + 1)/(K3 - X + 1);
    w[6] = V4*X*cell/(K4 + X);
}