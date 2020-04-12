#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_goldbeter2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = cell*kd;
    dwdx[1] = -Cvar*X*cell*vd/pow(Cvar + Kd, 2) + X*cell*vd/(Cvar + Kd);
    dwdx[2] = -Cvar*VM1*cell*(-Mvar + 1)/(pow(Cvar + Kc, 2)*(K1 - Mvar + 1)) + VM1*cell*(-Mvar + 1)/((Cvar + Kc)*(K1 - Mvar + 1));
    dwdx[3] = Cvar*VM1*cell*(-Mvar + 1)/((Cvar + Kc)*pow(K1 - Mvar + 1, 2)) - Cvar*VM1*cell/((Cvar + Kc)*(K1 - Mvar + 1));
    dwdx[4] = -Mvar*V2*cell/pow(K2 + Mvar, 2) + V2*cell/(K2 + Mvar);
    dwdx[5] = VM3*cell*(-X + 1)/(K3 - X + 1);
    dwdx[6] = Cvar*cell*vd/(Cvar + Kd);
    dwdx[7] = Mvar*VM3*cell*(-X + 1)/pow(K3 - X + 1, 2) - Mvar*VM3*cell/(K3 - X + 1);
    dwdx[8] = -V4*X*cell/pow(K4 + X, 2) + V4*cell/(K4 + X);
}