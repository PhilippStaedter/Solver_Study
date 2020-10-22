#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_hornberg1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -Vm8*x3p/(Ki8*Km8*pow(Inh/Ki8 + 1 + x3p/Km8, 2));
    dwdx[1] = -R*Vm1/pow(Km1 + R, 2) + Vm1/(Km1 + R);
    dwdx[2] = k3*x1/(Km3 + x1);
    dwdx[3] = -Rin*Vm2/pow(Km2 + Rin, 2) + Vm2/(Km2 + Rin);
    dwdx[4] = -R*k3*x1/pow(Km3 + x1, 2) + R*k3/(Km3 + x1);
    dwdx[5] = -Vm4*x1p/pow(Km4 + x1p, 2) + Vm4/(Km4 + x1p);
    dwdx[6] = k5*x2/(Km5 + x2);
    dwdx[7] = -k5*x1p*x2/pow(Km5 + x2, 2) + k5*x1p/(Km5 + x2);
    dwdx[8] = -Vm6*x2p/pow(Km6 + x2p, 2) + Vm6/(Km6 + x2p);
    dwdx[9] = k7*x3/(Km7 + x3);
    dwdx[10] = -k7*x2p*x3/pow(Km7 + x3, 2) + k7*x2p/(Km7 + x3);
    dwdx[11] = Vm8/(Km8*(Inh/Ki8 + 1 + x3p/Km8)) - Vm8*x3p/(pow(Km8, 2)*pow(Inh/Ki8 + 1 + x3p/Km8, 2));
}