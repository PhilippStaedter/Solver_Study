#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_hornberg1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = R*Vm1/(Km1 + R);
    w[1] = Rin*Vm2/(Km2 + Rin);
    w[2] = R*k3*x1/(Km3 + x1);
    w[3] = Vm4*x1p/(Km4 + x1p);
    w[4] = k5*x1p*x2/(Km5 + x2);
    w[5] = Vm6*x2p/(Km6 + x2p);
    w[6] = k7*x2p*x3/(Km7 + x3);
    w[7] = Vm8*x3p/(Km8*(Inh/Ki8 + 1 + x3p/Km8));
}