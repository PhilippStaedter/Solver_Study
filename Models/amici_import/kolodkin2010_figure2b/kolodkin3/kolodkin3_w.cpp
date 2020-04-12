#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_kolodkin3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = NRLn*RE*k1 - REL*k2;
    w[1] = (-Kapb*Ln_ + Kapf*Lc)/Vnucleus;
    w[2] = Ln_*NRn*k12 - NRLn*k22;
    w[3] = Vnucleus*(Ln_*NRm*k13 - NRLm*k23)/Vcytosol;
    w[4] = -Kapb4*NRm + Kapf4*NRc;
    w[5] = -Kapb5*NRLm + Kapf5*NRLc;
    w[6] = Lc*NRc*k16 - NRLc*k26;
}