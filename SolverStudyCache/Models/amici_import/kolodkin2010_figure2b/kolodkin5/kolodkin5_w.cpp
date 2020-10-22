#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_kolodkin5(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = NRLn*RE*k1 - REL*k2;
    w[1] = (-Kapb*Ln_ + Kapf*Lc)/Vnucleus;
    w[2] = Ln_*NR*k12 - NRLn*k22;
    w[3] = Lc*NRc*k13 - NRLc*k23;
    w[4] = (-Kapb4*NR + Kapf4*NRc)/Vcytosol;
    w[5] = (-Kapb5*NRLn + Kapf5*NRLc)/Vcytosol;
}