#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_kolodkin3(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = Kapf/Vnucleus;
    dwdx[1] = NRc*k16;
    dwdx[2] = -Kapb/Vnucleus;
    dwdx[3] = NRn*k12;
    dwdx[4] = NRm*Vnucleus*k13/Vcytosol;
    dwdx[5] = Kapf5;
    dwdx[6] = -k26;
    dwdx[7] = -Vnucleus*k23/Vcytosol;
    dwdx[8] = -Kapb5;
    dwdx[9] = RE*k1;
    dwdx[10] = -k22;
    dwdx[11] = Kapf4;
    dwdx[12] = Lc*k16;
    dwdx[13] = Ln_*Vnucleus*k13/Vcytosol;
    dwdx[14] = -Kapb4;
    dwdx[15] = Ln_*k12;
    dwdx[16] = NRLn*k1;
    dwdx[17] = -k2;
}