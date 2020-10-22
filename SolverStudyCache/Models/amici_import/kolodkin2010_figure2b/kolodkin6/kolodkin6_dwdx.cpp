#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_kolodkin6(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = Kapf/Vnucleus;
    dwdx[1] = NRc*k13;
    dwdx[2] = -Kapb/Vnucleus;
    dwdx[3] = NR*k12;
    dwdx[4] = Ln_*k12;
    dwdx[5] = -Kapb4/Vcytosol;
    dwdx[6] = -k23;
    dwdx[7] = Kapf5/Vcytosol;
    dwdx[8] = RE*k1;
    dwdx[9] = -k22;
    dwdx[10] = -Kapb5/Vcytosol;
    dwdx[11] = Lc*k13;
    dwdx[12] = Kapf4/Vcytosol;
    dwdx[13] = NRLn*k1;
    dwdx[14] = -k2;
}