#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_kolodkin1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 1:
            dwdp[1] = -L/Vnucleus;
            break;
        case 2:
            dwdp[1] = Lc/Vnucleus;
            break;
        case 3:
            dwdp[1] = Kapf/Vnucleus;
            break;
        case 5:
            dwdp[1] = -(-Kapb*L + Kapf*Lc)/pow(Vnucleus, 2);
            break;
        case 6:
            dwdp[0] = L*ReNR;
            break;
        case 7:
            dwdp[0] = -REL;
            break;
    }
}