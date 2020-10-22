#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_kolodkin6(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 1:
            dwdp[1] = -Ln_/Vnucleus;
            break;
        case 2:
            dwdp[4] = -NR/Vcytosol;
            break;
        case 3:
            dwdp[5] = -NRLn/Vcytosol;
            break;
        case 4:
            dwdp[1] = Lc/Vnucleus;
            break;
        case 5:
            dwdp[4] = NRc/Vcytosol;
            break;
        case 6:
            dwdp[5] = NRLc/Vcytosol;
            break;
        case 7:
            dwdp[4] = -(-Kapb4*NR + Kapf4*NRc)/pow(Vcytosol, 2);
            dwdp[5] = -(-Kapb5*NRLn + Kapf5*NRLc)/pow(Vcytosol, 2);
            break;
        case 8:
            dwdp[1] = -(-Kapb*Ln_ + Kapf*Lc)/pow(Vnucleus, 2);
            break;
        case 9:
            dwdp[0] = NRLn*RE;
            break;
        case 10:
            dwdp[2] = Ln_*NR;
            break;
        case 11:
            dwdp[3] = Lc*NRc;
            break;
        case 12:
            dwdp[0] = -REL;
            break;
        case 13:
            dwdp[2] = -NRLn;
            break;
        case 14:
            dwdp[3] = -NRLc;
            break;
    }
}