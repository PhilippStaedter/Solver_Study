#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model2_(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[2] = atp;
            break;
        case 1:
            dwdp[1] = -VMAXLG*fbp*phos*(atot - atp)/(pow(KMLGADP, 2)*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp)) + VMAXLG*fbp*phos*pow(atot - atp, 2)/(pow(KMLGADP, 3)*KMLGP*pow(1 + (atot - atp)/KMLGADP, 2)*(1 + phos/KMLGP)*(KMLGF16P + fbp));
            break;
        case 2:
            dwdp[1] = -VMAXLG*fbp*phos*(atot - atp)/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*pow(KMLGF16P + fbp, 2));
            dwdp[4] = -VMAXGLYC*fbp/pow(KMLGF16P + fbp, 2);
            break;
        case 3:
            dwdp[1] = -VMAXLG*fbp*phos*(atot - atp)/(KMLGADP*pow(KMLGP, 2)*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp)) + VMAXLG*fbp*pow(phos, 2)*(atot - atp)/(KMLGADP*pow(KMLGP, 3)*(1 + (atot - atp)/KMLGADP)*pow(1 + phos/KMLGP, 2)*(KMLGF16P + fbp));
            break;
        case 4:
            dwdp[0] = -VMAXPFK*atp/pow(KMPFKATP + atp*(1 + atp/KiPFKATP), 2);
            break;
        case 5:
            dwdp[0] = VMAXPFK*pow(atp, 3)/(pow(KiPFKATP, 2)*pow(KMPFKATP + atp*(1 + atp/KiPFKATP), 2));
            break;
        case 6:
            dwdp[4] = fbp/(KMLGF16P + fbp);
            break;
        case 7:
            dwdp[1] = fbp*phos*(atot - atp)/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp));
            break;
        case 8:
            dwdp[0] = atp/(KMPFKATP + atp*(1 + atp/KiPFKATP));
            break;
        case 9:
            dwdp[1] = VMAXLG*fbp*phos/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp)) - VMAXLG*fbp*phos*(atot - atp)/(pow(KMLGADP, 2)*KMLGP*pow(1 + (atot - atp)/KMLGADP, 2)*(1 + phos/KMLGP)*(KMLGF16P + fbp));
            break;
        case 10:
            dwdp[3] = pT - phos;
            break;
        case 11:
            dwdp[3] = k6;
            break;
    }
}