#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = VMAXPFK*atp*(-1 - 2*atp/KiPFKATP)/pow(KMPFKATP + atp*(1 + atp/KiPFKATP), 2) + VMAXPFK/(KMPFKATP + atp*(1 + atp/KiPFKATP));
    dwdx[1] = -VMAXLG*fbp*phos/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp)) + VMAXLG*fbp*phos*(atot - atp)/(pow(KMLGADP, 2)*KMLGP*pow(1 + (atot - atp)/KMLGADP, 2)*(1 + phos/KMLGP)*(KMLGF16P + fbp));
    dwdx[2] = KATPASE;
    dwdx[3] = -VMAXLG*fbp*phos*(atot - atp)/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*pow(KMLGF16P + fbp, 2)) + VMAXLG*phos*(atot - atp)/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp));
    dwdx[4] = -VMAXGLYC*fbp/pow(KMLGF16P + fbp, 2) + VMAXGLYC/(KMLGF16P + fbp);
    dwdx[5] = VMAXLG*fbp*(atot - atp)/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp)) - VMAXLG*fbp*phos*(atot - atp)/(KMLGADP*pow(KMLGP, 2)*(1 + (atot - atp)/KMLGADP)*pow(1 + phos/KMLGP, 2)*(KMLGF16P + fbp));
    dwdx[6] = -k6;
}