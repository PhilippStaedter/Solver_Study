#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model2_(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = VMAXPFK*atp/(KMPFKATP + atp*(1 + atp/KiPFKATP));
    w[1] = VMAXLG*fbp*phos*(atot - atp)/(KMLGADP*KMLGP*(1 + (atot - atp)/KMLGADP)*(1 + phos/KMLGP)*(KMLGF16P + fbp));
    w[2] = KATPASE*atp;
    w[3] = k6*(pT - phos);
    w[4] = VMAXGLYC*fbp/(KMLGF16P + fbp);
}