#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Ueda2001(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[1] = 1.0*flux_r10 - 1.0*flux_r18 - 1.0*flux_r30 - 1.0*flux_r6 + 1.0*flux_r7;
    xdot[2] = -1.0*flux_r20 - 1.0*flux_r31 + 1.0*flux_r6 - 1.0*flux_r7;
    xdot[3] = -1.0*flux_r10 + 1.0*flux_r13 - 1.0*flux_r19 - 1.0*flux_r29;
    xdot[4] = -1.0*flux_r28 + 1.0*flux_r4 - 1.0*flux_r5;
    xdot[5] = -1.0*flux_r11 + 1.0*flux_r14 - 1.0*flux_r15 - 1.0*flux_r23;
    xdot[6] = 1.0*flux_r0 - 1.0*flux_r1 - 1.0*flux_r22;
    xdot[7] = 1.0*flux_r11 - 1.0*flux_r16 - 1.0*flux_r26 + 1.0*flux_r8 - 1.0*flux_r9;
    xdot[8] = -1.0*flux_r17 - 1.0*flux_r27 - 1.0*flux_r8 + 1.0*flux_r9;
    xdot[9] = -1.0*flux_r11 + 1.0*flux_r12 - 1.0*flux_r21 - 1.0*flux_r25;
    xdot[10] = 1.0*flux_r2 - 1.0*flux_r24 - 1.0*flux_r3;
}