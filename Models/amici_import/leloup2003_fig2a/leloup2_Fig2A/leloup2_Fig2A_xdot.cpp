#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_leloup2_Fig2A(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 1.0*flux_r0 - 1.0*flux_r12 - 1.0*flux_r17;
    xdot[1] = 1.0*flux_r1 - 1.0*flux_r23 - 1.0*flux_r28 - 1.0*flux_r34 + 1.0*flux_r36;
    xdot[2] = -1.0*flux_r26 + 1.0*flux_r42 + 1.0*flux_r5 - 1.0*flux_r6 - 1.0*flux_r8;
    xdot[3] = -1.0*flux_r25 - 1.0*flux_r42 + 1.0*flux_r6 - 1.0*flux_r7;
    xdot[4] = -1.0*flux_r14 - 1.0*flux_r16 + 1.0*flux_r40 - 1.0*flux_r41 + 1.0*flux_r8;
    xdot[5] = -1.0*flux_r11 + 1.0*flux_r30 - 1.0*flux_r31 - 1.0*flux_r37;
    xdot[6] = 1.0*flux_r14 - 1.0*flux_r27 - 1.0*flux_r40 - 1.0*flux_r9;
    xdot[7] = 1.0*flux_r16 - 1.0*flux_r32 - 1.0*flux_r33;
    xdot[8] = -1.0*flux_r10 + 1.0*flux_r23 - 1.0*flux_r29 - 1.0*flux_r36;
    xdot[9] = -1.0*flux_r16 - 1.0*flux_r30 + 1.0*flux_r34 - 1.0*flux_r35 + 1.0*flux_r37;
    xdot[10] = -1.0*flux_r21 + 1.0*flux_r38 + 1.0*flux_r43 - 1.0*flux_r45 - 1.0*flux_r5;
    xdot[11] = 1.0*flux_r13 - 1.0*flux_r18 - 1.0*flux_r44;
    xdot[12] = -1.0*flux_r24 - 1.0*flux_r38 + 1.0*flux_r45 - 1.0*flux_r46;
    xdot[13] = -1.0*flux_r15 - 1.0*flux_r19 + 1.0*flux_r47;
    xdot[14] = 1.0*flux_r2 - 1.0*flux_r20 + 1.0*flux_r39 - 1.0*flux_r4 - 1.0*flux_r5;
    xdot[15] = -1.0*flux_r22 - 1.0*flux_r3 - 1.0*flux_r39 + 1.0*flux_r4;
}