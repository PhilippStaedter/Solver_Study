#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model0_bucher1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 70.422535211267601*flux_r0;
    xdot[1] = -70.422535211267601*flux_r0 - 70.422535211267601*flux_r13 + 70.422535211267601*flux_r19 - 70.422535211267601*flux_r24 + 70.422535211267601*flux_r28 - 70.422535211267601*flux_r8 - 70.422535211267601*flux_r9;
    xdot[2] = 0.5*flux_r13 - 0.5*flux_r19 - 0.5*flux_r25;
    xdot[3] = 70.422535211267601*flux_r1;
    xdot[4] = -70.422535211267601*flux_r1 - 70.422535211267601*flux_r14 + 70.422535211267601*flux_r20 - 70.422535211267601*flux_r6 + 70.422535211267601*flux_r8;
    xdot[5] = 0.5*flux_r14 - 0.5*flux_r20 - 0.5*flux_r26;
    xdot[6] = 70.422535211267601*flux_r2;
    xdot[7] = -70.422535211267601*flux_r15 - 70.422535211267601*flux_r2 + 70.422535211267601*flux_r21 - 70.422535211267601*flux_r7 + 70.422535211267601*flux_r9;
    xdot[8] = 0.5*flux_r15 - 0.5*flux_r21 - 0.5*flux_r27;
    xdot[9] = 70.422535211267601*flux_r3;
    xdot[10] = -70.422535211267601*flux_r10 - 70.422535211267601*flux_r11 - 70.422535211267601*flux_r12 + 70.422535211267601*flux_r18 + 70.422535211267601*flux_r24 - 70.422535211267601*flux_r28 - 70.422535211267601*flux_r3;
    xdot[11] = 0.5*flux_r12 - 0.5*flux_r18 + 0.5*flux_r25;
    xdot[12] = 70.422535211267601*flux_r4;
    xdot[13] = 70.422535211267601*flux_r10 - 70.422535211267601*flux_r16 + 70.422535211267601*flux_r22 - 70.422535211267601*flux_r4 + 70.422535211267601*flux_r6;
    xdot[14] = 0.5*flux_r16 - 0.5*flux_r22 + 0.5*flux_r26;
    xdot[15] = 70.422535211267601*flux_r5;
    xdot[16] = 70.422535211267601*flux_r11 - 70.422535211267601*flux_r17 + 70.422535211267601*flux_r23 - 70.422535211267601*flux_r5 + 70.422535211267601*flux_r7;
    xdot[17] = 0.5*flux_r17 - 0.5*flux_r23 + 0.5*flux_r27;
}