#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_fisher1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 3717472118959.1074*flux_r11 + 3717472118959.1074*flux_r16 + 3717472118959.1074*flux_r2 + 3717472118959.1074*flux_r4;
    xdot[1] = -8849557522123.8945*flux_r11 + 8849557522123.8945*flux_r13 + 8849557522123.8945*flux_r5 - 8849557522123.8945*flux_r9;
    xdot[2] = -11152416356877.322*flux_r4 - 3717472118959.1074*flux_r7;
    xdot[3] = -26548672566371.684*flux_r5 + 8849557522123.8945*flux_r7;
    xdot[4] = -3717472118959.1074*flux_r4 - 3717472118959.1074*flux_r6;
    xdot[5] = -8849557522123.8945*flux_r5 + 8849557522123.8945*flux_r6;
    xdot[6] = 3717472118959.1074*flux_r14 - 3717472118959.1074*flux_r15 - 3717472118959.1074*flux_r2;
    xdot[7] = -8849557522123.8945*flux_r12 - 8849557522123.8945*flux_r14 + 8849557522123.8945*flux_r9;
    xdot[8] = 3717472118959.1074*flux_r10 + 3717472118959.1074*flux_r2 + 3717472118959.1074*flux_r3;
    xdot[9] = 8849557522123.8945*flux_r0 - 8849557522123.8945*flux_r10 - 8849557522123.8945*flux_r9;
    xdot[10] = 3717472118959.1074*flux_r15 - 3717472118959.1074*flux_r16 - 3717472118959.1074*flux_r8;
    xdot[11] = 8849557522123.8945*flux_r12 - 8849557522123.8945*flux_r13 + 8849557522123.8945*flux_r8;
    xdot[12] = -3717472118959.1074*flux_r1 + 3717472118959.1074*flux_r16 - 3717472118959.1074*flux_r3;
    xdot[13] = -8849557522123.8945*flux_r0 + 8849557522123.8945*flux_r1 + 8849557522123.8945*flux_r13;
}