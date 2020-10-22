#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_fisher1(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = 8849557522123.8945;
    dxdotdw[1] = -8849557522123.8945;
    dxdotdw[2] = -3717472118959.1074;
    dxdotdw[3] = 8849557522123.8945;
    dxdotdw[4] = 3717472118959.1074;
    dxdotdw[5] = -3717472118959.1074;
    dxdotdw[6] = 3717472118959.1074;
    dxdotdw[7] = 3717472118959.1074;
    dxdotdw[8] = -3717472118959.1074;
    dxdotdw[9] = 3717472118959.1074;
    dxdotdw[10] = -11152416356877.322;
    dxdotdw[11] = -3717472118959.1074;
    dxdotdw[12] = 8849557522123.8945;
    dxdotdw[13] = -26548672566371.684;
    dxdotdw[14] = -8849557522123.8945;
    dxdotdw[15] = -3717472118959.1074;
    dxdotdw[16] = 8849557522123.8945;
    dxdotdw[17] = -3717472118959.1074;
    dxdotdw[18] = 8849557522123.8945;
    dxdotdw[19] = -3717472118959.1074;
    dxdotdw[20] = 8849557522123.8945;
    dxdotdw[21] = -8849557522123.8945;
    dxdotdw[22] = 8849557522123.8945;
    dxdotdw[23] = -8849557522123.8945;
    dxdotdw[24] = 3717472118959.1074;
    dxdotdw[25] = -8849557522123.8945;
    dxdotdw[26] = 3717472118959.1074;
    dxdotdw[27] = -8849557522123.8945;
    dxdotdw[28] = -8849557522123.8945;
    dxdotdw[29] = 8849557522123.8945;
    dxdotdw[30] = 8849557522123.8945;
    dxdotdw[31] = -8849557522123.8945;
    dxdotdw[32] = 8849557522123.8945;
    dxdotdw[33] = 3717472118959.1074;
    dxdotdw[34] = -8849557522123.8945;
    dxdotdw[35] = -3717472118959.1074;
    dxdotdw[36] = 3717472118959.1074;
    dxdotdw[37] = 3717472118959.1074;
    dxdotdw[38] = -3717472118959.1074;
    dxdotdw[39] = 3717472118959.1074;
}