#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_kholodenko1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = EGF;
    y[1] = GS;
    y[2] = Grb;
    y[3] = PLCg;
    y[4] = PLCgP;
    y[5] = PLCgl;
    y[6] = R;
    y[7] = R2;
    y[8] = RG;
    y[9] = RGS;
    y[10] = RP;
    y[11] = RPLCg;
    y[12] = RPLCgP;
    y[13] = RSh;
    y[14] = RShG;
    y[15] = RShGS;
    y[16] = RShP;
    y[17] = Ra;
    y[18] = SOS;
    y[19] = ShG;
    y[20] = ShGS;
    y[21] = ShP;
    y[22] = Shc;
}