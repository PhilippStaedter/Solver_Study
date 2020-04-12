#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kholodenko1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EGF;
    x_rdata[1] = GS;
    x_rdata[2] = Grb;
    x_rdata[3] = PLCg;
    x_rdata[4] = PLCgP;
    x_rdata[5] = PLCgl;
    x_rdata[6] = R;
    x_rdata[7] = R2;
    x_rdata[8] = RG;
    x_rdata[9] = RGS;
    x_rdata[10] = RP;
    x_rdata[11] = RPLCg;
    x_rdata[12] = RPLCgP;
    x_rdata[13] = RSh;
    x_rdata[14] = RShG;
    x_rdata[15] = RShGS;
    x_rdata[16] = RShP;
    x_rdata[17] = Ra;
    x_rdata[18] = SOS;
    x_rdata[19] = ShG;
    x_rdata[20] = ShGS;
    x_rdata[21] = ShP;
    x_rdata[22] = Shc;
}