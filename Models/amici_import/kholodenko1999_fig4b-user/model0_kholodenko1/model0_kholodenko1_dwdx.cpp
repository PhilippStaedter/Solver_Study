#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kholodenko1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = R*k1f;
    dwdx[1] = -RP*k11b;
    dwdx[2] = k12f;
    dwdx[3] = -ShP*k23b;
    dwdx[4] = RShP*k24f;
    dwdx[5] = -SOS*k12b;
    dwdx[6] = RShP*k17f;
    dwdx[7] = ShP*k21f;
    dwdx[8] = RP*k9f;
    dwdx[9] = RP*k5f;
    dwdx[10] = k25f;
    dwdx[11] = -RP*k7b;
    dwdx[12] = -PLCgP*V8/pow(K8 + PLCgP, 2) + V8/(K8 + PLCgP);
    dwdx[13] = -k25b;
    dwdx[14] = EGF*k1f;
    dwdx[15] = -k2b;
    dwdx[16] = k3f;
    dwdx[17] = SOS*k10f;
    dwdx[18] = -k9b;
    dwdx[19] = -k10b;
    dwdx[20] = k11f;
    dwdx[21] = -GS*k11b;
    dwdx[22] = Shc*k13f;
    dwdx[23] = -ShP*k15b;
    dwdx[24] = -ShG*k18b;
    dwdx[25] = -ShGS*k20b;
    dwdx[26] = -k3b;
    dwdx[27] = -RP*V4/pow(K4 + RP, 2) + V4/(K4 + RP);
    dwdx[28] = PLCg*k5f;
    dwdx[29] = -PLCgP*k7b;
    dwdx[30] = Grb*k9f;
    dwdx[31] = -k5b;
    dwdx[32] = ATP*k6f;
    dwdx[33] = -ADP*k6b;
    dwdx[34] = k7f;
    dwdx[35] = -k13b;
    dwdx[36] = k14f;
    dwdx[37] = -k17b;
    dwdx[38] = k18f;
    dwdx[39] = SOS*k19f;
    dwdx[40] = -k19b;
    dwdx[41] = k20f;
    dwdx[42] = -k24b;
    dwdx[43] = -k14b;
    dwdx[44] = k15f;
    dwdx[45] = Grb*k17f;
    dwdx[46] = GS*k24f;
    dwdx[47] = -k1b;
    dwdx[48] = 2*Ra*k2f;
    dwdx[49] = RG*k10f;
    dwdx[50] = -Grb*k12b;
    dwdx[51] = RShG*k19f;
    dwdx[52] = ShG*k22f;
    dwdx[53] = -RP*k18b;
    dwdx[54] = -k21b;
    dwdx[55] = SOS*k22f;
    dwdx[56] = -RP*k20b;
    dwdx[57] = -k22b;
    dwdx[58] = k23f;
    dwdx[59] = -RP*k15b;
    dwdx[60] = -ShP*V16/pow(K16 + ShP, 2) + V16/(K16 + ShP);
    dwdx[61] = Grb*k21f;
    dwdx[62] = -GS*k23b;
    dwdx[63] = RP*k13f;
}