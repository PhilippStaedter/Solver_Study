#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_kholodenko1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = EGF*R*k1f - Ra*k1b;
    w[1] = RG*SOS*k10f - RGS*k10b;
    w[2] = -GS*RP*k11b + RGS*k11f;
    w[3] = GS*k12f - Grb*SOS*k12b;
    w[4] = RP*Shc*k13f - RSh*k13b;
    w[5] = RSh*k14f - RShP*k14b;
    w[6] = -RP*ShP*k15b + RShP*k15f;
    w[7] = ShP*V16/(K16 + ShP);
    w[8] = Grb*RShP*k17f - RShG*k17b;
    w[9] = -RP*ShG*k18b + RShG*k18f;
    w[10] = RShG*SOS*k19f - RShGS*k19b;
    w[11] = -R2*k2b + pow(Ra, 2)*k2f;
    w[12] = -RP*ShGS*k20b + RShGS*k20f;
    w[13] = Grb*ShP*k21f - ShG*k21b;
    w[14] = SOS*ShG*k22f - ShGS*k22b;
    w[15] = -GS*ShP*k23b + ShGS*k23f;
    w[16] = GS*RShP*k24f - RShGS*k24b;
    w[17] = PLCgP*k25f - PLCgl*k25b;
    w[18] = R2*k3f - RP*k3b;
    w[19] = RP*V4/(K4 + RP);
    w[20] = PLCg*RP*k5f - RPLCg*k5b;
    w[21] = -ADP*RPLCgP*k6b + ATP*RPLCg*k6f;
    w[22] = -PLCgP*RP*k7b + RPLCgP*k7f;
    w[23] = PLCgP*V8/(K8 + PLCgP);
    w[24] = Grb*RP*k9f - RG*k9b;
}