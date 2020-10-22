#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kholodenko1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[21] = -RPLCgP*k6b;
            break;
        case 1:
            dwdp[21] = RPLCg*k6f;
            break;
        case 2:
            dwdp[7] = -ShP*V16/pow(K16 + ShP, 2);
            break;
        case 3:
            dwdp[19] = -RP*V4/pow(K4 + RP, 2);
            break;
        case 4:
            dwdp[23] = -PLCgP*V8/pow(K8 + PLCgP, 2);
            break;
        case 5:
            dwdp[7] = ShP/(K16 + ShP);
            break;
        case 6:
            dwdp[19] = RP/(K4 + RP);
            break;
        case 7:
            dwdp[23] = PLCgP/(K8 + PLCgP);
            break;
        case 8:
            dwdp[1] = -RGS;
            break;
        case 9:
            dwdp[1] = RG*SOS;
            break;
        case 10:
            dwdp[2] = -GS*RP;
            break;
        case 11:
            dwdp[2] = RGS;
            break;
        case 12:
            dwdp[3] = -Grb*SOS;
            break;
        case 13:
            dwdp[3] = GS;
            break;
        case 14:
            dwdp[4] = -RSh;
            break;
        case 15:
            dwdp[4] = RP*Shc;
            break;
        case 16:
            dwdp[5] = -RShP;
            break;
        case 17:
            dwdp[5] = RSh;
            break;
        case 18:
            dwdp[6] = -RP*ShP;
            break;
        case 19:
            dwdp[6] = RShP;
            break;
        case 20:
            dwdp[8] = -RShG;
            break;
        case 21:
            dwdp[8] = Grb*RShP;
            break;
        case 22:
            dwdp[9] = -RP*ShG;
            break;
        case 23:
            dwdp[9] = RShG;
            break;
        case 24:
            dwdp[10] = -RShGS;
            break;
        case 25:
            dwdp[10] = RShG*SOS;
            break;
        case 26:
            dwdp[0] = -Ra;
            break;
        case 27:
            dwdp[0] = EGF*R;
            break;
        case 28:
            dwdp[12] = -RP*ShGS;
            break;
        case 29:
            dwdp[12] = RShGS;
            break;
        case 30:
            dwdp[13] = -ShG;
            break;
        case 31:
            dwdp[13] = Grb*ShP;
            break;
        case 32:
            dwdp[14] = -ShGS;
            break;
        case 33:
            dwdp[14] = SOS*ShG;
            break;
        case 34:
            dwdp[15] = -GS*ShP;
            break;
        case 35:
            dwdp[15] = ShGS;
            break;
        case 36:
            dwdp[16] = -RShGS;
            break;
        case 37:
            dwdp[16] = GS*RShP;
            break;
        case 38:
            dwdp[17] = -PLCgl;
            break;
        case 39:
            dwdp[17] = PLCgP;
            break;
        case 40:
            dwdp[11] = -R2;
            break;
        case 41:
            dwdp[11] = pow(Ra, 2);
            break;
        case 42:
            dwdp[18] = -RP;
            break;
        case 43:
            dwdp[18] = R2;
            break;
        case 44:
            dwdp[20] = -RPLCg;
            break;
        case 45:
            dwdp[20] = PLCg*RP;
            break;
        case 46:
            dwdp[21] = -ADP*RPLCgP;
            break;
        case 47:
            dwdp[21] = ATP*RPLCg;
            break;
        case 48:
            dwdp[22] = -PLCgP*RP;
            break;
        case 49:
            dwdp[22] = RPLCgP;
            break;
        case 50:
            dwdp[24] = -RG;
            break;
        case 51:
            dwdp[24] = Grb*RP;
            break;
    }
}