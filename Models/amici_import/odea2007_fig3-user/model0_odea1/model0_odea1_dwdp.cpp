#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_odea1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*IKK;
            break;
        case 1:
            dwdp[1] = 1.0*IkBa_cytoplasm;
            break;
        case 2:
            dwdp[2] = 1.0*IkBa_nucleus;
            break;
        case 3:
            dwdp[3] = 1.0*IkBaIKK;
            break;
        case 4:
            dwdp[4] = 1.0*IkBaIKKNFkB;
            break;
        case 5:
            dwdp[5] = 1.0*IkBaNFkB_cytoplasm;
            break;
        case 6:
            dwdp[6] = 1.0*IkBaNFkB_nucleus;
            break;
        case 7:
            dwdp[7] = 1.0*IkBb_cytoplasm;
            break;
        case 8:
            dwdp[8] = 1.0*IkBb_nucleus;
            break;
        case 9:
            dwdp[9] = 1.0*IkBbIKK;
            break;
        case 10:
            dwdp[10] = 1.0*IkBbIKKNFkB;
            break;
        case 11:
            dwdp[11] = 1.0*IkBbNFkB_cytoplasm;
            break;
        case 12:
            dwdp[12] = 1.0*IkBbNFkB_nucleus;
            break;
        case 13:
            dwdp[13] = 1.0*IkBe_cytoplasm;
            break;
        case 14:
            dwdp[14] = 1.0*IkBe_nucleus;
            break;
        case 15:
            dwdp[15] = 1.0*IkBeIKK;
            break;
        case 16:
            dwdp[16] = 1.0*IkBeIKKNFkB;
            break;
        case 17:
            dwdp[17] = 1.0*IkBeNFkB_cytoplasm;
            break;
        case 18:
            dwdp[18] = 1.0*IkBeNFkB_nucleus;
            break;
        case 19:
            dwdp[19] = -1.0*IkBaIKKNFkB;
            break;
        case 20:
            dwdp[19] = 1.0*IkBaIKK*NFkB_cytoplasm;
            break;
        case 21:
            dwdp[20] = -1.0*IkBaIKKNFkB;
            break;
        case 22:
            dwdp[20] = 1.0*IKK*IkBaNFkB_cytoplasm;
            break;
        case 23:
            dwdp[21] = -1.0*IkBbIKKNFkB;
            break;
        case 24:
            dwdp[21] = 1.0*IkBbIKK*NFkB_cytoplasm;
            break;
        case 25:
            dwdp[22] = -1.0*IkBbIKKNFkB;
            break;
        case 26:
            dwdp[22] = 1.0*IKK*IkBbNFkB_cytoplasm;
            break;
        case 27:
            dwdp[23] = -1.0*IkBeIKKNFkB;
            break;
        case 28:
            dwdp[23] = 1.0*IkBeIKK*NFkB_cytoplasm;
            break;
        case 29:
            dwdp[24] = -1.0*IkBeIKKNFkB;
            break;
        case 30:
            dwdp[24] = 1.0*IKK*IkBeNFkB_cytoplasm;
            break;
        case 31:
            dwdp[25] = -1.0*IkBaIKK;
            break;
        case 32:
            dwdp[25] = 1.0*IKK*IkBa_cytoplasm;
            break;
        case 33:
            dwdp[26] = -1.0*IkBaNFkB_cytoplasm;
            break;
        case 34:
            dwdp[26] = 1.0*IkBa_cytoplasm*NFkB_cytoplasm;
            break;
        case 35:
            dwdp[27] = -1.0*IkBaNFkB_nucleus;
            break;
        case 36:
            dwdp[27] = 1.0*IkBa_nucleus*NFkB_nucleus;
            break;
        case 37:
            dwdp[28] = -1.0*IkBbIKK;
            break;
        case 38:
            dwdp[28] = 1.0*IKK*IkBb_cytoplasm;
            break;
        case 39:
            dwdp[29] = -1.0*IkBbNFkB_cytoplasm;
            break;
        case 40:
            dwdp[29] = 1.0*IkBb_cytoplasm*NFkB_cytoplasm;
            break;
        case 41:
            dwdp[30] = -1.0*IkBbNFkB_nucleus;
            break;
        case 42:
            dwdp[30] = 1.0*IkBb_nucleus*NFkB_nucleus;
            break;
        case 43:
            dwdp[31] = -1.0*IkBeIKK;
            break;
        case 44:
            dwdp[31] = 1.0*IKK*IkBe_cytoplasm;
            break;
        case 45:
            dwdp[32] = -1.0*IkBeNFkB_cytoplasm;
            break;
        case 46:
            dwdp[32] = 1.0*IkBe_cytoplasm*NFkB_cytoplasm;
            break;
        case 47:
            dwdp[33] = -1.0*IkBeNFkB_nucleus;
            break;
        case 48:
            dwdp[33] = 1.0*IkBe_nucleus*NFkB_nucleus;
            break;
        case 49:
            dwdp[34] = 1.0*pow(NFkB_nucleus, 2);
            break;
        case 50:
            dwdp[35] = -1.0*IkBa_nucleus;
            break;
        case 51:
            dwdp[35] = 1.0*IkBa_cytoplasm;
            break;
        case 52:
            dwdp[36] = 1.0*IkBaNFkB_nucleus;
            break;
        case 53:
            dwdp[37] = -1.0*IkBb_nucleus;
            break;
        case 54:
            dwdp[37] = 1.0*IkBb_cytoplasm;
            break;
        case 55:
            dwdp[38] = 1.0*IkBbNFkB_nucleus;
            break;
        case 56:
            dwdp[39] = -1.0*IkBe_nucleus;
            break;
        case 57:
            dwdp[39] = 1.0*IkBe_cytoplasm;
            break;
        case 58:
            dwdp[40] = 1.0*IkBeNFkB_nucleus;
            break;
        case 59:
            dwdp[41] = -1.0*NFkB_nucleus;
            break;
        case 60:
            dwdp[41] = 1.0*NFkB_cytoplasm;
            break;
        case 61:
            dwdp[42] = 1.0*IkBa_mRNA;
            break;
        case 62:
            dwdp[43] = 1.0*IkBb_mRNA;
            break;
        case 63:
            dwdp[44] = 1.0*IkBe_mRNA;
            break;
        case 64:
            dwdp[45] = 1.0*IkBa_mRNA;
            break;
        case 65:
            dwdp[46] = 1.0*IkBb_mRNA;
            break;
        case 66:
            dwdp[47] = 1.0*IkBe_mRNA;
            break;
        case 67:
            dwdp[48] = 1.0;
            break;
        case 68:
            dwdp[49] = 1.0;
            break;
        case 69:
            dwdp[50] = 1.0;
            break;
    }
}