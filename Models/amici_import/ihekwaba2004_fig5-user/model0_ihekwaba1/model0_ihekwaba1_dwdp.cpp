#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_ihekwaba1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[17] = IkBa*NFkB;
            break;
        case 1:
            dwdp[61] = IKKIkBb*NFkB;
            break;
        case 2:
            dwdp[63] = IKKIkBbNFkB;
            break;
        case 3:
            dwdp[62] = IKKIkBbNFkB;
            break;
        case 4:
            dwdp[4] = IKKIkBe*NFkB;
            break;
        case 5:
            dwdp[6] = IKKIkBeNFkB;
            break;
        case 6:
            dwdp[5] = IKKIkBeNFkB;
            break;
        case 7:
            dwdp[24] = IkBaNFkB;
            break;
        case 8:
            dwdp[35] = IkBbNFkB;
            break;
        case 9:
            dwdp[46] = IkBeNFkB;
            break;
        case 10:
            dwdp[52] = NFkB;
            break;
        case 11:
            dwdp[19] = IkBaNFkB;
            break;
        case 12:
            dwdp[51] = NFkBn;
            break;
        case 13:
            dwdp[20] = IkBan*NFkBn;
            break;
        case 14:
            dwdp[21] = IkBanNFkBn;
            break;
        case 15:
            dwdp[31] = IkBbn*NFkBn;
            break;
        case 16:
            dwdp[32] = IkBbnNFkBn;
            break;
        case 17:
            dwdp[42] = IkBen*NFkBn;
            break;
        case 18:
            dwdp[43] = IkBenNFkBn;
            break;
        case 19:
            dwdp[53] = source;
            break;
        case 20:
            dwdp[18] = pow(NFkBn, 2);
            break;
        case 21:
            dwdp[25] = IkBat;
            break;
        case 22:
            dwdp[29] = IkBb*NFkB;
            break;
        case 23:
            dwdp[54] = source;
            break;
        case 24:
            dwdp[36] = IkBbt;
            break;
        case 25:
            dwdp[56] = source;
            break;
        case 26:
            dwdp[47] = IkBet;
            break;
        case 27:
            dwdp[33] = IKK*IkBa;
            break;
        case 28:
            dwdp[55] = IKKIkBa;
            break;
        case 29:
            dwdp[28] = IkBat;
            break;
        case 30:
            dwdp[57] = IkBa;
            break;
        case 31:
            dwdp[27] = IkBa;
            break;
        case 32:
            dwdp[26] = IkBan;
            break;
        case 33:
            dwdp[30] = IkBbNFkB;
            break;
        case 34:
            dwdp[1] = IKK*IkBb;
            break;
        case 35:
            dwdp[3] = IKKIkBb;
            break;
        case 36:
            dwdp[39] = IkBbt;
            break;
        case 37:
            dwdp[58] = IkBb;
            break;
        case 38:
            dwdp[38] = IkBb;
            break;
        case 39:
            dwdp[37] = IkBbn;
            break;
        case 40:
            dwdp[7] = IKK*IkBe;
            break;
        case 41:
            dwdp[9] = IKKIkBe;
            break;
        case 42:
            dwdp[50] = IkBet;
            break;
        case 43:
            dwdp[59] = IkBe;
            break;
        case 44:
            dwdp[40] = IkBe*NFkB;
            break;
        case 45:
            dwdp[49] = IkBe;
            break;
        case 46:
            dwdp[48] = IkBen;
            break;
        case 47:
            dwdp[10] = IKK*IkBaNFkB;
            break;
        case 48:
            dwdp[12] = IKKIkBaNFkB;
            break;
        case 49:
            dwdp[23] = IkBanNFkBn;
            break;
        case 50:
            dwdp[13] = IKK*IkBbNFkB;
            break;
        case 51:
            dwdp[14] = IKKIkBbNFkB;
            break;
        case 52:
            dwdp[34] = IkBbnNFkBn;
            break;
        case 53:
            dwdp[15] = IKK*IkBeNFkB;
            break;
        case 54:
            dwdp[16] = IKKIkBeNFkB;
            break;
        case 55:
            dwdp[41] = IkBeNFkB;
            break;
        case 56:
            dwdp[45] = IkBenNFkBn;
            break;
        case 57:
            dwdp[60] = IKK;
            break;
        case 58:
            dwdp[44] = IKKIkBa;
            break;
        case 59:
            dwdp[2] = IKKIkBb;
            break;
        case 60:
            dwdp[8] = IKKIkBe;
            break;
        case 61:
            dwdp[0] = IKKIkBa*NFkB;
            break;
        case 62:
            dwdp[22] = IKKIkBaNFkB;
            break;
        case 63:
            dwdp[11] = IKKIkBaNFkB;
            break;
    }
}