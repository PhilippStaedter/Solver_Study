#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Singh2006(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*kf0*x1;
            break;
        case 1:
            dwdp[0] = 1.0*IL6*x1;
            break;
        case 2:
            dwdp[1] = 1.0*x2;
            break;
        case 3:
            dwdp[3] = 1.0*x3*x4;
            break;
        case 4:
            dwdp[3] = -1.0*x5;
            break;
        case 5:
            dwdp[2] = -1.0*x2*x5;
            break;
        case 6:
            dwdp[2] = 1.0*x6;
            break;
        case 7:
            dwdp[4] = 1.0*pow(x6, 2);
            break;
        case 8:
            dwdp[5] = 1.0*x7;
            break;
        case 9:
            dwdp[6] = 1.0*x7;
            break;
        case 10:
            dwdp[7] = 1.0*x16;
            dwdp[20] = 1.0*x32;
            break;
        case 11:
            dwdp[8] = 1.0*x8*x9;
            dwdp[19] = 1.0*x30*x9;
            break;
        case 12:
            dwdp[8] = -1.0*x11;
            dwdp[19] = -1.0*x31;
            break;
        case 13:
            dwdp[9] = 1.0*x11;
            break;
        case 14:
            dwdp[10] = 1.0*x10*x8;
            break;
        case 15:
            dwdp[10] = -1.0*x12;
            break;
        case 16:
            dwdp[11] = 1.0*x15*x8;
            dwdp[26] = 1.0*x15*x31;
            break;
        case 17:
            dwdp[11] = -1.0*x16;
            dwdp[26] = -1.0*x32;
            break;
        case 18:
            dwdp[12] = 1.0*x29*x8;
            break;
        case 19:
            dwdp[12] = -1.0*x30;
            break;
        case 20:
            dwdp[13] = 1.0*x39;
            break;
        case 21:
            dwdp[13] = -1.0*x46*x8;
            break;
        case 22:
            dwdp[14] = 1.0*x40;
            break;
        case 23:
            dwdp[14] = -1.0*x45*x8;
            break;
        case 24:
            dwdp[15] = 1.0*x41;
            break;
        case 25:
            dwdp[15] = -1.0*x44*x8;
            break;
        case 26:
            dwdp[16] = 1.0*x18;
            dwdp[25] = 1.0*x19;
            break;
        case 27:
            dwdp[17] = 1.0*x10*x9;
            dwdp[32] = -1.0*x21*x22;
            break;
        case 28:
            dwdp[17] = -1.0*x14;
            dwdp[32] = 1.0*x24;
            break;
        case 29:
            dwdp[18] = 1.0*x22;
            break;
        case 30:
            dwdp[21] = 2.0*pow(x10, 2);
            dwdp[29] = -1.0*pow(x21, 2);
            break;
        case 31:
            dwdp[21] = -2.0*x13;
            dwdp[29] = 1.0*x20;
            break;
        case 32:
            dwdp[22] = 1.0*x10*x17;
            dwdp[23] = 1.0*x13*x17;
            break;
        case 33:
            dwdp[22] = -1.0*x18;
            dwdp[23] = -1.0*x19;
            break;
        case 34:
            dwdp[24] = 1.0*x13;
            break;
        case 35:
            dwdp[27] = 1.0*x46/(Km + x46);
            break;
        case 36:
            dwdp[27] = -1.0*Vm*x46/pow(Km + x46, 2);
            break;
        case 37:
            dwdp[28] = 1.0*x16;
            break;
        case 38:
            dwdp[28] = -1.0*x39;
            break;
        case 41:
            dwdp[30] = 1.0*x20*x23;
            dwdp[31] = 1.0*x21*x23;
            break;
        case 42:
            dwdp[30] = -1.0*x27;
            dwdp[31] = -1.0*x28;
            break;
        case 45:
            dwdp[33] = 1.0*x28;
            dwdp[34] = 1.0*x27;
            break;
        case 46:
            dwdp[35] = 1.0*x20/(k18b + x20);
            break;
        case 47:
            dwdp[35] = -1.0*k18a*x20/pow(k18b + x20, 2);
            break;
        case 48:
            dwdp[36] = 1.0*x25;
            break;
        case 49:
            dwdp[37] = 1.0*x26;
            break;
        case 50:
            dwdp[38] = 1.0*x26;
            break;
        case 51:
            dwdp[39] = 1.0*x29;
            dwdp[40] = 1.0*x32;
            break;
        case 52:
            dwdp[41] = 1.0*x34*x46;
            break;
        case 53:
            dwdp[41] = -1.0*x45;
            break;
        case 54:
            dwdp[42] = 1.0*x38;
            break;
        case 55:
            dwdp[42] = -1.0*x34*x35;
            break;
        case 56:
            dwdp[43] = 1.0*x34*x39;
            break;
        case 57:
            dwdp[43] = -1.0*x40;
            break;
        case 58:
            dwdp[44] = 1.0*x35*x40;
            break;
        case 59:
            dwdp[44] = -1.0*x41;
            break;
        case 60:
            dwdp[45] = 1.0*x35*x45;
            break;
        case 61:
            dwdp[45] = -1.0*x44;
            break;
        case 62:
            dwdp[46] = 1.0*x36*x41;
            break;
        case 63:
            dwdp[46] = -1.0*x42;
            break;
        case 64:
            dwdp[47] = 1.0*x43;
            break;
        case 65:
            dwdp[47] = -1.0*x36*x41;
            break;
        case 66:
            dwdp[48] = 1.0*x42;
            break;
        case 67:
            dwdp[48] = -1.0*x37*x41;
            break;
        case 68:
            dwdp[49] = 1.0*x37*x47;
            break;
        case 69:
            dwdp[49] = -1.0*x48;
            break;
        case 70:
            dwdp[50] = 1.0*x38*x39;
            break;
        case 71:
            dwdp[50] = -1.0*x41;
            break;
        case 72:
            dwdp[51] = 1.0*x44;
            break;
        case 73:
            dwdp[51] = -1.0*x38*x46;
            break;
        case 74:
            dwdp[52] = 1.0*x41*x49;
            break;
        case 75:
            dwdp[52] = -1.0*x43;
            break;
        case 76:
            dwdp[53] = 1.0*x52;
            break;
        case 77:
            dwdp[54] = 1.0*x48;
            break;
        case 78:
            dwdp[54] = -1.0*x49*x51;
            break;
        case 79:
            dwdp[55] = 1.0*x50*x51;
            break;
        case 80:
            dwdp[55] = -1.0*x52;
            break;
        case 81:
            dwdp[56] = 1.0*x51*x53;
            break;
        case 82:
            dwdp[56] = -1.0*x54;
            break;
        case 83:
            dwdp[57] = 1.0*x54;
            break;
        case 84:
            dwdp[58] = 1.0*x51*x55;
            break;
        case 85:
            dwdp[58] = -1.0*x56;
            break;
        case 86:
            dwdp[59] = 1.0*x60;
            break;
        case 87:
            dwdp[60] = 1.0*x58;
            break;
        case 88:
            dwdp[61] = 1.0*x55*x59;
            break;
        case 89:
            dwdp[61] = -1.0*x60;
            break;
        case 90:
            dwdp[62] = 1.0*x56;
            break;
        case 91:
            dwdp[63] = 1.0*x57*x59;
            break;
        case 92:
            dwdp[63] = -1.0*x58;
            break;
        case 93:
            dwdp[64] = 1.0*x57*x61;
            break;
        case 94:
            dwdp[64] = -1.0*x62;
            break;
        case 95:
            dwdp[65] = 1.0*x62;
            break;
        case 96:
            dwdp[66] = 1.0*x64;
            break;
        case 97:
            dwdp[66] = -1.0*x57*x63;
            break;
        case 98:
            dwdp[67] = 1.0*x64;
            break;
        case 99:
            dwdp[68] = 1.0*x68;
            break;
        case 100:
            dwdp[69] = 1.0*x67;
            break;
        case 101:
            dwdp[70] = 1.0*x63*x66;
            break;
        case 102:
            dwdp[70] = -1.0*x68;
            break;
        case 103:
            dwdp[71] = 1.0*x65*x66;
            break;
        case 104:
            dwdp[71] = -1.0*x67;
            break;
    }
}