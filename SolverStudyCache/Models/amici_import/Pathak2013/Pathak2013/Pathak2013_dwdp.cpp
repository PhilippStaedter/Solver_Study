#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Pathak2013(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*s7;
            break;
        case 1:
            dwdp[0] = 1.0*s1;
            break;
        case 2:
            dwdp[1] = -1.0*s7;
            break;
        case 3:
            dwdp[1] = 1.0*s2;
            break;
        case 4:
            dwdp[2] = -1.0*s8;
            break;
        case 5:
            dwdp[2] = 1.0*s2;
            break;
        case 6:
            dwdp[3] = -1.0*s9;
            break;
        case 7:
            dwdp[3] = 1.0*s2;
            break;
        case 8:
            dwdp[4] = -1.0*s12;
            break;
        case 9:
            dwdp[4] = 1.0*s6;
            break;
        case 10:
            dwdp[5] = -1.0*s11;
            break;
        case 11:
            dwdp[5] = 1.0*s6;
            break;
        case 12:
            dwdp[6] = -1.0*s10;
            break;
        case 13:
            dwdp[6] = 1.0*s6;
            break;
        case 14:
            dwdp[7] = -1.0*s9;
            break;
        case 15:
            dwdp[7] = 1.0*s3;
            break;
        case 16:
            dwdp[8] = -1.0*s7;
            break;
        case 17:
            dwdp[8] = 1.0*s3;
            break;
        case 18:
            dwdp[9] = -1.0*s7;
            break;
        case 19:
            dwdp[9] = 1.0*s4;
            break;
        case 20:
            dwdp[10] = -1.0*s7;
            break;
        case 21:
            dwdp[10] = 1.0*s5;
            break;
        case 22:
            dwdp[11] = -1.0*s14;
            break;
        case 23:
            dwdp[11] = 1.0*s13;
            break;
        case 24:
            dwdp[12] = -1.0*s13;
            break;
        case 25:
            dwdp[12] = 1.0*s7;
            break;
        case 26:
            dwdp[13] = -1.0*s13;
            break;
        case 27:
            dwdp[13] = 1.0*s8;
            break;
        case 28:
            dwdp[14] = -1.0*s13;
            break;
        case 29:
            dwdp[14] = 1.0*s9;
            break;
        case 30:
            dwdp[15] = -1.0*s13;
            break;
        case 31:
            dwdp[15] = 1.0*s10;
            break;
        case 32:
            dwdp[16] = -1.0*s15;
            break;
        case 33:
            dwdp[16] = 1.0*s14;
            break;
        case 34:
            dwdp[17] = -1.0*s15;
            break;
        case 35:
            dwdp[17] = 1.0*s7;
            break;
        case 36:
            dwdp[18] = -1.0*s16;
            break;
        case 37:
            dwdp[18] = 1.0*s14;
            break;
        case 38:
            dwdp[19] = -1.0*s16;
            break;
        case 39:
            dwdp[19] = 1.0*s11;
            break;
        case 40:
            dwdp[20] = -1.0*s16;
            break;
        case 41:
            dwdp[20] = 1.0*s12;
            break;
        case 42:
            dwdp[21] = -1.0*s18;
            break;
        case 43:
            dwdp[21] = 1.0*s17;
            break;
        case 44:
            dwdp[22] = -1.0*s17;
            break;
        case 45:
            dwdp[22] = 1.0*s14;
            break;
        case 46:
            dwdp[23] = -1.0*s19;
            break;
        case 47:
            dwdp[23] = 1.0*s18;
            break;
        case 48:
            dwdp[24] = -1.0*s20;
            break;
        case 49:
            dwdp[24] = 1.0*s18;
            break;
        case 50:
            dwdp[25] = -1.0*s21;
            break;
        case 51:
            dwdp[25] = 1.0*s18;
            break;
        case 52:
            dwdp[26] = -1.0*s22;
            break;
        case 53:
            dwdp[26] = 1.0*s18;
            break;
        case 54:
            dwdp[27] = -1.0*s23;
            break;
        case 55:
            dwdp[27] = 1.0*s18;
            break;
        case 56:
            dwdp[28] = -1.0*s24;
            break;
        case 57:
            dwdp[28] = 1.0*s18;
            break;
        case 58:
            dwdp[29] = -1.0*s25;
            break;
        case 59:
            dwdp[29] = 1.0*s18;
            break;
        case 60:
            dwdp[30] = -1.0*s26;
            break;
        case 61:
            dwdp[30] = 1.0*s18;
            break;
        case 62:
            dwdp[31] = -1.0*s28;
            break;
        case 63:
            dwdp[31] = 1.0*s27;
            break;
        case 64:
            dwdp[32] = -1.0*s27;
            break;
        case 65:
            dwdp[32] = 1.0*s18;
            break;
        case 66:
            dwdp[33] = -1.0*s19;
            break;
        case 67:
            dwdp[33] = 1.0*s15;
            break;
        case 68:
            dwdp[34] = -1.0*s20;
            break;
        case 69:
            dwdp[34] = 1.0*s15;
            break;
        case 70:
            dwdp[35] = -1.0*s26;
            break;
        case 71:
            dwdp[35] = 1.0*s16;
            break;
        case 72:
            dwdp[36] = -1.0*s29;
            break;
        case 73:
            dwdp[36] = 1.0*s28;
            break;
        case 74:
            dwdp[37] = -1.0*s30;
            break;
        case 75:
            dwdp[37] = 1.0*s28;
            break;
        case 76:
            dwdp[38] = -1.0*s31;
            break;
        case 77:
            dwdp[38] = 1.0*s28;
            break;
        case 78:
            dwdp[39] = -1.0*s32;
            break;
        case 79:
            dwdp[39] = 1.0*s28;
            break;
        case 80:
            dwdp[40] = -1.0*s30;
            break;
        case 81:
            dwdp[40] = 1.0*s20;
            break;
        case 82:
            dwdp[41] = -1.0*s31;
            break;
        case 83:
            dwdp[41] = 1.0*s20;
            break;
        case 84:
            dwdp[42] = -1.0*s32;
            break;
        case 85:
            dwdp[42] = 1.0*s20;
            break;
        case 86:
            dwdp[43] = -1.0*s30;
            break;
        case 87:
            dwdp[43] = 1.0*s26;
            break;
        case 88:
            dwdp[44] = -1.0*s34;
            break;
        case 89:
            dwdp[44] = 1.0*s33;
            break;
        case 90:
            dwdp[45] = -1.0*s36;
            break;
        case 91:
            dwdp[45] = 1.0*s35;
            break;
        case 92:
            dwdp[46] = -1.0*s38;
            break;
        case 93:
            dwdp[46] = 1.0*s37;
            break;
        case 94:
            dwdp[47] = -1.0*s40;
            break;
        case 95:
            dwdp[47] = 1.0*s39;
            break;
        case 96:
            dwdp[48] = -1.0*s42;
            break;
        case 97:
            dwdp[48] = 1.0*s41;
            break;
        case 98:
            dwdp[49] = -1.0*s44;
            break;
        case 99:
            dwdp[49] = 1.0*s43;
            break;
        case 100:
            dwdp[50] = -1.0*s46;
            break;
        case 101:
            dwdp[50] = 1.0*s45;
            break;
        case 102:
            dwdp[51] = -1.0*s48;
            break;
        case 103:
            dwdp[51] = 1.0*s47;
            break;
        case 104:
            dwdp[52] = -1.0*s50;
            break;
        case 105:
            dwdp[52] = 1.0*s49;
            break;
        case 106:
            dwdp[53] = -1.0*s52;
            break;
        case 107:
            dwdp[53] = 1.0*s51;
            break;
        case 108:
            dwdp[54] = -1.0*s37;
            break;
        case 109:
            dwdp[54] = 1.0*s29;
            break;
        case 110:
            dwdp[55] = -1.0*s33;
            break;
        case 111:
            dwdp[55] = 1.0*s29;
            break;
        case 112:
            dwdp[56] = -1.0*s35;
            break;
        case 113:
            dwdp[56] = 1.0*s30;
            break;
        case 114:
            dwdp[57] = -1.0*s41;
            break;
        case 115:
            dwdp[57] = 1.0*s30;
            break;
        case 116:
            dwdp[58] = -1.0*s47;
            break;
        case 117:
            dwdp[58] = 1.0*s30;
            break;
        case 118:
            dwdp[59] = -1.0*s33;
            break;
        case 119:
            dwdp[59] = 1.0*s31;
            break;
        case 120:
            dwdp[60] = -1.0*s45;
            break;
        case 121:
            dwdp[60] = 1.0*s31;
            break;
        case 122:
            dwdp[61] = -1.0*s39;
            break;
        case 123:
            dwdp[61] = 1.0*s31;
            break;
        case 124:
            dwdp[62] = -1.0*s47;
            break;
        case 125:
            dwdp[62] = 1.0*s32;
            break;
        case 126:
            dwdp[63] = -1.0*s45;
            break;
        case 127:
            dwdp[63] = 1.0*s32;
            break;
        case 128:
            dwdp[64] = -1.0*s35;
            break;
        case 129:
            dwdp[64] = 1.0*s32;
            break;
        case 130:
            dwdp[65] = -1.0*s56;
            break;
        case 131:
            dwdp[65] = 1.0*s28;
            break;
        case 132:
            dwdp[66] = -1.0*s49;
            break;
        case 133:
            dwdp[66] = 1.0*s28;
            break;
        case 134:
            dwdp[67] = -1.0*s51;
            break;
        case 135:
            dwdp[67] = 1.0*s28;
            break;
        case 136:
            dwdp[68] = -1.0*s53;
            break;
        case 137:
            dwdp[68] = 1.0*s28;
            break;
        case 138:
            dwdp[69] = -1.0*s54;
            break;
        case 139:
            dwdp[69] = 1.0*s28;
            break;
        case 140:
            dwdp[70] = -1.0*s55;
            break;
        case 141:
            dwdp[70] = 1.0*s28;
            break;
        case 142:
            dwdp[71] = -1.0*s57;
            break;
        case 143:
            dwdp[71] = 1.0*s40;
            break;
        case 144:
            dwdp[72] = -1.0*s57;
            break;
        case 145:
            dwdp[72] = 1.0*s53;
            break;
        case 146:
            dwdp[73] = -1.0*s57;
            break;
        case 147:
            dwdp[73] = 1.0*s54;
            break;
        case 148:
            dwdp[74] = -1.0*s57;
            break;
        case 149:
            dwdp[74] = 1.0*s52;
            break;
        case 150:
            dwdp[75] = -1.0*s57;
            break;
        case 151:
            dwdp[75] = 1.0*s50;
            break;
        case 152:
            dwdp[76] = -1.0*s57;
            break;
        case 153:
            dwdp[76] = 1.0*s56;
            break;
        case 154:
            dwdp[77] = -1.0*s57;
            break;
        case 155:
            dwdp[77] = 1.0*s48;
            break;
        case 156:
            dwdp[78] = -1.0*s43;
            break;
        case 157:
            dwdp[78] = 1.0*s30;
            break;
        case 158:
            dwdp[79] = -1.0*s57;
            break;
        case 159:
            dwdp[79] = 1.0*s55;
            break;
        case 160:
            dwdp[80] = -1.0*s57;
            break;
        case 161:
            dwdp[80] = 1.0*s42;
            break;
        case 162:
            dwdp[81] = -1.0*s57;
            break;
        case 163:
            dwdp[81] = 1.0*s44;
            break;
        case 164:
            dwdp[82] = -1.0*s57;
            break;
        case 165:
            dwdp[82] = 1.0*s38;
            break;
        case 166:
            dwdp[83] = -1.0*s57;
            break;
        case 167:
            dwdp[83] = 1.0*s36;
            break;
        case 168:
            dwdp[84] = -1.0*s57;
            break;
        case 169:
            dwdp[84] = 1.0*s34;
            break;
        case 170:
            dwdp[85] = -1.0*s57;
            break;
        case 171:
            dwdp[85] = 1.0*s46;
            break;
    }
}