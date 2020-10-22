#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Ung2008(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*species_2;
            break;
        case 1:
            dwdp[0] = 1.0*species_0*species_1;
            break;
        case 2:
            dwdp[1] = -1.0*species_3;
            break;
        case 3:
            dwdp[1] = 1.0*pow(species_2, 2);
            break;
        case 4:
            dwdp[2] = 1.0*species_3;
            break;
        case 5:
            dwdp[3] = -1.0*species_6;
            break;
        case 6:
            dwdp[3] = 1.0*species_4*species_5;
            break;
        case 7:
            dwdp[4] = 1.0*species_6;
            break;
        case 8:
            dwdp[5] = -1.0*species_8;
            break;
        case 9:
            dwdp[5] = 1.0*species_4*species_7;
            break;
        case 10:
            dwdp[6] = 1.0*species_8;
            break;
        case 11:
            dwdp[7] = -1.0*species_10*species_4;
            break;
        case 12:
            dwdp[7] = 1.0*species_9;
            break;
        case 13:
            dwdp[8] = -1.0*species_11;
            break;
        case 14:
            dwdp[8] = 1.0*species_10*species_5;
            break;
        case 15:
            dwdp[9] = 1.0*species_11;
            break;
        case 16:
            dwdp[10] = -1.0*species_13;
            break;
        case 17:
            dwdp[10] = 1.0*species_12*species_9;
            break;
        case 18:
            dwdp[11] = -1.0*species_15;
            break;
        case 19:
            dwdp[11] = 1.0*species_13*species_14;
            break;
        case 20:
            dwdp[12] = -1.0*species_16;
            break;
        case 21:
            dwdp[12] = 1.0*species_12*species_14;
            break;
        case 22:
            dwdp[13] = -1.0*species_15;
            break;
        case 23:
            dwdp[13] = 1.0*species_16*species_9;
            break;
        case 24:
            dwdp[14] = -1.0*species_18;
            break;
        case 25:
            dwdp[14] = 1.0*species_15*species_17;
            break;
        case 26:
            dwdp[15] = 1.0*species_18;
            break;
        case 27:
            dwdp[16] = -1.0*species_20;
            break;
        case 28:
            dwdp[16] = 1.0*species_12*species_4;
            break;
        case 29:
            dwdp[17] = -1.0*species_21;
            break;
        case 30:
            dwdp[17] = 1.0*species_14*species_20;
            break;
        case 31:
            dwdp[18] = -1.0*species_21;
            break;
        case 32:
            dwdp[18] = 1.0*species_16*species_4;
            break;
        case 33:
            dwdp[19] = -1.0*species_22;
            break;
        case 34:
            dwdp[19] = 1.0*species_17*species_21;
            break;
        case 35:
            dwdp[20] = 1.0*species_22;
            break;
        case 36:
            dwdp[21] = -1.0*species_24;
            break;
        case 37:
            dwdp[21] = 1.0*species_19*species_23;
            break;
        case 38:
            dwdp[22] = 1.0*species_24;
            break;
        case 39:
            dwdp[23] = -1.0*species_27;
            break;
        case 40:
            dwdp[23] = 1.0*species_25*species_26;
            break;
        case 41:
            dwdp[24] = 1.0*species_27;
            break;
        case 42:
            dwdp[25] = -1.0*species_29;
            break;
        case 43:
            dwdp[25] = 1.0*species_25*species_28;
            break;
        case 44:
            dwdp[26] = 1.0*species_29;
            break;
        case 45:
            dwdp[27] = -1.0*species_32;
            break;
        case 46:
            dwdp[27] = 1.0*species_30*species_31;
            break;
        case 47:
            dwdp[28] = 1.0*species_32;
            break;
        case 48:
            dwdp[29] = -1.0*species_34;
            break;
        case 49:
            dwdp[29] = 1.0*species_30*species_33;
            break;
        case 50:
            dwdp[30] = 1.0*species_34;
            break;
        case 51:
            dwdp[31] = -1.0*species_37;
            break;
        case 52:
            dwdp[31] = 1.0*species_25*species_36;
            break;
        case 53:
            dwdp[32] = 1.0*species_37;
            break;
        case 54:
            dwdp[33] = -1.0*species_39;
            break;
        case 55:
            dwdp[33] = 1.0*species_30*species_38;
            break;
        case 56:
            dwdp[34] = 1.0*species_39;
            break;
        case 57:
            dwdp[35] = -1.0*species_40;
            break;
        case 58:
            dwdp[35] = 1.0*species_28*species_38;
            break;
        case 59:
            dwdp[36] = 1.0*species_40;
            break;
        case 60:
            dwdp[37] = -1.0*species_42;
            break;
        case 61:
            dwdp[37] = 1.0*species_35*species_41;
            break;
        case 62:
            dwdp[38] = 1.0*species_42;
            break;
        case 63:
            dwdp[39] = -1.0*species_43;
            break;
        case 64:
            dwdp[39] = 1.0*species_33*species_41;
            break;
        case 65:
            dwdp[40] = 1.0*species_43;
            break;
        case 66:
            dwdp[41] = 1.0*species_19;
            break;
        case 67:
            dwdp[42] = -1.0*species_45;
            break;
        case 68:
            dwdp[42] = 1.0*species_19*species_44;
            break;
        case 69:
            dwdp[43] = 1.0*species_45;
            break;
        case 70:
            dwdp[44] = -1.0*species_46;
            break;
        case 71:
            dwdp[44] = 1.0*species_15*species_35;
            break;
        case 72:
            dwdp[45] = 1.0*species_46;
            break;
        case 73:
            dwdp[46] = -1.0*species_48;
            break;
        case 74:
            dwdp[46] = 1.0*species_21*species_35;
            break;
        case 75:
            dwdp[47] = 1.0*species_48;
            break;
        case 76:
            dwdp[48] = 1.0*species_47;
            break;
        case 77:
            dwdp[49] = -1.0*species_50;
            break;
        case 78:
            dwdp[49] = 1.0*species_4*species_49;
            break;
        case 79:
            dwdp[50] = 1.0*species_50;
            break;
        case 80:
            dwdp[51] = -1.0*species_54;
            break;
        case 81:
            dwdp[51] = 1.0*species_52*species_53;
            break;
        case 82:
            dwdp[52] = 1.0*species_54;
            break;
        case 83:
            dwdp[53] = -1.0*species_56;
            break;
        case 84:
            dwdp[53] = 1.0*species_52*species_55;
            break;
        case 85:
            dwdp[54] = 1.0*species_56;
            break;
        case 86:
            dwdp[55] = -1.0*species_59;
            break;
        case 87:
            dwdp[55] = 1.0*species_57*species_58;
            break;
        case 88:
            dwdp[56] = -1.0*species_61;
            break;
        case 89:
            dwdp[56] = 1.0*species_59*species_60;
            break;
        case 90:
            dwdp[57] = 1.0*species_61;
            break;
        case 91:
            dwdp[58] = -1.0*species_57*species_63;
            break;
        case 92:
            dwdp[58] = 1.0*species_62;
            break;
        case 93:
            dwdp[59] = -1.0*species_65;
            break;
        case 94:
            dwdp[59] = 1.0*species_62*species_64;
            break;
        case 95:
            dwdp[60] = 1.0*species_65;
            break;
        case 96:
            dwdp[61] = -1.0*species_66;
            break;
        case 97:
            dwdp[61] = 1.0*species_25*species_62;
            break;
        case 98:
            dwdp[62] = 1.0*species_66;
            break;
        case 99:
            dwdp[63] = 1.0*species_67;
            break;
        case 100:
            dwdp[64] = -1.0*species_70;
            break;
        case 101:
            dwdp[64] = 1.0*species_68*species_69;
            break;
        case 102:
            dwdp[65] = 1.0*species_70;
            break;
        case 103:
            dwdp[66] = -1.0*species_72;
            break;
        case 104:
            dwdp[66] = 1.0*species_57*species_71;
            break;
        case 105:
            dwdp[67] = 1.0*species_72;
            break;
        case 106:
            dwdp[68] = 1.0*species_71;
            break;
        case 107:
            dwdp[69] = 1.0*species_57;
            break;
        case 108:
            dwdp[70] = -1.0*species_74;
            break;
        case 109:
            dwdp[70] = 1.0*species_57*species_73;
            break;
        case 110:
            dwdp[71] = -1.0*species_76;
            break;
        case 111:
            dwdp[71] = 1.0*species_74*species_75;
            break;
        case 112:
            dwdp[72] = 1.0*species_76;
            break;
        case 113:
            dwdp[73] = -1.0*species_79;
            break;
        case 114:
            dwdp[73] = 1.0*species_75*species_78;
            break;
        case 115:
            dwdp[74] = 1.0*species_77;
            break;
        case 116:
            dwdp[75] = -1.0*species_81;
            break;
        case 117:
            dwdp[75] = 1.0*species_77*species_80;
            break;
        case 118:
            dwdp[76] = 1.0*species_81;
            break;
        case 119:
            dwdp[77] = -1.0*species_83;
            break;
        case 120:
            dwdp[77] = 1.0*species_78*species_82;
            break;
        case 121:
            dwdp[78] = -1.0*species_85;
            break;
        case 122:
            dwdp[78] = 1.0*species_82*species_84;
            break;
        case 123:
            dwdp[79] = 1.0*species_85;
            break;
        case 124:
            dwdp[80] = 1.0*species_86;
            break;
        case 125:
            dwdp[81] = -1.0*species_87;
            break;
        case 126:
            dwdp[81] = 1.0*species_4*species_44;
            break;
        case 127:
            dwdp[82] = -1.0*species_88;
            break;
        case 128:
            dwdp[82] = 1.0*species_19*species_87;
            break;
        case 129:
            dwdp[83] = 1.0*species_88;
            break;
        case 130:
            dwdp[84] = -1.0*species_91;
            break;
        case 131:
            dwdp[84] = 1.0*species_13*species_90;
            break;
        case 132:
            dwdp[85] = -1.0*species_92;
            break;
        case 133:
            dwdp[85] = 1.0*species_20*species_90;
            break;
        case 134:
            dwdp[86] = 1.0*species_91;
            break;
        case 135:
            dwdp[87] = 1.0*species_92;
            break;
        case 136:
            dwdp[88] = -1.0*species_93;
            break;
        case 137:
            dwdp[88] = 1.0*species_84*species_91;
            break;
        case 138:
            dwdp[89] = 1.0*species_93;
            break;
        case 139:
            dwdp[90] = -1.0*species_96;
            break;
        case 140:
            dwdp[90] = 1.0*species_91*species_95;
            break;
        case 141:
            dwdp[91] = 1.0*species_96;
            break;
        case 142:
            dwdp[92] = -1.0*species_98;
            break;
        case 143:
            dwdp[92] = 1.0*species_84*species_92;
            break;
        case 144:
            dwdp[93] = 1.0*species_98;
            break;
        case 145:
            dwdp[94] = -1.0*species_99;
            break;
        case 146:
            dwdp[94] = 1.0*species_92*species_95;
            break;
        case 147:
            dwdp[95] = 1.0*species_99;
            break;
        case 148:
            dwdp[96] = -1.0*species_100;
            break;
        case 149:
            dwdp[96] = 1.0*species_87*species_90;
            break;
        case 150:
            dwdp[97] = 1.0*species_100;
            break;
        case 151:
            dwdp[98] = -1.0*species_102;
            break;
        case 152:
            dwdp[98] = 1.0*species_101*species_94;
            break;
        case 153:
            dwdp[99] = 1.0*species_102;
            break;
        case 154:
            dwdp[100] = 1.0*species_84;
            break;
        case 155:
            dwdp[101] = -1.0*species_103;
            break;
        case 156:
            dwdp[101] = 1.0*species_101*species_97;
            break;
        case 157:
            dwdp[102] = 1.0*species_103;
            break;
        case 158:
            dwdp[103] = -1.0*species_104;
            break;
        case 159:
            dwdp[103] = 1.0*species_86*species_95;
            break;
        case 160:
            dwdp[104] = 1.0*species_104;
            break;
        case 161:
            dwdp[105] = -1.0*species_106;
            break;
        case 162:
            dwdp[105] = 1.0*species_105*species_86;
            break;
        case 163:
            dwdp[106] = 1.0*species_106;
            break;
        case 164:
            dwdp[107] = -1.0*species_108;
            break;
        case 165:
            dwdp[107] = 1.0*species_107*species_4;
            break;
        case 166:
            dwdp[108] = 1.0*species_108;
            break;
        case 167:
            dwdp[109] = -1.0*species_101*species_110;
            break;
        case 168:
            dwdp[109] = 1.0*species_4;
            break;
        case 169:
            dwdp[110] = -1.0*species_112;
            break;
        case 170:
            dwdp[110] = 1.0*species_101*species_111;
            break;
        case 171:
            dwdp[111] = 1.0*species_112;
            break;
        case 172:
            dwdp[112] = -1.0*species_107*species_111;
            break;
        case 173:
            dwdp[112] = 1.0*species_113;
            break;
        case 174:
            dwdp[113] = -1.0*species_115;
            break;
        case 175:
            dwdp[113] = 1.0*species_114*species_15;
            break;
        case 176:
            dwdp[114] = -1.0*species_116;
            break;
        case 177:
            dwdp[114] = 1.0*species_114*species_21;
            break;
        case 178:
            dwdp[115] = -1.0*species_118;
            break;
        case 179:
            dwdp[115] = 1.0*species_115*species_117;
            break;
        case 180:
            dwdp[116] = -1.0*species_119;
            break;
        case 181:
            dwdp[116] = 1.0*species_116*species_117;
            break;
        case 182:
            dwdp[117] = 1.0*species_118;
            break;
        case 183:
            dwdp[118] = 1.0*species_119;
            break;
        case 184:
            dwdp[119] = 1.0*species_122;
            break;
        case 185:
            dwdp[120] = -1.0*species_123;
            break;
        case 186:
            dwdp[120] = 1.0*species_117*species_68;
            break;
        case 187:
            dwdp[121] = 1.0*species_123;
            break;
        case 188:
            dwdp[122] = -1.0*species_126;
            break;
        case 189:
            dwdp[122] = 1.0*species_124*species_125;
            break;
        case 190:
            dwdp[123] = 1.0*species_127;
            break;
        case 191:
            dwdp[124] = -1.0*species_128;
            break;
        case 192:
            dwdp[124] = 1.0*species_19*species_94;
            break;
        case 193:
            dwdp[125] = 1.0*species_128;
            break;
        case 194:
            dwdp[126] = -1.0*species_129;
            break;
        case 195:
            dwdp[126] = 1.0*species_35*species_68;
            break;
        case 196:
            dwdp[127] = 1.0*species_129;
            break;
        case 197:
            dwdp[128] = -1.0*species_131;
            break;
        case 198:
            dwdp[128] = 1.0*species_12*species_130;
            break;
        case 199:
            dwdp[129] = -1.0*species_132;
            break;
        case 200:
            dwdp[129] = 1.0*species_131*species_9;
            break;
        case 201:
            dwdp[130] = -1.0*species_133;
            break;
        case 202:
            dwdp[130] = 1.0*species_131*species_4;
            break;
        case 203:
            dwdp[131] = -1.0*species_134;
            break;
        case 204:
            dwdp[131] = 1.0*species_132*species_26;
            break;
        case 205:
            dwdp[132] = 1.0*species_134;
            break;
        case 206:
            dwdp[133] = 1.0*species_135;
            break;
        case 207:
            dwdp[134] = -1.0*species_132*species_30;
            break;
        case 208:
            dwdp[134] = 1.0*species_136;
            break;
        case 209:
            dwdp[135] = -1.0*species_137;
            break;
        case 210:
            dwdp[135] = 1.0*species_133*species_26;
            break;
        case 211:
            dwdp[136] = 1.0*species_137;
            break;
        case 212:
            dwdp[137] = 1.0*species_138;
            break;
        case 213:
            dwdp[138] = -1.0*species_133*species_30;
            break;
        case 214:
            dwdp[138] = 1.0*species_139;
            break;
        case 215:
            dwdp[139] = -1.0*species_140;
            break;
        case 216:
            dwdp[139] = 1.0*species_132*species_23;
            break;
        case 217:
            dwdp[140] = -1.0*species_141;
            break;
        case 218:
            dwdp[140] = 1.0*species_140*species_26;
            break;
        case 219:
            dwdp[141] = 1.0*species_141;
            break;
        case 220:
            dwdp[142] = 1.0*species_142;
            break;
        case 221:
            dwdp[143] = -1.0*species_140*species_30;
            break;
        case 222:
            dwdp[143] = 1.0*species_143;
            break;
        case 223:
            dwdp[144] = -1.0*species_144;
            break;
        case 224:
            dwdp[144] = 1.0*species_133*species_23;
            break;
        case 225:
            dwdp[145] = -1.0*species_145;
            break;
        case 226:
            dwdp[145] = 1.0*species_144*species_26;
            break;
        case 227:
            dwdp[146] = 1.0*species_145;
            break;
        case 228:
            dwdp[147] = 1.0*species_146;
            break;
        case 229:
            dwdp[148] = -1.0*species_144*species_30;
            break;
        case 230:
            dwdp[148] = 1.0*species_147;
            break;
        case 231:
            dwdp[149] = -1.0*species_148;
            break;
        case 232:
            dwdp[149] = 1.0*species_132*species_19;
            break;
        case 233:
            dwdp[150] = -1.0*species_149;
            break;
        case 234:
            dwdp[150] = 1.0*species_148*species_26;
            break;
        case 235:
            dwdp[151] = 1.0*species_149;
            break;
        case 236:
            dwdp[152] = 1.0*species_150;
            break;
        case 237:
            dwdp[153] = -1.0*species_148*species_30;
            break;
        case 238:
            dwdp[153] = 1.0*species_151;
            break;
        case 239:
            dwdp[154] = -1.0*species_152;
            break;
        case 240:
            dwdp[154] = 1.0*species_133*species_19;
            break;
        case 241:
            dwdp[155] = -1.0*species_153;
            break;
        case 242:
            dwdp[155] = 1.0*species_152*species_26;
            break;
        case 243:
            dwdp[156] = 1.0*species_153;
            break;
        case 244:
            dwdp[157] = 1.0*species_154;
            break;
        case 245:
            dwdp[158] = -1.0*species_152*species_30;
            break;
        case 246:
            dwdp[158] = 1.0*species_155;
            break;
        case 247:
            dwdp[159] = -1.0*species_156;
            break;
        case 248:
            dwdp[159] = 1.0*species_151*species_31;
            break;
        case 249:
            dwdp[160] = 1.0*species_156;
            break;
        case 250:
            dwdp[161] = 1.0*species_157;
            break;
        case 251:
            dwdp[162] = 1.0*species_158;
            break;
        case 252:
            dwdp[163] = -1.0*species_159;
            break;
        case 253:
            dwdp[163] = 1.0*species_155*species_31;
            break;
        case 254:
            dwdp[164] = 1.0*species_159;
            break;
        case 255:
            dwdp[165] = 1.0*species_160;
            break;
        case 256:
            dwdp[166] = 1.0*species_161;
            break;
        case 257:
            dwdp[167] = -1.0*species_162;
            break;
        case 258:
            dwdp[167] = 1.0*species_132*species_25;
            break;
        case 259:
            dwdp[168] = -1.0*species_163;
            break;
        case 260:
            dwdp[168] = 1.0*species_162*species_26;
            break;
        case 261:
            dwdp[169] = 1.0*species_163;
            break;
        case 262:
            dwdp[170] = 1.0*species_164;
            break;
        case 263:
            dwdp[171] = -1.0*species_162*species_30;
            break;
        case 264:
            dwdp[171] = 1.0*species_165;
            break;
        case 265:
            dwdp[172] = -1.0*species_166;
            break;
        case 266:
            dwdp[172] = 1.0*species_133*species_25;
            break;
        case 267:
            dwdp[173] = -1.0*species_167;
            break;
        case 268:
            dwdp[173] = 1.0*species_166*species_26;
            break;
        case 269:
            dwdp[174] = 1.0*species_167;
            break;
        case 270:
            dwdp[175] = 1.0*species_168;
            break;
        case 271:
            dwdp[176] = -1.0*species_166*species_30;
            break;
        case 272:
            dwdp[176] = 1.0*species_169;
            break;
        case 273:
            dwdp[177] = -1.0*species_170;
            break;
        case 274:
            dwdp[177] = 1.0*species_134*species_86;
            break;
        case 275:
            dwdp[178] = -1.0*species_172;
            break;
        case 276:
            dwdp[178] = 1.0*species_137*species_86;
            break;
        case 277:
            dwdp[179] = -1.0*species_174;
            break;
        case 278:
            dwdp[179] = 1.0*species_141*species_86;
            break;
        case 279:
            dwdp[180] = -1.0*species_176;
            break;
        case 280:
            dwdp[180] = 1.0*species_145*species_86;
            break;
        case 281:
            dwdp[181] = -1.0*species_178;
            break;
        case 282:
            dwdp[181] = 1.0*species_136*species_31;
            break;
        case 283:
            dwdp[182] = 1.0*species_178;
            break;
        case 284:
            dwdp[183] = 1.0*species_179;
            break;
        case 285:
            dwdp[184] = 1.0*species_180;
            break;
        case 286:
            dwdp[185] = -1.0*species_181;
            break;
        case 287:
            dwdp[185] = 1.0*species_139*species_31;
            break;
        case 288:
            dwdp[186] = 1.0*species_181;
            break;
        case 289:
            dwdp[187] = 1.0*species_182;
            break;
        case 290:
            dwdp[188] = 1.0*species_183;
            break;
        case 291:
            dwdp[189] = -1.0*species_184;
            break;
        case 292:
            dwdp[189] = 1.0*species_143*species_31;
            break;
        case 293:
            dwdp[190] = 1.0*species_184;
            break;
        case 294:
            dwdp[191] = 1.0*species_185;
            break;
        case 295:
            dwdp[192] = 1.0*species_186;
            break;
        case 296:
            dwdp[193] = -1.0*species_187;
            break;
        case 297:
            dwdp[193] = 1.0*species_147*species_31;
            break;
        case 298:
            dwdp[194] = 1.0*species_187;
            break;
        case 299:
            dwdp[195] = 1.0*species_188;
            break;
        case 300:
            dwdp[196] = 1.0*species_189;
            break;
        case 301:
            dwdp[197] = -1.0*species_190;
            break;
        case 302:
            dwdp[197] = 1.0*species_171*species_95;
            break;
        case 303:
            dwdp[198] = 1.0*species_190;
            break;
        case 304:
            dwdp[199] = -1.0*species_191;
            break;
        case 305:
            dwdp[199] = 1.0*species_173*species_95;
            break;
        case 306:
            dwdp[200] = 1.0*species_191;
            break;
        case 307:
            dwdp[201] = -1.0*species_192;
            break;
        case 308:
            dwdp[201] = 1.0*species_175*species_95;
            break;
        case 309:
            dwdp[202] = 1.0*species_192;
            break;
        case 310:
            dwdp[203] = -1.0*species_193;
            break;
        case 311:
            dwdp[203] = 1.0*species_177*species_95;
            break;
        case 312:
            dwdp[204] = 1.0*species_193;
            break;
    }
}