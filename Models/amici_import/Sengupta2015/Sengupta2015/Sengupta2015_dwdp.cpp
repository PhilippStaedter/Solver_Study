#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Sengupta2015(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = s11;
            break;
        case 1:
            dwdp[1] = s340;
            break;
        case 2:
            dwdp[2] = s341;
            break;
        case 3:
            dwdp[3] = s343;
            break;
        case 4:
            dwdp[4] = s73;
            break;
        case 5:
            dwdp[5] = s342*s345;
            break;
        case 6:
            dwdp[6] = s23*s345*s38;
            break;
        case 7:
            dwdp[7] = s293;
            break;
        case 8:
            dwdp[8] = s296*s67;
            break;
        case 9:
            dwdp[9] = s3*s63;
            break;
        case 10:
            dwdp[10] = s63*s8;
            break;
        case 11:
            dwdp[11] = s35*s44*s93;
            break;
        case 12:
            dwdp[12] = s73;
            break;
        case 13:
            dwdp[13] = s335;
            break;
        case 14:
            dwdp[14] = s262;
            break;
        case 15:
            dwdp[15] = s256;
            break;
        case 16:
            dwdp[16] = s347;
            break;
        case 17:
            dwdp[17] = s358;
            break;
        case 18:
            dwdp[18] = s374;
            break;
        case 19:
            dwdp[19] = s358;
            break;
        case 20:
            dwdp[20] = s263;
            break;
        case 21:
            dwdp[21] = s388;
            break;
        case 22:
            dwdp[22] = s388;
            break;
        case 23:
            dwdp[23] = s358;
            break;
        case 24:
            dwdp[24] = s353*s377;
            break;
        case 25:
            dwdp[25] = s381;
            break;
        case 26:
            dwdp[26] = s51*s53;
            break;
        case 27:
            dwdp[27] = s11;
            break;
        case 28:
            dwdp[28] = s307;
            break;
        case 29:
            dwdp[29] = s322;
            break;
        case 30:
            dwdp[30] = s378;
            break;
        case 32:
            dwdp[31] = s342;
            break;
        case 33:
            dwdp[32] = s33*s34;
            break;
        case 34:
            dwdp[33] = s39;
            break;
        case 35:
            dwdp[34] = s343*s347;
            break;
        case 36:
            dwdp[35] = s37;
            break;
        case 37:
            dwdp[36] = s66;
            break;
        case 38:
            dwdp[37] = s116;
            break;
        case 39:
            dwdp[38] = s80;
            break;
        case 40:
            dwdp[39] = s111;
            break;
        case 41:
            dwdp[40] = s108;
            break;
        case 42:
            dwdp[41] = s89;
            break;
        case 43:
            dwdp[42] = s100;
            break;
        case 44:
            dwdp[43] = s103;
            break;
        case 45:
            dwdp[44] = s353;
            break;
        case 46:
            dwdp[45] = s269;
            break;
        case 47:
            dwdp[46] = s355;
            break;
        case 48:
            dwdp[47] = s355;
            break;
        case 49:
            dwdp[48] = s355;
            break;
        case 50:
            dwdp[49] = s355;
            break;
        case 51:
            dwdp[50] = s356;
            break;
        case 52:
            dwdp[51] = s265;
            break;
        case 53:
            dwdp[52] = s269;
            break;
        case 54:
            dwdp[53] = s377;
            break;
        case 55:
            dwdp[54] = s14*s21;
            break;
        case 56:
            dwdp[55] = s12*s17;
            break;
        case 57:
            dwdp[56] = s351;
            break;
        case 58:
            dwdp[57] = s351;
            break;
        case 59:
            dwdp[58] = s351;
            break;
        case 60:
            dwdp[59] = s351;
            break;
        case 61:
            dwdp[60] = s267;
            break;
        case 62:
            dwdp[61] = s259;
            break;
        case 63:
            dwdp[62] = s269;
            break;
        case 64:
            dwdp[63] = s388;
            break;
        case 65:
            dwdp[64] = s388;
            break;
        case 66:
            dwdp[65] = -re89_v1*s302/pow(re89_k1 + s302, 2);
            break;
        case 67:
            dwdp[65] = s302/(re89_k1 + s302);
            break;
        case 68:
            dwdp[66] = -re25_v1*s71/pow(re25_k1 + s71, 2);
            break;
        case 69:
            dwdp[66] = s71/(re25_k1 + s71);
            break;
        case 70:
            dwdp[67] = -re4_v1*s340/pow(re4_k1 + s340, 2);
            break;
        case 71:
            dwdp[67] = s340/(re4_k1 + s340);
            break;
        case 72:
            dwdp[68] = s379*s50;
            break;
        case 73:
            dwdp[69] = -re20_v1*s7/pow(re20_k1 + s7, 2);
            break;
        case 74:
            dwdp[69] = s7/(re20_k1 + s7);
            break;
        case 75:
            dwdp[70] = -re21_v1*s40/pow(re21_k1 + s40, 2);
            break;
        case 76:
            dwdp[70] = s40/(re21_k1 + s40);
            break;
        case 77:
            dwdp[71] = -re22_v1*s9/pow(re22_k1 + s9, 2);
            break;
        case 78:
            dwdp[71] = s9/(re22_k1 + s9);
            break;
        case 79:
            dwdp[72] = -re33_v1*s73/pow(re33_k1 + s73, 2);
            break;
        case 80:
            dwdp[72] = s73/(re33_k1 + s73);
            break;
        case 81:
            dwdp[73] = -re39_v1*s35/pow(re39_k1 + s35, 2);
            break;
        case 82:
            dwdp[73] = s35/(re39_k1 + s35);
            break;
        case 83:
            dwdp[74] = -re6_v1*s25/pow(re6_k1 + s25, 2);
            break;
        case 84:
            dwdp[74] = s25/(re6_k1 + s25);
            break;
        case 85:
            dwdp[75] = -re29_v1*s8/pow(re29_k1 + s8, 2);
            break;
        case 86:
            dwdp[75] = s8/(re29_k1 + s8);
            break;
        case 87:
            dwdp[76] = -re31_v1*s10/pow(re31_k1 + s10, 2);
            break;
        case 88:
            dwdp[76] = s10/(re31_k1 + s10);
            break;
        case 89:
            dwdp[77] = -re15_v1*s348/pow(re15_k1 + s348, 2);
            break;
        case 90:
            dwdp[77] = s348/(re15_k1 + s348);
            break;
        case 91:
            dwdp[78] = -re23_v1*s52/pow(re23_k1 + s52, 2);
            break;
        case 92:
            dwdp[78] = s52/(re23_k1 + s52);
            break;
        case 93:
            dwdp[79] = -re18_v1*s5/pow(re18_k1 + s5, 2);
            break;
        case 94:
            dwdp[79] = s5/(re18_k1 + s5);
            break;
        case 95:
            dwdp[80] = -re19_v1*s6/pow(re19_k1 + s6, 2);
            break;
        case 96:
            dwdp[80] = s6/(re19_k1 + s6);
            break;
        case 97:
            dwdp[81] = -re32_v1*s234/pow(re32_k1 + s234, 2);
            break;
        case 98:
            dwdp[81] = s234/(re32_k1 + s234);
            break;
        case 99:
            dwdp[82] = -re26_v1*s234/pow(re26_k1 + s234, 2);
            break;
        case 100:
            dwdp[82] = s234/(re26_k1 + s234);
            break;
        case 101:
            dwdp[83] = -re41_v1*s74/pow(re41_k1 + s74, 2);
            break;
        case 102:
            dwdp[83] = s74/(re41_k1 + s74);
            break;
        case 103:
            dwdp[84] = -re28_v1*s74/pow(re28_k1 + s74, 2);
            break;
        case 104:
            dwdp[84] = s74/(re28_k1 + s74);
            break;
        case 105:
            dwdp[85] = -re30_v1*s75/pow(re30_k1 + s75, 2);
            break;
        case 106:
            dwdp[85] = s75/(re30_k1 + s75);
            break;
        case 107:
            dwdp[86] = -re35_v1*s81/pow(re35_k1 + s81, 2);
            break;
        case 108:
            dwdp[86] = s81/(re35_k1 + s81);
            break;
        case 109:
            dwdp[87] = -re34_v1*s11/pow(re34_k1 + s11, 2);
            break;
        case 110:
            dwdp[87] = s11/(re34_k1 + s11);
            break;
        case 111:
            dwdp[88] = -re42_v1*s11/pow(re42_k1 + s11, 2);
            break;
        case 112:
            dwdp[88] = s11/(re42_k1 + s11);
            break;
        case 113:
            dwdp[89] = s127;
            break;
        case 114:
            dwdp[90] = s135;
            break;
        case 115:
            dwdp[91] = s133;
            break;
        case 116:
            dwdp[92] = s131;
            break;
        case 117:
            dwdp[93] = s129;
            break;
        case 118:
            dwdp[94] = s125*s126;
            break;
        case 119:
            dwdp[95] = s126*s128;
            break;
        case 120:
            dwdp[96] = s126*s130;
            break;
        case 121:
            dwdp[97] = s126*s132;
            break;
        case 122:
            dwdp[98] = s126*s134;
            break;
        case 123:
            dwdp[99] = s367;
            break;
        case 124:
            dwdp[100] = s371;
            break;
        case 125:
            dwdp[101] = s370;
            break;
        case 126:
            dwdp[102] = s369;
            break;
        case 127:
            dwdp[103] = s368;
            break;
        case 128:
            dwdp[104] = -re43_v1*s234/pow(re43_k1 + s234, 2);
            break;
        case 129:
            dwdp[104] = s234/(re43_k1 + s234);
            break;
        case 130:
            dwdp[105] = -re44_v1*s181/pow(re44_k1 + s181, 2);
            break;
        case 131:
            dwdp[105] = s181/(re44_k1 + s181);
            break;
        case 132:
            dwdp[106] = -re45_v1*s182/pow(re45_k1 + s182, 2);
            break;
        case 133:
            dwdp[106] = s182/(re45_k1 + s182);
            break;
        case 134:
            dwdp[107] = -re46_v1*s183/pow(re46_k1 + s183, 2);
            break;
        case 135:
            dwdp[107] = s183/(re46_k1 + s183);
            break;
        case 136:
            dwdp[108] = -re47_v1*s187/pow(re47_k1 + s187, 2);
            break;
        case 137:
            dwdp[108] = s187/(re47_k1 + s187);
            break;
        case 138:
            dwdp[109] = -re48_v1*s195/pow(re48_k1 + s195, 2);
            break;
        case 139:
            dwdp[109] = s195/(re48_k1 + s195);
            break;
        case 140:
            dwdp[110] = -re53_v1*s35/pow(re53_k1 + s35, 2);
            break;
        case 141:
            dwdp[110] = s35/(re53_k1 + s35);
            break;
        case 142:
            dwdp[111] = -re54_v1*s234/pow(re54_k1 + s234, 2);
            break;
        case 143:
            dwdp[111] = s234/(re54_k1 + s234);
            break;
        case 144:
            dwdp[112] = -re55_v1*s198/pow(re55_k1 + s198, 2);
            break;
        case 145:
            dwdp[112] = s198/(re55_k1 + s198);
            break;
        case 146:
            dwdp[113] = -re56_v1*s252/pow(re56_k1 + s252, 2);
            break;
        case 147:
            dwdp[113] = s252/(re56_k1 + s252);
            break;
        case 148:
            dwdp[114] = -re57_v1*s253/pow(re57_k1 + s253, 2);
            break;
        case 149:
            dwdp[114] = s253/(re57_k1 + s253);
            break;
        case 150:
            dwdp[115] = -re80_v1*s284/pow(re80_k1 + s284, 2);
            break;
        case 151:
            dwdp[115] = s284/(re80_k1 + s284);
            break;
        case 152:
            dwdp[116] = -re79_v1*s286/pow(re79_k1 + s286, 2);
            break;
        case 153:
            dwdp[116] = s286/(re79_k1 + s286);
            break;
        case 154:
            dwdp[117] = -re83_v1*s285/pow(re83_k1 + s285, 2);
            break;
        case 155:
            dwdp[117] = s285/(re83_k1 + s285);
            break;
        case 156:
            dwdp[118] = -re84_v1*s293/pow(re84_k1 + s293, 2);
            break;
        case 157:
            dwdp[118] = s293/(re84_k1 + s293);
            break;
        case 158:
            dwdp[119] = -re87_v1*s297/pow(re87_k1 + s297, 2);
            break;
        case 159:
            dwdp[119] = s297/(re87_k1 + s297);
            break;
        case 160:
            dwdp[120] = -re88_v1*s306/pow(re88_k1 + s306, 2);
            break;
        case 161:
            dwdp[120] = s306/(re88_k1 + s306);
            break;
        case 162:
            dwdp[121] = -re91_v1*s335/pow(re91_k1 + s335, 2);
            break;
        case 163:
            dwdp[121] = s335/(re91_k1 + s335);
            break;
        case 164:
            dwdp[122] = -re92_v1*s321/pow(re92_k1 + s321, 2);
            break;
        case 165:
            dwdp[122] = s321/(re92_k1 + s321);
            break;
        case 166:
            dwdp[123] = -re94_v1*s323/pow(re94_k1 + s323, 2);
            break;
        case 167:
            dwdp[123] = s323/(re94_k1 + s323);
            break;
        case 168:
            dwdp[124] = -re95_v1*s323/pow(re95_k1 + s323, 2);
            break;
        case 169:
            dwdp[124] = s323/(re95_k1 + s323);
            break;
        case 170:
            dwdp[125] = -re98_v1*s323/pow(re98_k1 + s323, 2);
            break;
        case 171:
            dwdp[125] = s323/(re98_k1 + s323);
            break;
        case 172:
            dwdp[126] = -re50_v1*s197/pow(re50_k1 + s197, 2);
            break;
        case 173:
            dwdp[126] = s197/(re50_k1 + s197);
            break;
        case 174:
            dwdp[127] = -re52_v1*s188/pow(re52_k1 + s188, 2);
            break;
        case 175:
            dwdp[127] = s188/(re52_k1 + s188);
            break;
        case 176:
            dwdp[128] = -re51_v1*s188/pow(re51_k1 + s188, 2);
            break;
        case 177:
            dwdp[128] = s188/(re51_k1 + s188);
            break;
        case 178:
            dwdp[129] = -re16_v1*s2/pow(re16_k1 + s2, 2);
            break;
        case 179:
            dwdp[129] = s2/(re16_k1 + s2);
            break;
        case 180:
            dwdp[130] = -re17_v1*s4/pow(re17_k1 + s4, 2);
            break;
        case 181:
            dwdp[130] = s4/(re17_k1 + s4);
            break;
        case 182:
            dwdp[131] = -re113_v1*s37/pow(re113_k1 + s37, 2);
            break;
        case 183:
            dwdp[131] = s37/(re113_k1 + s37);
            break;
        case 184:
            dwdp[132] = -re116_v1*s66/pow(re116_k1 + s66, 2);
            break;
        case 185:
            dwdp[132] = s66/(re116_k1 + s66);
            break;
        case 186:
            dwdp[133] = -re120_v1*s80/pow(re120_k1 + s80, 2);
            break;
        case 187:
            dwdp[133] = s80/(re120_k1 + s80);
            break;
        case 188:
            dwdp[134] = -re123_v1*s89/pow(re123_k1 + s89, 2);
            break;
        case 189:
            dwdp[134] = s89/(re123_k1 + s89);
            break;
        case 190:
            dwdp[135] = -re128_v1*s100/pow(re128_k1 + s100, 2);
            break;
        case 191:
            dwdp[135] = s100/(re128_k1 + s100);
            break;
        case 192:
            dwdp[136] = -re131_v1*s103/pow(re131_k1 + s103, 2);
            break;
        case 193:
            dwdp[136] = s103/(re131_k1 + s103);
            break;
        case 194:
            dwdp[137] = -re135_v1*s108/pow(re135_k1 + s108, 2);
            break;
        case 195:
            dwdp[137] = s108/(re135_k1 + s108);
            break;
        case 196:
            dwdp[138] = -re138_v2*s111/pow(re138_k2 + s111, 2);
            break;
        case 197:
            dwdp[138] = s111/(re138_k2 + s111);
            break;
        case 200:
            dwdp[139] = -re144_v1*s116/pow(re144_k1 + s116, 2);
            break;
        case 201:
            dwdp[139] = s116/(re144_k1 + s116);
            break;
        case 202:
            dwdp[140] = -re97_v1*s348/pow(re97_k1 + s348, 2);
            break;
        case 203:
            dwdp[140] = s348/(re97_k1 + s348);
            break;
        case 204:
            dwdp[141] = -re96_v1*s321/pow(re96_k1 + s321, 2);
            break;
        case 205:
            dwdp[141] = s321/(re96_k1 + s321);
            break;
        case 206:
            dwdp[142] = -re147_v1*s119/pow(re147_k1 + s119, 2);
            break;
        case 207:
            dwdp[142] = s119/(re147_k1 + s119);
            break;
        case 208:
            dwdp[143] = -re110_v1*s15/pow(re110_k1 + s15, 2);
            break;
        case 209:
            dwdp[143] = s15/(re110_k1 + s15);
            break;
        case 210:
            dwdp[144] = -re117_v1*s42/pow(re117_k1 + s42, 2);
            break;
        case 211:
            dwdp[144] = s42/(re117_k1 + s42);
            break;
        case 212:
            dwdp[145] = -re125_v1*s65/pow(re125_k1 + s65, 2);
            break;
        case 213:
            dwdp[145] = s65/(re125_k1 + s65);
            break;
        case 214:
            dwdp[146] = -re124_v1*s83/pow(re124_k1 + s83, 2);
            break;
        case 215:
            dwdp[146] = s83/(re124_k1 + s83);
            break;
        case 216:
            dwdp[147] = -re141_v1*s84/pow(re141_k1 + s84, 2);
            break;
        case 217:
            dwdp[147] = s84/(re141_k1 + s84);
            break;
        case 218:
            dwdp[148] = -re132_v1*s101/pow(re132_k1 + s101, 2);
            break;
        case 219:
            dwdp[148] = s101/(re132_k1 + s101);
            break;
        case 220:
            dwdp[149] = -re140_v1*s102/pow(re140_k1 + s102, 2);
            break;
        case 221:
            dwdp[149] = s102/(re140_k1 + s102);
            break;
        case 222:
            dwdp[150] = -re139_v1*s109/pow(re139_k1 + s109, 2);
            break;
        case 223:
            dwdp[150] = s109/(re139_k1 + s109);
            break;
        case 224:
            dwdp[151] = -re149_v1*s110/pow(re149_k1 + s110, 2);
            break;
        case 225:
            dwdp[151] = s110/(re149_k1 + s110);
            break;
        case 226:
            dwdp[152] = -re148_v1*s117/pow(re148_k1 + s117, 2);
            break;
        case 227:
            dwdp[152] = s117/(re148_k1 + s117);
            break;
        case 228:
            dwdp[153] = -re111_v1*s19/pow(re111_k1 + s19, 2);
            break;
        case 229:
            dwdp[153] = s19/(re111_k1 + s19);
            break;
        case 230:
            dwdp[154] = -re114_v1*s69/pow(re114_k1 + s69, 2);
            break;
        case 231:
            dwdp[154] = s69/(re114_k1 + s69);
            break;
        case 232:
            dwdp[155] = -re118_v1*s70/pow(re118_k1 + s70, 2);
            break;
        case 233:
            dwdp[155] = s70/(re118_k1 + s70);
            break;
        case 234:
            dwdp[156] = -re121_v1*s97/pow(re121_k1 + s97, 2);
            break;
        case 235:
            dwdp[156] = s97/(re121_k1 + s97);
            break;
        case 236:
            dwdp[157] = -re126_v1*s98/pow(re126_k1 + s98, 2);
            break;
        case 237:
            dwdp[157] = s98/(re126_k1 + s98);
            break;
        case 238:
            dwdp[158] = -re129_v1*s105/pow(re129_k1 + s105, 2);
            break;
        case 239:
            dwdp[158] = s105/(re129_k1 + s105);
            break;
        case 240:
            dwdp[159] = -re133_v1*s106/pow(re133_k1 + s106, 2);
            break;
        case 241:
            dwdp[159] = s106/(re133_k1 + s106);
            break;
        case 242:
            dwdp[160] = -re136_v1*s113/pow(re136_k1 + s113, 2);
            break;
        case 243:
            dwdp[160] = s113/(re136_k1 + s113);
            break;
        case 244:
            dwdp[161] = -re142_v1*s114/pow(re142_k1 + s114, 2);
            break;
        case 245:
            dwdp[161] = s114/(re142_k1 + s114);
            break;
        case 246:
            dwdp[162] = -re145_v1*s121/pow(re145_k1 + s121, 2);
            break;
        case 247:
            dwdp[162] = s121/(re145_k1 + s121);
            break;
        case 248:
            dwdp[163] = -re112_v1*s22/pow(re112_k1 + s22, 2);
            break;
        case 249:
            dwdp[163] = s22/(re112_k1 + s22);
            break;
        case 250:
            dwdp[164] = -re115_v1*s68/pow(re115_k1 + s68, 2);
            break;
        case 251:
            dwdp[164] = s68/(re115_k1 + s68);
            break;
        case 252:
            dwdp[165] = -re119_v1*s72/pow(re119_k1 + s72, 2);
            break;
        case 253:
            dwdp[165] = s72/(re119_k1 + s72);
            break;
        case 254:
            dwdp[166] = -re122_v1*s96/pow(re122_k1 + s96, 2);
            break;
        case 255:
            dwdp[166] = s96/(re122_k1 + s96);
            break;
        case 256:
            dwdp[167] = -re127_v1*s99/pow(re127_k1 + s99, 2);
            break;
        case 257:
            dwdp[167] = s99/(re127_k1 + s99);
            break;
        case 258:
            dwdp[168] = -re130_v1*s104/pow(re130_k1 + s104, 2);
            break;
        case 259:
            dwdp[168] = s104/(re130_k1 + s104);
            break;
        case 260:
            dwdp[169] = -re134_v1*s107/pow(re134_k1 + s107, 2);
            break;
        case 261:
            dwdp[169] = s107/(re134_k1 + s107);
            break;
        case 262:
            dwdp[170] = -re137_v1*s112/pow(re137_k1 + s112, 2);
            break;
        case 263:
            dwdp[170] = s112/(re137_k1 + s112);
            break;
        case 264:
            dwdp[171] = -re143_v1*s115/pow(re143_k1 + s115, 2);
            break;
        case 265:
            dwdp[171] = s115/(re143_k1 + s115);
            break;
        case 266:
            dwdp[172] = -re146_v1*s120/pow(re146_k1 + s120, 2);
            break;
        case 267:
            dwdp[172] = s120/(re146_k1 + s120);
            break;
        case 268:
            dwdp[173] = -re175_v1*s124/pow(re175_k1 + s124, 2);
            break;
        case 269:
            dwdp[173] = s124/(re175_k1 + s124);
            break;
        case 270:
            dwdp[174] = -re182_v1*s123/pow(re182_k1 + s123, 2);
            break;
        case 271:
            dwdp[174] = s123/(re182_k1 + s123);
            break;
    }
}