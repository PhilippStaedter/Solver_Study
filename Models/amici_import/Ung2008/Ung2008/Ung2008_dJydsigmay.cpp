#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydsigmay_Ung2008(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigmay0 - 1.0*pow(-my0 + y0, 2)/pow(sigmay0, 3);
            break;
        case 1:
            dJydsigmay[1] = 1.0/sigmay1 - 1.0*pow(-my1 + y1, 2)/pow(sigmay1, 3);
            break;
        case 2:
            dJydsigmay[2] = 1.0/sigmay2 - 1.0*pow(-my2 + y2, 2)/pow(sigmay2, 3);
            break;
        case 3:
            dJydsigmay[3] = 1.0/sigmay3 - 1.0*pow(-my3 + y3, 2)/pow(sigmay3, 3);
            break;
        case 4:
            dJydsigmay[4] = 1.0/sigmay4 - 1.0*pow(-my4 + y4, 2)/pow(sigmay4, 3);
            break;
        case 5:
            dJydsigmay[5] = 1.0/sigmay5 - 1.0*pow(-my5 + y5, 2)/pow(sigmay5, 3);
            break;
        case 6:
            dJydsigmay[6] = 1.0/sigmay6 - 1.0*pow(-my6 + y6, 2)/pow(sigmay6, 3);
            break;
        case 7:
            dJydsigmay[7] = 1.0/sigmay7 - 1.0*pow(-my7 + y7, 2)/pow(sigmay7, 3);
            break;
        case 8:
            dJydsigmay[8] = 1.0/sigmay8 - 1.0*pow(-my8 + y8, 2)/pow(sigmay8, 3);
            break;
        case 9:
            dJydsigmay[9] = 1.0/sigmay9 - 1.0*pow(-my9 + y9, 2)/pow(sigmay9, 3);
            break;
        case 10:
            dJydsigmay[10] = 1.0/sigmay10 - 1.0*pow(-my10 + y10, 2)/pow(sigmay10, 3);
            break;
        case 11:
            dJydsigmay[11] = 1.0/sigmay11 - 1.0*pow(-my11 + y11, 2)/pow(sigmay11, 3);
            break;
        case 12:
            dJydsigmay[12] = 1.0/sigmay12 - 1.0*pow(-my12 + y12, 2)/pow(sigmay12, 3);
            break;
        case 13:
            dJydsigmay[13] = 1.0/sigmay13 - 1.0*pow(-my13 + y13, 2)/pow(sigmay13, 3);
            break;
        case 14:
            dJydsigmay[14] = 1.0/sigmay14 - 1.0*pow(-my14 + y14, 2)/pow(sigmay14, 3);
            break;
        case 15:
            dJydsigmay[15] = 1.0/sigmay15 - 1.0*pow(-my15 + y15, 2)/pow(sigmay15, 3);
            break;
        case 16:
            dJydsigmay[16] = 1.0/sigmay16 - 1.0*pow(-my16 + y16, 2)/pow(sigmay16, 3);
            break;
        case 17:
            dJydsigmay[17] = 1.0/sigmay17 - 1.0*pow(-my17 + y17, 2)/pow(sigmay17, 3);
            break;
        case 18:
            dJydsigmay[18] = 1.0/sigmay18 - 1.0*pow(-my18 + y18, 2)/pow(sigmay18, 3);
            break;
        case 19:
            dJydsigmay[19] = 1.0/sigmay19 - 1.0*pow(-my19 + y19, 2)/pow(sigmay19, 3);
            break;
        case 20:
            dJydsigmay[20] = 1.0/sigmay20 - 1.0*pow(-my20 + y20, 2)/pow(sigmay20, 3);
            break;
        case 21:
            dJydsigmay[21] = 1.0/sigmay21 - 1.0*pow(-my21 + y21, 2)/pow(sigmay21, 3);
            break;
        case 22:
            dJydsigmay[22] = 1.0/sigmay22 - 1.0*pow(-my22 + y22, 2)/pow(sigmay22, 3);
            break;
        case 23:
            dJydsigmay[23] = 1.0/sigmay23 - 1.0*pow(-my23 + y23, 2)/pow(sigmay23, 3);
            break;
        case 24:
            dJydsigmay[24] = 1.0/sigmay24 - 1.0*pow(-my24 + y24, 2)/pow(sigmay24, 3);
            break;
        case 25:
            dJydsigmay[25] = 1.0/sigmay25 - 1.0*pow(-my25 + y25, 2)/pow(sigmay25, 3);
            break;
        case 26:
            dJydsigmay[26] = 1.0/sigmay26 - 1.0*pow(-my26 + y26, 2)/pow(sigmay26, 3);
            break;
        case 27:
            dJydsigmay[27] = 1.0/sigmay27 - 1.0*pow(-my27 + y27, 2)/pow(sigmay27, 3);
            break;
        case 28:
            dJydsigmay[28] = 1.0/sigmay28 - 1.0*pow(-my28 + y28, 2)/pow(sigmay28, 3);
            break;
        case 29:
            dJydsigmay[29] = 1.0/sigmay29 - 1.0*pow(-my29 + y29, 2)/pow(sigmay29, 3);
            break;
        case 30:
            dJydsigmay[30] = 1.0/sigmay30 - 1.0*pow(-my30 + y30, 2)/pow(sigmay30, 3);
            break;
        case 31:
            dJydsigmay[31] = 1.0/sigmay31 - 1.0*pow(-my31 + y31, 2)/pow(sigmay31, 3);
            break;
        case 32:
            dJydsigmay[32] = 1.0/sigmay32 - 1.0*pow(-my32 + y32, 2)/pow(sigmay32, 3);
            break;
        case 33:
            dJydsigmay[33] = 1.0/sigmay33 - 1.0*pow(-my33 + y33, 2)/pow(sigmay33, 3);
            break;
        case 34:
            dJydsigmay[34] = 1.0/sigmay34 - 1.0*pow(-my34 + y34, 2)/pow(sigmay34, 3);
            break;
        case 35:
            dJydsigmay[35] = 1.0/sigmay35 - 1.0*pow(-my35 + y35, 2)/pow(sigmay35, 3);
            break;
        case 36:
            dJydsigmay[36] = 1.0/sigmay36 - 1.0*pow(-my36 + y36, 2)/pow(sigmay36, 3);
            break;
        case 37:
            dJydsigmay[37] = 1.0/sigmay37 - 1.0*pow(-my37 + y37, 2)/pow(sigmay37, 3);
            break;
        case 38:
            dJydsigmay[38] = 1.0/sigmay38 - 1.0*pow(-my38 + y38, 2)/pow(sigmay38, 3);
            break;
        case 39:
            dJydsigmay[39] = 1.0/sigmay39 - 1.0*pow(-my39 + y39, 2)/pow(sigmay39, 3);
            break;
        case 40:
            dJydsigmay[40] = 1.0/sigmay40 - 1.0*pow(-my40 + y40, 2)/pow(sigmay40, 3);
            break;
        case 41:
            dJydsigmay[41] = 1.0/sigmay41 - 1.0*pow(-my41 + y41, 2)/pow(sigmay41, 3);
            break;
        case 42:
            dJydsigmay[42] = 1.0/sigmay42 - 1.0*pow(-my42 + y42, 2)/pow(sigmay42, 3);
            break;
        case 43:
            dJydsigmay[43] = 1.0/sigmay43 - 1.0*pow(-my43 + y43, 2)/pow(sigmay43, 3);
            break;
        case 44:
            dJydsigmay[44] = 1.0/sigmay44 - 1.0*pow(-my44 + y44, 2)/pow(sigmay44, 3);
            break;
        case 45:
            dJydsigmay[45] = 1.0/sigmay45 - 1.0*pow(-my45 + y45, 2)/pow(sigmay45, 3);
            break;
        case 46:
            dJydsigmay[46] = 1.0/sigmay46 - 1.0*pow(-my46 + y46, 2)/pow(sigmay46, 3);
            break;
        case 47:
            dJydsigmay[47] = 1.0/sigmay47 - 1.0*pow(-my47 + y47, 2)/pow(sigmay47, 3);
            break;
        case 48:
            dJydsigmay[48] = 1.0/sigmay48 - 1.0*pow(-my48 + y48, 2)/pow(sigmay48, 3);
            break;
        case 49:
            dJydsigmay[49] = 1.0/sigmay49 - 1.0*pow(-my49 + y49, 2)/pow(sigmay49, 3);
            break;
        case 50:
            dJydsigmay[50] = 1.0/sigmay50 - 1.0*pow(-my50 + y50, 2)/pow(sigmay50, 3);
            break;
        case 51:
            dJydsigmay[51] = 1.0/sigmay51 - 1.0*pow(-my51 + y51, 2)/pow(sigmay51, 3);
            break;
        case 52:
            dJydsigmay[52] = 1.0/sigmay52 - 1.0*pow(-my52 + y52, 2)/pow(sigmay52, 3);
            break;
        case 53:
            dJydsigmay[53] = 1.0/sigmay53 - 1.0*pow(-my53 + y53, 2)/pow(sigmay53, 3);
            break;
        case 54:
            dJydsigmay[54] = 1.0/sigmay54 - 1.0*pow(-my54 + y54, 2)/pow(sigmay54, 3);
            break;
        case 55:
            dJydsigmay[55] = 1.0/sigmay55 - 1.0*pow(-my55 + y55, 2)/pow(sigmay55, 3);
            break;
        case 56:
            dJydsigmay[56] = 1.0/sigmay56 - 1.0*pow(-my56 + y56, 2)/pow(sigmay56, 3);
            break;
        case 57:
            dJydsigmay[57] = 1.0/sigmay57 - 1.0*pow(-my57 + y57, 2)/pow(sigmay57, 3);
            break;
        case 58:
            dJydsigmay[58] = 1.0/sigmay58 - 1.0*pow(-my58 + y58, 2)/pow(sigmay58, 3);
            break;
        case 59:
            dJydsigmay[59] = 1.0/sigmay59 - 1.0*pow(-my59 + y59, 2)/pow(sigmay59, 3);
            break;
        case 60:
            dJydsigmay[60] = 1.0/sigmay60 - 1.0*pow(-my60 + y60, 2)/pow(sigmay60, 3);
            break;
        case 61:
            dJydsigmay[61] = 1.0/sigmay61 - 1.0*pow(-my61 + y61, 2)/pow(sigmay61, 3);
            break;
        case 62:
            dJydsigmay[62] = 1.0/sigmay62 - 1.0*pow(-my62 + y62, 2)/pow(sigmay62, 3);
            break;
        case 63:
            dJydsigmay[63] = 1.0/sigmay63 - 1.0*pow(-my63 + y63, 2)/pow(sigmay63, 3);
            break;
        case 64:
            dJydsigmay[64] = 1.0/sigmay64 - 1.0*pow(-my64 + y64, 2)/pow(sigmay64, 3);
            break;
        case 65:
            dJydsigmay[65] = 1.0/sigmay65 - 1.0*pow(-my65 + y65, 2)/pow(sigmay65, 3);
            break;
        case 66:
            dJydsigmay[66] = 1.0/sigmay66 - 1.0*pow(-my66 + y66, 2)/pow(sigmay66, 3);
            break;
        case 67:
            dJydsigmay[67] = 1.0/sigmay67 - 1.0*pow(-my67 + y67, 2)/pow(sigmay67, 3);
            break;
        case 68:
            dJydsigmay[68] = 1.0/sigmay68 - 1.0*pow(-my68 + y68, 2)/pow(sigmay68, 3);
            break;
        case 69:
            dJydsigmay[69] = 1.0/sigmay69 - 1.0*pow(-my69 + y69, 2)/pow(sigmay69, 3);
            break;
        case 70:
            dJydsigmay[70] = 1.0/sigmay70 - 1.0*pow(-my70 + y70, 2)/pow(sigmay70, 3);
            break;
        case 71:
            dJydsigmay[71] = 1.0/sigmay71 - 1.0*pow(-my71 + y71, 2)/pow(sigmay71, 3);
            break;
        case 72:
            dJydsigmay[72] = 1.0/sigmay72 - 1.0*pow(-my72 + y72, 2)/pow(sigmay72, 3);
            break;
        case 73:
            dJydsigmay[73] = 1.0/sigmay73 - 1.0*pow(-my73 + y73, 2)/pow(sigmay73, 3);
            break;
        case 74:
            dJydsigmay[74] = 1.0/sigmay74 - 1.0*pow(-my74 + y74, 2)/pow(sigmay74, 3);
            break;
        case 75:
            dJydsigmay[75] = 1.0/sigmay75 - 1.0*pow(-my75 + y75, 2)/pow(sigmay75, 3);
            break;
        case 76:
            dJydsigmay[76] = 1.0/sigmay76 - 1.0*pow(-my76 + y76, 2)/pow(sigmay76, 3);
            break;
        case 77:
            dJydsigmay[77] = 1.0/sigmay77 - 1.0*pow(-my77 + y77, 2)/pow(sigmay77, 3);
            break;
        case 78:
            dJydsigmay[78] = 1.0/sigmay78 - 1.0*pow(-my78 + y78, 2)/pow(sigmay78, 3);
            break;
        case 79:
            dJydsigmay[79] = 1.0/sigmay79 - 1.0*pow(-my79 + y79, 2)/pow(sigmay79, 3);
            break;
        case 80:
            dJydsigmay[80] = 1.0/sigmay80 - 1.0*pow(-my80 + y80, 2)/pow(sigmay80, 3);
            break;
        case 81:
            dJydsigmay[81] = 1.0/sigmay81 - 1.0*pow(-my81 + y81, 2)/pow(sigmay81, 3);
            break;
        case 82:
            dJydsigmay[82] = 1.0/sigmay82 - 1.0*pow(-my82 + y82, 2)/pow(sigmay82, 3);
            break;
        case 83:
            dJydsigmay[83] = 1.0/sigmay83 - 1.0*pow(-my83 + y83, 2)/pow(sigmay83, 3);
            break;
        case 84:
            dJydsigmay[84] = 1.0/sigmay84 - 1.0*pow(-my84 + y84, 2)/pow(sigmay84, 3);
            break;
        case 85:
            dJydsigmay[85] = 1.0/sigmay85 - 1.0*pow(-my85 + y85, 2)/pow(sigmay85, 3);
            break;
        case 86:
            dJydsigmay[86] = 1.0/sigmay86 - 1.0*pow(-my86 + y86, 2)/pow(sigmay86, 3);
            break;
        case 87:
            dJydsigmay[87] = 1.0/sigmay87 - 1.0*pow(-my87 + y87, 2)/pow(sigmay87, 3);
            break;
        case 88:
            dJydsigmay[88] = 1.0/sigmay88 - 1.0*pow(-my88 + y88, 2)/pow(sigmay88, 3);
            break;
        case 89:
            dJydsigmay[89] = 1.0/sigmay89 - 1.0*pow(-my89 + y89, 2)/pow(sigmay89, 3);
            break;
        case 90:
            dJydsigmay[90] = 1.0/sigmay90 - 1.0*pow(-my90 + y90, 2)/pow(sigmay90, 3);
            break;
        case 91:
            dJydsigmay[91] = 1.0/sigmay91 - 1.0*pow(-my91 + y91, 2)/pow(sigmay91, 3);
            break;
        case 92:
            dJydsigmay[92] = 1.0/sigmay92 - 1.0*pow(-my92 + y92, 2)/pow(sigmay92, 3);
            break;
        case 93:
            dJydsigmay[93] = 1.0/sigmay93 - 1.0*pow(-my93 + y93, 2)/pow(sigmay93, 3);
            break;
        case 94:
            dJydsigmay[94] = 1.0/sigmay94 - 1.0*pow(-my94 + y94, 2)/pow(sigmay94, 3);
            break;
        case 95:
            dJydsigmay[95] = 1.0/sigmay95 - 1.0*pow(-my95 + y95, 2)/pow(sigmay95, 3);
            break;
        case 96:
            dJydsigmay[96] = 1.0/sigmay96 - 1.0*pow(-my96 + y96, 2)/pow(sigmay96, 3);
            break;
        case 97:
            dJydsigmay[97] = 1.0/sigmay97 - 1.0*pow(-my97 + y97, 2)/pow(sigmay97, 3);
            break;
        case 98:
            dJydsigmay[98] = 1.0/sigmay98 - 1.0*pow(-my98 + y98, 2)/pow(sigmay98, 3);
            break;
        case 99:
            dJydsigmay[99] = 1.0/sigmay99 - 1.0*pow(-my99 + y99, 2)/pow(sigmay99, 3);
            break;
        case 100:
            dJydsigmay[100] = 1.0/sigmay100 - 1.0*pow(-my100 + y100, 2)/pow(sigmay100, 3);
            break;
        case 101:
            dJydsigmay[101] = 1.0/sigmay101 - 1.0*pow(-my101 + y101, 2)/pow(sigmay101, 3);
            break;
        case 102:
            dJydsigmay[102] = 1.0/sigmay102 - 1.0*pow(-my102 + y102, 2)/pow(sigmay102, 3);
            break;
        case 103:
            dJydsigmay[103] = 1.0/sigmay103 - 1.0*pow(-my103 + y103, 2)/pow(sigmay103, 3);
            break;
        case 104:
            dJydsigmay[104] = 1.0/sigmay104 - 1.0*pow(-my104 + y104, 2)/pow(sigmay104, 3);
            break;
        case 105:
            dJydsigmay[105] = 1.0/sigmay105 - 1.0*pow(-my105 + y105, 2)/pow(sigmay105, 3);
            break;
        case 106:
            dJydsigmay[106] = 1.0/sigmay106 - 1.0*pow(-my106 + y106, 2)/pow(sigmay106, 3);
            break;
        case 107:
            dJydsigmay[107] = 1.0/sigmay107 - 1.0*pow(-my107 + y107, 2)/pow(sigmay107, 3);
            break;
        case 108:
            dJydsigmay[108] = 1.0/sigmay108 - 1.0*pow(-my108 + y108, 2)/pow(sigmay108, 3);
            break;
        case 109:
            dJydsigmay[109] = 1.0/sigmay109 - 1.0*pow(-my109 + y109, 2)/pow(sigmay109, 3);
            break;
        case 110:
            dJydsigmay[110] = 1.0/sigmay110 - 1.0*pow(-my110 + y110, 2)/pow(sigmay110, 3);
            break;
        case 111:
            dJydsigmay[111] = 1.0/sigmay111 - 1.0*pow(-my111 + y111, 2)/pow(sigmay111, 3);
            break;
        case 112:
            dJydsigmay[112] = 1.0/sigmay112 - 1.0*pow(-my112 + y112, 2)/pow(sigmay112, 3);
            break;
        case 113:
            dJydsigmay[113] = 1.0/sigmay113 - 1.0*pow(-my113 + y113, 2)/pow(sigmay113, 3);
            break;
        case 114:
            dJydsigmay[114] = 1.0/sigmay114 - 1.0*pow(-my114 + y114, 2)/pow(sigmay114, 3);
            break;
        case 115:
            dJydsigmay[115] = 1.0/sigmay115 - 1.0*pow(-my115 + y115, 2)/pow(sigmay115, 3);
            break;
        case 116:
            dJydsigmay[116] = 1.0/sigmay116 - 1.0*pow(-my116 + y116, 2)/pow(sigmay116, 3);
            break;
        case 117:
            dJydsigmay[117] = 1.0/sigmay117 - 1.0*pow(-my117 + y117, 2)/pow(sigmay117, 3);
            break;
        case 118:
            dJydsigmay[118] = 1.0/sigmay118 - 1.0*pow(-my118 + y118, 2)/pow(sigmay118, 3);
            break;
        case 119:
            dJydsigmay[119] = 1.0/sigmay119 - 1.0*pow(-my119 + y119, 2)/pow(sigmay119, 3);
            break;
        case 120:
            dJydsigmay[120] = 1.0/sigmay120 - 1.0*pow(-my120 + y120, 2)/pow(sigmay120, 3);
            break;
        case 121:
            dJydsigmay[121] = 1.0/sigmay121 - 1.0*pow(-my121 + y121, 2)/pow(sigmay121, 3);
            break;
        case 122:
            dJydsigmay[122] = 1.0/sigmay122 - 1.0*pow(-my122 + y122, 2)/pow(sigmay122, 3);
            break;
        case 123:
            dJydsigmay[123] = 1.0/sigmay123 - 1.0*pow(-my123 + y123, 2)/pow(sigmay123, 3);
            break;
        case 124:
            dJydsigmay[124] = 1.0/sigmay124 - 1.0*pow(-my124 + y124, 2)/pow(sigmay124, 3);
            break;
        case 125:
            dJydsigmay[125] = 1.0/sigmay125 - 1.0*pow(-my125 + y125, 2)/pow(sigmay125, 3);
            break;
        case 126:
            dJydsigmay[126] = 1.0/sigmay126 - 1.0*pow(-my126 + y126, 2)/pow(sigmay126, 3);
            break;
        case 127:
            dJydsigmay[127] = 1.0/sigmay127 - 1.0*pow(-my127 + y127, 2)/pow(sigmay127, 3);
            break;
        case 128:
            dJydsigmay[128] = 1.0/sigmay128 - 1.0*pow(-my128 + y128, 2)/pow(sigmay128, 3);
            break;
        case 129:
            dJydsigmay[129] = 1.0/sigmay129 - 1.0*pow(-my129 + y129, 2)/pow(sigmay129, 3);
            break;
        case 130:
            dJydsigmay[130] = 1.0/sigmay130 - 1.0*pow(-my130 + y130, 2)/pow(sigmay130, 3);
            break;
        case 131:
            dJydsigmay[131] = 1.0/sigmay131 - 1.0*pow(-my131 + y131, 2)/pow(sigmay131, 3);
            break;
        case 132:
            dJydsigmay[132] = 1.0/sigmay132 - 1.0*pow(-my132 + y132, 2)/pow(sigmay132, 3);
            break;
        case 133:
            dJydsigmay[133] = 1.0/sigmay133 - 1.0*pow(-my133 + y133, 2)/pow(sigmay133, 3);
            break;
        case 134:
            dJydsigmay[134] = 1.0/sigmay134 - 1.0*pow(-my134 + y134, 2)/pow(sigmay134, 3);
            break;
        case 135:
            dJydsigmay[135] = 1.0/sigmay135 - 1.0*pow(-my135 + y135, 2)/pow(sigmay135, 3);
            break;
        case 136:
            dJydsigmay[136] = 1.0/sigmay136 - 1.0*pow(-my136 + y136, 2)/pow(sigmay136, 3);
            break;
        case 137:
            dJydsigmay[137] = 1.0/sigmay137 - 1.0*pow(-my137 + y137, 2)/pow(sigmay137, 3);
            break;
        case 138:
            dJydsigmay[138] = 1.0/sigmay138 - 1.0*pow(-my138 + y138, 2)/pow(sigmay138, 3);
            break;
        case 139:
            dJydsigmay[139] = 1.0/sigmay139 - 1.0*pow(-my139 + y139, 2)/pow(sigmay139, 3);
            break;
        case 140:
            dJydsigmay[140] = 1.0/sigmay140 - 1.0*pow(-my140 + y140, 2)/pow(sigmay140, 3);
            break;
        case 141:
            dJydsigmay[141] = 1.0/sigmay141 - 1.0*pow(-my141 + y141, 2)/pow(sigmay141, 3);
            break;
        case 142:
            dJydsigmay[142] = 1.0/sigmay142 - 1.0*pow(-my142 + y142, 2)/pow(sigmay142, 3);
            break;
        case 143:
            dJydsigmay[143] = 1.0/sigmay143 - 1.0*pow(-my143 + y143, 2)/pow(sigmay143, 3);
            break;
        case 144:
            dJydsigmay[144] = 1.0/sigmay144 - 1.0*pow(-my144 + y144, 2)/pow(sigmay144, 3);
            break;
        case 145:
            dJydsigmay[145] = 1.0/sigmay145 - 1.0*pow(-my145 + y145, 2)/pow(sigmay145, 3);
            break;
        case 146:
            dJydsigmay[146] = 1.0/sigmay146 - 1.0*pow(-my146 + y146, 2)/pow(sigmay146, 3);
            break;
        case 147:
            dJydsigmay[147] = 1.0/sigmay147 - 1.0*pow(-my147 + y147, 2)/pow(sigmay147, 3);
            break;
        case 148:
            dJydsigmay[148] = 1.0/sigmay148 - 1.0*pow(-my148 + y148, 2)/pow(sigmay148, 3);
            break;
        case 149:
            dJydsigmay[149] = 1.0/sigmay149 - 1.0*pow(-my149 + y149, 2)/pow(sigmay149, 3);
            break;
        case 150:
            dJydsigmay[150] = 1.0/sigmay150 - 1.0*pow(-my150 + y150, 2)/pow(sigmay150, 3);
            break;
        case 151:
            dJydsigmay[151] = 1.0/sigmay151 - 1.0*pow(-my151 + y151, 2)/pow(sigmay151, 3);
            break;
        case 152:
            dJydsigmay[152] = 1.0/sigmay152 - 1.0*pow(-my152 + y152, 2)/pow(sigmay152, 3);
            break;
        case 153:
            dJydsigmay[153] = 1.0/sigmay153 - 1.0*pow(-my153 + y153, 2)/pow(sigmay153, 3);
            break;
        case 154:
            dJydsigmay[154] = 1.0/sigmay154 - 1.0*pow(-my154 + y154, 2)/pow(sigmay154, 3);
            break;
        case 155:
            dJydsigmay[155] = 1.0/sigmay155 - 1.0*pow(-my155 + y155, 2)/pow(sigmay155, 3);
            break;
        case 156:
            dJydsigmay[156] = 1.0/sigmay156 - 1.0*pow(-my156 + y156, 2)/pow(sigmay156, 3);
            break;
        case 157:
            dJydsigmay[157] = 1.0/sigmay157 - 1.0*pow(-my157 + y157, 2)/pow(sigmay157, 3);
            break;
        case 158:
            dJydsigmay[158] = 1.0/sigmay158 - 1.0*pow(-my158 + y158, 2)/pow(sigmay158, 3);
            break;
        case 159:
            dJydsigmay[159] = 1.0/sigmay159 - 1.0*pow(-my159 + y159, 2)/pow(sigmay159, 3);
            break;
        case 160:
            dJydsigmay[160] = 1.0/sigmay160 - 1.0*pow(-my160 + y160, 2)/pow(sigmay160, 3);
            break;
        case 161:
            dJydsigmay[161] = 1.0/sigmay161 - 1.0*pow(-my161 + y161, 2)/pow(sigmay161, 3);
            break;
        case 162:
            dJydsigmay[162] = 1.0/sigmay162 - 1.0*pow(-my162 + y162, 2)/pow(sigmay162, 3);
            break;
        case 163:
            dJydsigmay[163] = 1.0/sigmay163 - 1.0*pow(-my163 + y163, 2)/pow(sigmay163, 3);
            break;
        case 164:
            dJydsigmay[164] = 1.0/sigmay164 - 1.0*pow(-my164 + y164, 2)/pow(sigmay164, 3);
            break;
        case 165:
            dJydsigmay[165] = 1.0/sigmay165 - 1.0*pow(-my165 + y165, 2)/pow(sigmay165, 3);
            break;
        case 166:
            dJydsigmay[166] = 1.0/sigmay166 - 1.0*pow(-my166 + y166, 2)/pow(sigmay166, 3);
            break;
        case 167:
            dJydsigmay[167] = 1.0/sigmay167 - 1.0*pow(-my167 + y167, 2)/pow(sigmay167, 3);
            break;
        case 168:
            dJydsigmay[168] = 1.0/sigmay168 - 1.0*pow(-my168 + y168, 2)/pow(sigmay168, 3);
            break;
        case 169:
            dJydsigmay[169] = 1.0/sigmay169 - 1.0*pow(-my169 + y169, 2)/pow(sigmay169, 3);
            break;
        case 170:
            dJydsigmay[170] = 1.0/sigmay170 - 1.0*pow(-my170 + y170, 2)/pow(sigmay170, 3);
            break;
        case 171:
            dJydsigmay[171] = 1.0/sigmay171 - 1.0*pow(-my171 + y171, 2)/pow(sigmay171, 3);
            break;
        case 172:
            dJydsigmay[172] = 1.0/sigmay172 - 1.0*pow(-my172 + y172, 2)/pow(sigmay172, 3);
            break;
        case 173:
            dJydsigmay[173] = 1.0/sigmay173 - 1.0*pow(-my173 + y173, 2)/pow(sigmay173, 3);
            break;
        case 174:
            dJydsigmay[174] = 1.0/sigmay174 - 1.0*pow(-my174 + y174, 2)/pow(sigmay174, 3);
            break;
        case 175:
            dJydsigmay[175] = 1.0/sigmay175 - 1.0*pow(-my175 + y175, 2)/pow(sigmay175, 3);
            break;
        case 176:
            dJydsigmay[176] = 1.0/sigmay176 - 1.0*pow(-my176 + y176, 2)/pow(sigmay176, 3);
            break;
        case 177:
            dJydsigmay[177] = 1.0/sigmay177 - 1.0*pow(-my177 + y177, 2)/pow(sigmay177, 3);
            break;
        case 178:
            dJydsigmay[178] = 1.0/sigmay178 - 1.0*pow(-my178 + y178, 2)/pow(sigmay178, 3);
            break;
        case 179:
            dJydsigmay[179] = 1.0/sigmay179 - 1.0*pow(-my179 + y179, 2)/pow(sigmay179, 3);
            break;
        case 180:
            dJydsigmay[180] = 1.0/sigmay180 - 1.0*pow(-my180 + y180, 2)/pow(sigmay180, 3);
            break;
        case 181:
            dJydsigmay[181] = 1.0/sigmay181 - 1.0*pow(-my181 + y181, 2)/pow(sigmay181, 3);
            break;
        case 182:
            dJydsigmay[182] = 1.0/sigmay182 - 1.0*pow(-my182 + y182, 2)/pow(sigmay182, 3);
            break;
        case 183:
            dJydsigmay[183] = 1.0/sigmay183 - 1.0*pow(-my183 + y183, 2)/pow(sigmay183, 3);
            break;
        case 184:
            dJydsigmay[184] = 1.0/sigmay184 - 1.0*pow(-my184 + y184, 2)/pow(sigmay184, 3);
            break;
        case 185:
            dJydsigmay[185] = 1.0/sigmay185 - 1.0*pow(-my185 + y185, 2)/pow(sigmay185, 3);
            break;
        case 186:
            dJydsigmay[186] = 1.0/sigmay186 - 1.0*pow(-my186 + y186, 2)/pow(sigmay186, 3);
            break;
        case 187:
            dJydsigmay[187] = 1.0/sigmay187 - 1.0*pow(-my187 + y187, 2)/pow(sigmay187, 3);
            break;
        case 188:
            dJydsigmay[188] = 1.0/sigmay188 - 1.0*pow(-my188 + y188, 2)/pow(sigmay188, 3);
            break;
        case 189:
            dJydsigmay[189] = 1.0/sigmay189 - 1.0*pow(-my189 + y189, 2)/pow(sigmay189, 3);
            break;
        case 190:
            dJydsigmay[190] = 1.0/sigmay190 - 1.0*pow(-my190 + y190, 2)/pow(sigmay190, 3);
            break;
        case 191:
            dJydsigmay[191] = 1.0/sigmay191 - 1.0*pow(-my191 + y191, 2)/pow(sigmay191, 3);
            break;
        case 192:
            dJydsigmay[192] = 1.0/sigmay192 - 1.0*pow(-my192 + y192, 2)/pow(sigmay192, 3);
            break;
        case 193:
            dJydsigmay[193] = 1.0/sigmay193 - 1.0*pow(-my193 + y193, 2)/pow(sigmay193, 3);
            break;
    }
}