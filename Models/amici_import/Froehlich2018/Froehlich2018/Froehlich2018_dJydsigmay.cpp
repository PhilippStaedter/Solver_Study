#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydsigmay_Froehlich2018(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
        case 194:
            dJydsigmay[194] = 1.0/sigmay194 - 1.0*pow(-my194 + y194, 2)/pow(sigmay194, 3);
            break;
        case 195:
            dJydsigmay[195] = 1.0/sigmay195 - 1.0*pow(-my195 + y195, 2)/pow(sigmay195, 3);
            break;
        case 196:
            dJydsigmay[196] = 1.0/sigmay196 - 1.0*pow(-my196 + y196, 2)/pow(sigmay196, 3);
            break;
        case 197:
            dJydsigmay[197] = 1.0/sigmay197 - 1.0*pow(-my197 + y197, 2)/pow(sigmay197, 3);
            break;
        case 198:
            dJydsigmay[198] = 1.0/sigmay198 - 1.0*pow(-my198 + y198, 2)/pow(sigmay198, 3);
            break;
        case 199:
            dJydsigmay[199] = 1.0/sigmay199 - 1.0*pow(-my199 + y199, 2)/pow(sigmay199, 3);
            break;
        case 200:
            dJydsigmay[200] = 1.0/sigmay200 - 1.0*pow(-my200 + y200, 2)/pow(sigmay200, 3);
            break;
        case 201:
            dJydsigmay[201] = 1.0/sigmay201 - 1.0*pow(-my201 + y201, 2)/pow(sigmay201, 3);
            break;
        case 202:
            dJydsigmay[202] = 1.0/sigmay202 - 1.0*pow(-my202 + y202, 2)/pow(sigmay202, 3);
            break;
        case 203:
            dJydsigmay[203] = 1.0/sigmay203 - 1.0*pow(-my203 + y203, 2)/pow(sigmay203, 3);
            break;
        case 204:
            dJydsigmay[204] = 1.0/sigmay204 - 1.0*pow(-my204 + y204, 2)/pow(sigmay204, 3);
            break;
        case 205:
            dJydsigmay[205] = 1.0/sigmay205 - 1.0*pow(-my205 + y205, 2)/pow(sigmay205, 3);
            break;
        case 206:
            dJydsigmay[206] = 1.0/sigmay206 - 1.0*pow(-my206 + y206, 2)/pow(sigmay206, 3);
            break;
        case 207:
            dJydsigmay[207] = 1.0/sigmay207 - 1.0*pow(-my207 + y207, 2)/pow(sigmay207, 3);
            break;
        case 208:
            dJydsigmay[208] = 1.0/sigmay208 - 1.0*pow(-my208 + y208, 2)/pow(sigmay208, 3);
            break;
        case 209:
            dJydsigmay[209] = 1.0/sigmay209 - 1.0*pow(-my209 + y209, 2)/pow(sigmay209, 3);
            break;
        case 210:
            dJydsigmay[210] = 1.0/sigmay210 - 1.0*pow(-my210 + y210, 2)/pow(sigmay210, 3);
            break;
        case 211:
            dJydsigmay[211] = 1.0/sigmay211 - 1.0*pow(-my211 + y211, 2)/pow(sigmay211, 3);
            break;
        case 212:
            dJydsigmay[212] = 1.0/sigmay212 - 1.0*pow(-my212 + y212, 2)/pow(sigmay212, 3);
            break;
        case 213:
            dJydsigmay[213] = 1.0/sigmay213 - 1.0*pow(-my213 + y213, 2)/pow(sigmay213, 3);
            break;
        case 214:
            dJydsigmay[214] = 1.0/sigmay214 - 1.0*pow(-my214 + y214, 2)/pow(sigmay214, 3);
            break;
        case 215:
            dJydsigmay[215] = 1.0/sigmay215 - 1.0*pow(-my215 + y215, 2)/pow(sigmay215, 3);
            break;
        case 216:
            dJydsigmay[216] = 1.0/sigmay216 - 1.0*pow(-my216 + y216, 2)/pow(sigmay216, 3);
            break;
        case 217:
            dJydsigmay[217] = 1.0/sigmay217 - 1.0*pow(-my217 + y217, 2)/pow(sigmay217, 3);
            break;
        case 218:
            dJydsigmay[218] = 1.0/sigmay218 - 1.0*pow(-my218 + y218, 2)/pow(sigmay218, 3);
            break;
        case 219:
            dJydsigmay[219] = 1.0/sigmay219 - 1.0*pow(-my219 + y219, 2)/pow(sigmay219, 3);
            break;
        case 220:
            dJydsigmay[220] = 1.0/sigmay220 - 1.0*pow(-my220 + y220, 2)/pow(sigmay220, 3);
            break;
        case 221:
            dJydsigmay[221] = 1.0/sigmay221 - 1.0*pow(-my221 + y221, 2)/pow(sigmay221, 3);
            break;
        case 222:
            dJydsigmay[222] = 1.0/sigmay222 - 1.0*pow(-my222 + y222, 2)/pow(sigmay222, 3);
            break;
        case 223:
            dJydsigmay[223] = 1.0/sigmay223 - 1.0*pow(-my223 + y223, 2)/pow(sigmay223, 3);
            break;
        case 224:
            dJydsigmay[224] = 1.0/sigmay224 - 1.0*pow(-my224 + y224, 2)/pow(sigmay224, 3);
            break;
        case 225:
            dJydsigmay[225] = 1.0/sigmay225 - 1.0*pow(-my225 + y225, 2)/pow(sigmay225, 3);
            break;
        case 226:
            dJydsigmay[226] = 1.0/sigmay226 - 1.0*pow(-my226 + y226, 2)/pow(sigmay226, 3);
            break;
        case 227:
            dJydsigmay[227] = 1.0/sigmay227 - 1.0*pow(-my227 + y227, 2)/pow(sigmay227, 3);
            break;
        case 228:
            dJydsigmay[228] = 1.0/sigmay228 - 1.0*pow(-my228 + y228, 2)/pow(sigmay228, 3);
            break;
        case 229:
            dJydsigmay[229] = 1.0/sigmay229 - 1.0*pow(-my229 + y229, 2)/pow(sigmay229, 3);
            break;
        case 230:
            dJydsigmay[230] = 1.0/sigmay230 - 1.0*pow(-my230 + y230, 2)/pow(sigmay230, 3);
            break;
        case 231:
            dJydsigmay[231] = 1.0/sigmay231 - 1.0*pow(-my231 + y231, 2)/pow(sigmay231, 3);
            break;
        case 232:
            dJydsigmay[232] = 1.0/sigmay232 - 1.0*pow(-my232 + y232, 2)/pow(sigmay232, 3);
            break;
        case 233:
            dJydsigmay[233] = 1.0/sigmay233 - 1.0*pow(-my233 + y233, 2)/pow(sigmay233, 3);
            break;
        case 234:
            dJydsigmay[234] = 1.0/sigmay234 - 1.0*pow(-my234 + y234, 2)/pow(sigmay234, 3);
            break;
        case 235:
            dJydsigmay[235] = 1.0/sigmay235 - 1.0*pow(-my235 + y235, 2)/pow(sigmay235, 3);
            break;
        case 236:
            dJydsigmay[236] = 1.0/sigmay236 - 1.0*pow(-my236 + y236, 2)/pow(sigmay236, 3);
            break;
        case 237:
            dJydsigmay[237] = 1.0/sigmay237 - 1.0*pow(-my237 + y237, 2)/pow(sigmay237, 3);
            break;
        case 238:
            dJydsigmay[238] = 1.0/sigmay238 - 1.0*pow(-my238 + y238, 2)/pow(sigmay238, 3);
            break;
        case 239:
            dJydsigmay[239] = 1.0/sigmay239 - 1.0*pow(-my239 + y239, 2)/pow(sigmay239, 3);
            break;
        case 240:
            dJydsigmay[240] = 1.0/sigmay240 - 1.0*pow(-my240 + y240, 2)/pow(sigmay240, 3);
            break;
        case 241:
            dJydsigmay[241] = 1.0/sigmay241 - 1.0*pow(-my241 + y241, 2)/pow(sigmay241, 3);
            break;
        case 242:
            dJydsigmay[242] = 1.0/sigmay242 - 1.0*pow(-my242 + y242, 2)/pow(sigmay242, 3);
            break;
        case 243:
            dJydsigmay[243] = 1.0/sigmay243 - 1.0*pow(-my243 + y243, 2)/pow(sigmay243, 3);
            break;
        case 244:
            dJydsigmay[244] = 1.0/sigmay244 - 1.0*pow(-my244 + y244, 2)/pow(sigmay244, 3);
            break;
        case 245:
            dJydsigmay[245] = 1.0/sigmay245 - 1.0*pow(-my245 + y245, 2)/pow(sigmay245, 3);
            break;
        case 246:
            dJydsigmay[246] = 1.0/sigmay246 - 1.0*pow(-my246 + y246, 2)/pow(sigmay246, 3);
            break;
        case 247:
            dJydsigmay[247] = 1.0/sigmay247 - 1.0*pow(-my247 + y247, 2)/pow(sigmay247, 3);
            break;
        case 248:
            dJydsigmay[248] = 1.0/sigmay248 - 1.0*pow(-my248 + y248, 2)/pow(sigmay248, 3);
            break;
        case 249:
            dJydsigmay[249] = 1.0/sigmay249 - 1.0*pow(-my249 + y249, 2)/pow(sigmay249, 3);
            break;
        case 250:
            dJydsigmay[250] = 1.0/sigmay250 - 1.0*pow(-my250 + y250, 2)/pow(sigmay250, 3);
            break;
        case 251:
            dJydsigmay[251] = 1.0/sigmay251 - 1.0*pow(-my251 + y251, 2)/pow(sigmay251, 3);
            break;
        case 252:
            dJydsigmay[252] = 1.0/sigmay252 - 1.0*pow(-my252 + y252, 2)/pow(sigmay252, 3);
            break;
        case 253:
            dJydsigmay[253] = 1.0/sigmay253 - 1.0*pow(-my253 + y253, 2)/pow(sigmay253, 3);
            break;
        case 254:
            dJydsigmay[254] = 1.0/sigmay254 - 1.0*pow(-my254 + y254, 2)/pow(sigmay254, 3);
            break;
        case 255:
            dJydsigmay[255] = 1.0/sigmay255 - 1.0*pow(-my255 + y255, 2)/pow(sigmay255, 3);
            break;
        case 256:
            dJydsigmay[256] = 1.0/sigmay256 - 1.0*pow(-my256 + y256, 2)/pow(sigmay256, 3);
            break;
        case 257:
            dJydsigmay[257] = 1.0/sigmay257 - 1.0*pow(-my257 + y257, 2)/pow(sigmay257, 3);
            break;
        case 258:
            dJydsigmay[258] = 1.0/sigmay258 - 1.0*pow(-my258 + y258, 2)/pow(sigmay258, 3);
            break;
        case 259:
            dJydsigmay[259] = 1.0/sigmay259 - 1.0*pow(-my259 + y259, 2)/pow(sigmay259, 3);
            break;
        case 260:
            dJydsigmay[260] = 1.0/sigmay260 - 1.0*pow(-my260 + y260, 2)/pow(sigmay260, 3);
            break;
        case 261:
            dJydsigmay[261] = 1.0/sigmay261 - 1.0*pow(-my261 + y261, 2)/pow(sigmay261, 3);
            break;
        case 262:
            dJydsigmay[262] = 1.0/sigmay262 - 1.0*pow(-my262 + y262, 2)/pow(sigmay262, 3);
            break;
        case 263:
            dJydsigmay[263] = 1.0/sigmay263 - 1.0*pow(-my263 + y263, 2)/pow(sigmay263, 3);
            break;
        case 264:
            dJydsigmay[264] = 1.0/sigmay264 - 1.0*pow(-my264 + y264, 2)/pow(sigmay264, 3);
            break;
        case 265:
            dJydsigmay[265] = 1.0/sigmay265 - 1.0*pow(-my265 + y265, 2)/pow(sigmay265, 3);
            break;
        case 266:
            dJydsigmay[266] = 1.0/sigmay266 - 1.0*pow(-my266 + y266, 2)/pow(sigmay266, 3);
            break;
        case 267:
            dJydsigmay[267] = 1.0/sigmay267 - 1.0*pow(-my267 + y267, 2)/pow(sigmay267, 3);
            break;
        case 268:
            dJydsigmay[268] = 1.0/sigmay268 - 1.0*pow(-my268 + y268, 2)/pow(sigmay268, 3);
            break;
        case 269:
            dJydsigmay[269] = 1.0/sigmay269 - 1.0*pow(-my269 + y269, 2)/pow(sigmay269, 3);
            break;
        case 270:
            dJydsigmay[270] = 1.0/sigmay270 - 1.0*pow(-my270 + y270, 2)/pow(sigmay270, 3);
            break;
        case 271:
            dJydsigmay[271] = 1.0/sigmay271 - 1.0*pow(-my271 + y271, 2)/pow(sigmay271, 3);
            break;
        case 272:
            dJydsigmay[272] = 1.0/sigmay272 - 1.0*pow(-my272 + y272, 2)/pow(sigmay272, 3);
            break;
        case 273:
            dJydsigmay[273] = 1.0/sigmay273 - 1.0*pow(-my273 + y273, 2)/pow(sigmay273, 3);
            break;
        case 274:
            dJydsigmay[274] = 1.0/sigmay274 - 1.0*pow(-my274 + y274, 2)/pow(sigmay274, 3);
            break;
        case 275:
            dJydsigmay[275] = 1.0/sigmay275 - 1.0*pow(-my275 + y275, 2)/pow(sigmay275, 3);
            break;
        case 276:
            dJydsigmay[276] = 1.0/sigmay276 - 1.0*pow(-my276 + y276, 2)/pow(sigmay276, 3);
            break;
        case 277:
            dJydsigmay[277] = 1.0/sigmay277 - 1.0*pow(-my277 + y277, 2)/pow(sigmay277, 3);
            break;
        case 278:
            dJydsigmay[278] = 1.0/sigmay278 - 1.0*pow(-my278 + y278, 2)/pow(sigmay278, 3);
            break;
        case 279:
            dJydsigmay[279] = 1.0/sigmay279 - 1.0*pow(-my279 + y279, 2)/pow(sigmay279, 3);
            break;
        case 280:
            dJydsigmay[280] = 1.0/sigmay280 - 1.0*pow(-my280 + y280, 2)/pow(sigmay280, 3);
            break;
        case 281:
            dJydsigmay[281] = 1.0/sigmay281 - 1.0*pow(-my281 + y281, 2)/pow(sigmay281, 3);
            break;
        case 282:
            dJydsigmay[282] = 1.0/sigmay282 - 1.0*pow(-my282 + y282, 2)/pow(sigmay282, 3);
            break;
        case 283:
            dJydsigmay[283] = 1.0/sigmay283 - 1.0*pow(-my283 + y283, 2)/pow(sigmay283, 3);
            break;
        case 284:
            dJydsigmay[284] = 1.0/sigmay284 - 1.0*pow(-my284 + y284, 2)/pow(sigmay284, 3);
            break;
        case 285:
            dJydsigmay[285] = 1.0/sigmay285 - 1.0*pow(-my285 + y285, 2)/pow(sigmay285, 3);
            break;
        case 286:
            dJydsigmay[286] = 1.0/sigmay286 - 1.0*pow(-my286 + y286, 2)/pow(sigmay286, 3);
            break;
        case 287:
            dJydsigmay[287] = 1.0/sigmay287 - 1.0*pow(-my287 + y287, 2)/pow(sigmay287, 3);
            break;
        case 288:
            dJydsigmay[288] = 1.0/sigmay288 - 1.0*pow(-my288 + y288, 2)/pow(sigmay288, 3);
            break;
        case 289:
            dJydsigmay[289] = 1.0/sigmay289 - 1.0*pow(-my289 + y289, 2)/pow(sigmay289, 3);
            break;
        case 290:
            dJydsigmay[290] = 1.0/sigmay290 - 1.0*pow(-my290 + y290, 2)/pow(sigmay290, 3);
            break;
        case 291:
            dJydsigmay[291] = 1.0/sigmay291 - 1.0*pow(-my291 + y291, 2)/pow(sigmay291, 3);
            break;
        case 292:
            dJydsigmay[292] = 1.0/sigmay292 - 1.0*pow(-my292 + y292, 2)/pow(sigmay292, 3);
            break;
        case 293:
            dJydsigmay[293] = 1.0/sigmay293 - 1.0*pow(-my293 + y293, 2)/pow(sigmay293, 3);
            break;
        case 294:
            dJydsigmay[294] = 1.0/sigmay294 - 1.0*pow(-my294 + y294, 2)/pow(sigmay294, 3);
            break;
        case 295:
            dJydsigmay[295] = 1.0/sigmay295 - 1.0*pow(-my295 + y295, 2)/pow(sigmay295, 3);
            break;
        case 296:
            dJydsigmay[296] = 1.0/sigmay296 - 1.0*pow(-my296 + y296, 2)/pow(sigmay296, 3);
            break;
        case 297:
            dJydsigmay[297] = 1.0/sigmay297 - 1.0*pow(-my297 + y297, 2)/pow(sigmay297, 3);
            break;
        case 298:
            dJydsigmay[298] = 1.0/sigmay298 - 1.0*pow(-my298 + y298, 2)/pow(sigmay298, 3);
            break;
        case 299:
            dJydsigmay[299] = 1.0/sigmay299 - 1.0*pow(-my299 + y299, 2)/pow(sigmay299, 3);
            break;
        case 300:
            dJydsigmay[300] = 1.0/sigmay300 - 1.0*pow(-my300 + y300, 2)/pow(sigmay300, 3);
            break;
        case 301:
            dJydsigmay[301] = 1.0/sigmay301 - 1.0*pow(-my301 + y301, 2)/pow(sigmay301, 3);
            break;
        case 302:
            dJydsigmay[302] = 1.0/sigmay302 - 1.0*pow(-my302 + y302, 2)/pow(sigmay302, 3);
            break;
        case 303:
            dJydsigmay[303] = 1.0/sigmay303 - 1.0*pow(-my303 + y303, 2)/pow(sigmay303, 3);
            break;
        case 304:
            dJydsigmay[304] = 1.0/sigmay304 - 1.0*pow(-my304 + y304, 2)/pow(sigmay304, 3);
            break;
        case 305:
            dJydsigmay[305] = 1.0/sigmay305 - 1.0*pow(-my305 + y305, 2)/pow(sigmay305, 3);
            break;
        case 306:
            dJydsigmay[306] = 1.0/sigmay306 - 1.0*pow(-my306 + y306, 2)/pow(sigmay306, 3);
            break;
        case 307:
            dJydsigmay[307] = 1.0/sigmay307 - 1.0*pow(-my307 + y307, 2)/pow(sigmay307, 3);
            break;
        case 308:
            dJydsigmay[308] = 1.0/sigmay308 - 1.0*pow(-my308 + y308, 2)/pow(sigmay308, 3);
            break;
        case 309:
            dJydsigmay[309] = 1.0/sigmay309 - 1.0*pow(-my309 + y309, 2)/pow(sigmay309, 3);
            break;
        case 310:
            dJydsigmay[310] = 1.0/sigmay310 - 1.0*pow(-my310 + y310, 2)/pow(sigmay310, 3);
            break;
        case 311:
            dJydsigmay[311] = 1.0/sigmay311 - 1.0*pow(-my311 + y311, 2)/pow(sigmay311, 3);
            break;
        case 312:
            dJydsigmay[312] = 1.0/sigmay312 - 1.0*pow(-my312 + y312, 2)/pow(sigmay312, 3);
            break;
        case 313:
            dJydsigmay[313] = 1.0/sigmay313 - 1.0*pow(-my313 + y313, 2)/pow(sigmay313, 3);
            break;
        case 314:
            dJydsigmay[314] = 1.0/sigmay314 - 1.0*pow(-my314 + y314, 2)/pow(sigmay314, 3);
            break;
        case 315:
            dJydsigmay[315] = 1.0/sigmay315 - 1.0*pow(-my315 + y315, 2)/pow(sigmay315, 3);
            break;
        case 316:
            dJydsigmay[316] = 1.0/sigmay316 - 1.0*pow(-my316 + y316, 2)/pow(sigmay316, 3);
            break;
        case 317:
            dJydsigmay[317] = 1.0/sigmay317 - 1.0*pow(-my317 + y317, 2)/pow(sigmay317, 3);
            break;
        case 318:
            dJydsigmay[318] = 1.0/sigmay318 - 1.0*pow(-my318 + y318, 2)/pow(sigmay318, 3);
            break;
        case 319:
            dJydsigmay[319] = 1.0/sigmay319 - 1.0*pow(-my319 + y319, 2)/pow(sigmay319, 3);
            break;
        case 320:
            dJydsigmay[320] = 1.0/sigmay320 - 1.0*pow(-my320 + y320, 2)/pow(sigmay320, 3);
            break;
        case 321:
            dJydsigmay[321] = 1.0/sigmay321 - 1.0*pow(-my321 + y321, 2)/pow(sigmay321, 3);
            break;
        case 322:
            dJydsigmay[322] = 1.0/sigmay322 - 1.0*pow(-my322 + y322, 2)/pow(sigmay322, 3);
            break;
        case 323:
            dJydsigmay[323] = 1.0/sigmay323 - 1.0*pow(-my323 + y323, 2)/pow(sigmay323, 3);
            break;
        case 324:
            dJydsigmay[324] = 1.0/sigmay324 - 1.0*pow(-my324 + y324, 2)/pow(sigmay324, 3);
            break;
        case 325:
            dJydsigmay[325] = 1.0/sigmay325 - 1.0*pow(-my325 + y325, 2)/pow(sigmay325, 3);
            break;
        case 326:
            dJydsigmay[326] = 1.0/sigmay326 - 1.0*pow(-my326 + y326, 2)/pow(sigmay326, 3);
            break;
        case 327:
            dJydsigmay[327] = 1.0/sigmay327 - 1.0*pow(-my327 + y327, 2)/pow(sigmay327, 3);
            break;
        case 328:
            dJydsigmay[328] = 1.0/sigmay328 - 1.0*pow(-my328 + y328, 2)/pow(sigmay328, 3);
            break;
        case 329:
            dJydsigmay[329] = 1.0/sigmay329 - 1.0*pow(-my329 + y329, 2)/pow(sigmay329, 3);
            break;
        case 330:
            dJydsigmay[330] = 1.0/sigmay330 - 1.0*pow(-my330 + y330, 2)/pow(sigmay330, 3);
            break;
        case 331:
            dJydsigmay[331] = 1.0/sigmay331 - 1.0*pow(-my331 + y331, 2)/pow(sigmay331, 3);
            break;
        case 332:
            dJydsigmay[332] = 1.0/sigmay332 - 1.0*pow(-my332 + y332, 2)/pow(sigmay332, 3);
            break;
        case 333:
            dJydsigmay[333] = 1.0/sigmay333 - 1.0*pow(-my333 + y333, 2)/pow(sigmay333, 3);
            break;
        case 334:
            dJydsigmay[334] = 1.0/sigmay334 - 1.0*pow(-my334 + y334, 2)/pow(sigmay334, 3);
            break;
        case 335:
            dJydsigmay[335] = 1.0/sigmay335 - 1.0*pow(-my335 + y335, 2)/pow(sigmay335, 3);
            break;
        case 336:
            dJydsigmay[336] = 1.0/sigmay336 - 1.0*pow(-my336 + y336, 2)/pow(sigmay336, 3);
            break;
        case 337:
            dJydsigmay[337] = 1.0/sigmay337 - 1.0*pow(-my337 + y337, 2)/pow(sigmay337, 3);
            break;
        case 338:
            dJydsigmay[338] = 1.0/sigmay338 - 1.0*pow(-my338 + y338, 2)/pow(sigmay338, 3);
            break;
        case 339:
            dJydsigmay[339] = 1.0/sigmay339 - 1.0*pow(-my339 + y339, 2)/pow(sigmay339, 3);
            break;
        case 340:
            dJydsigmay[340] = 1.0/sigmay340 - 1.0*pow(-my340 + y340, 2)/pow(sigmay340, 3);
            break;
        case 341:
            dJydsigmay[341] = 1.0/sigmay341 - 1.0*pow(-my341 + y341, 2)/pow(sigmay341, 3);
            break;
        case 342:
            dJydsigmay[342] = 1.0/sigmay342 - 1.0*pow(-my342 + y342, 2)/pow(sigmay342, 3);
            break;
        case 343:
            dJydsigmay[343] = 1.0/sigmay343 - 1.0*pow(-my343 + y343, 2)/pow(sigmay343, 3);
            break;
        case 344:
            dJydsigmay[344] = 1.0/sigmay344 - 1.0*pow(-my344 + y344, 2)/pow(sigmay344, 3);
            break;
        case 345:
            dJydsigmay[345] = 1.0/sigmay345 - 1.0*pow(-my345 + y345, 2)/pow(sigmay345, 3);
            break;
        case 346:
            dJydsigmay[346] = 1.0/sigmay346 - 1.0*pow(-my346 + y346, 2)/pow(sigmay346, 3);
            break;
        case 347:
            dJydsigmay[347] = 1.0/sigmay347 - 1.0*pow(-my347 + y347, 2)/pow(sigmay347, 3);
            break;
        case 348:
            dJydsigmay[348] = 1.0/sigmay348 - 1.0*pow(-my348 + y348, 2)/pow(sigmay348, 3);
            break;
        case 349:
            dJydsigmay[349] = 1.0/sigmay349 - 1.0*pow(-my349 + y349, 2)/pow(sigmay349, 3);
            break;
        case 350:
            dJydsigmay[350] = 1.0/sigmay350 - 1.0*pow(-my350 + y350, 2)/pow(sigmay350, 3);
            break;
        case 351:
            dJydsigmay[351] = 1.0/sigmay351 - 1.0*pow(-my351 + y351, 2)/pow(sigmay351, 3);
            break;
        case 352:
            dJydsigmay[352] = 1.0/sigmay352 - 1.0*pow(-my352 + y352, 2)/pow(sigmay352, 3);
            break;
        case 353:
            dJydsigmay[353] = 1.0/sigmay353 - 1.0*pow(-my353 + y353, 2)/pow(sigmay353, 3);
            break;
        case 354:
            dJydsigmay[354] = 1.0/sigmay354 - 1.0*pow(-my354 + y354, 2)/pow(sigmay354, 3);
            break;
        case 355:
            dJydsigmay[355] = 1.0/sigmay355 - 1.0*pow(-my355 + y355, 2)/pow(sigmay355, 3);
            break;
        case 356:
            dJydsigmay[356] = 1.0/sigmay356 - 1.0*pow(-my356 + y356, 2)/pow(sigmay356, 3);
            break;
        case 357:
            dJydsigmay[357] = 1.0/sigmay357 - 1.0*pow(-my357 + y357, 2)/pow(sigmay357, 3);
            break;
        case 358:
            dJydsigmay[358] = 1.0/sigmay358 - 1.0*pow(-my358 + y358, 2)/pow(sigmay358, 3);
            break;
        case 359:
            dJydsigmay[359] = 1.0/sigmay359 - 1.0*pow(-my359 + y359, 2)/pow(sigmay359, 3);
            break;
        case 360:
            dJydsigmay[360] = 1.0/sigmay360 - 1.0*pow(-my360 + y360, 2)/pow(sigmay360, 3);
            break;
        case 361:
            dJydsigmay[361] = 1.0/sigmay361 - 1.0*pow(-my361 + y361, 2)/pow(sigmay361, 3);
            break;
        case 362:
            dJydsigmay[362] = 1.0/sigmay362 - 1.0*pow(-my362 + y362, 2)/pow(sigmay362, 3);
            break;
        case 363:
            dJydsigmay[363] = 1.0/sigmay363 - 1.0*pow(-my363 + y363, 2)/pow(sigmay363, 3);
            break;
        case 364:
            dJydsigmay[364] = 1.0/sigmay364 - 1.0*pow(-my364 + y364, 2)/pow(sigmay364, 3);
            break;
        case 365:
            dJydsigmay[365] = 1.0/sigmay365 - 1.0*pow(-my365 + y365, 2)/pow(sigmay365, 3);
            break;
        case 366:
            dJydsigmay[366] = 1.0/sigmay366 - 1.0*pow(-my366 + y366, 2)/pow(sigmay366, 3);
            break;
        case 367:
            dJydsigmay[367] = 1.0/sigmay367 - 1.0*pow(-my367 + y367, 2)/pow(sigmay367, 3);
            break;
        case 368:
            dJydsigmay[368] = 1.0/sigmay368 - 1.0*pow(-my368 + y368, 2)/pow(sigmay368, 3);
            break;
        case 369:
            dJydsigmay[369] = 1.0/sigmay369 - 1.0*pow(-my369 + y369, 2)/pow(sigmay369, 3);
            break;
        case 370:
            dJydsigmay[370] = 1.0/sigmay370 - 1.0*pow(-my370 + y370, 2)/pow(sigmay370, 3);
            break;
        case 371:
            dJydsigmay[371] = 1.0/sigmay371 - 1.0*pow(-my371 + y371, 2)/pow(sigmay371, 3);
            break;
        case 372:
            dJydsigmay[372] = 1.0/sigmay372 - 1.0*pow(-my372 + y372, 2)/pow(sigmay372, 3);
            break;
        case 373:
            dJydsigmay[373] = 1.0/sigmay373 - 1.0*pow(-my373 + y373, 2)/pow(sigmay373, 3);
            break;
        case 374:
            dJydsigmay[374] = 1.0/sigmay374 - 1.0*pow(-my374 + y374, 2)/pow(sigmay374, 3);
            break;
        case 375:
            dJydsigmay[375] = 1.0/sigmay375 - 1.0*pow(-my375 + y375, 2)/pow(sigmay375, 3);
            break;
        case 376:
            dJydsigmay[376] = 1.0/sigmay376 - 1.0*pow(-my376 + y376, 2)/pow(sigmay376, 3);
            break;
        case 377:
            dJydsigmay[377] = 1.0/sigmay377 - 1.0*pow(-my377 + y377, 2)/pow(sigmay377, 3);
            break;
        case 378:
            dJydsigmay[378] = 1.0/sigmay378 - 1.0*pow(-my378 + y378, 2)/pow(sigmay378, 3);
            break;
        case 379:
            dJydsigmay[379] = 1.0/sigmay379 - 1.0*pow(-my379 + y379, 2)/pow(sigmay379, 3);
            break;
        case 380:
            dJydsigmay[380] = 1.0/sigmay380 - 1.0*pow(-my380 + y380, 2)/pow(sigmay380, 3);
            break;
        case 381:
            dJydsigmay[381] = 1.0/sigmay381 - 1.0*pow(-my381 + y381, 2)/pow(sigmay381, 3);
            break;
        case 382:
            dJydsigmay[382] = 1.0/sigmay382 - 1.0*pow(-my382 + y382, 2)/pow(sigmay382, 3);
            break;
        case 383:
            dJydsigmay[383] = 1.0/sigmay383 - 1.0*pow(-my383 + y383, 2)/pow(sigmay383, 3);
            break;
        case 384:
            dJydsigmay[384] = 1.0/sigmay384 - 1.0*pow(-my384 + y384, 2)/pow(sigmay384, 3);
            break;
        case 385:
            dJydsigmay[385] = 1.0/sigmay385 - 1.0*pow(-my385 + y385, 2)/pow(sigmay385, 3);
            break;
        case 386:
            dJydsigmay[386] = 1.0/sigmay386 - 1.0*pow(-my386 + y386, 2)/pow(sigmay386, 3);
            break;
        case 387:
            dJydsigmay[387] = 1.0/sigmay387 - 1.0*pow(-my387 + y387, 2)/pow(sigmay387, 3);
            break;
        case 388:
            dJydsigmay[388] = 1.0/sigmay388 - 1.0*pow(-my388 + y388, 2)/pow(sigmay388, 3);
            break;
        case 389:
            dJydsigmay[389] = 1.0/sigmay389 - 1.0*pow(-my389 + y389, 2)/pow(sigmay389, 3);
            break;
        case 390:
            dJydsigmay[390] = 1.0/sigmay390 - 1.0*pow(-my390 + y390, 2)/pow(sigmay390, 3);
            break;
        case 391:
            dJydsigmay[391] = 1.0/sigmay391 - 1.0*pow(-my391 + y391, 2)/pow(sigmay391, 3);
            break;
        case 392:
            dJydsigmay[392] = 1.0/sigmay392 - 1.0*pow(-my392 + y392, 2)/pow(sigmay392, 3);
            break;
        case 393:
            dJydsigmay[393] = 1.0/sigmay393 - 1.0*pow(-my393 + y393, 2)/pow(sigmay393, 3);
            break;
        case 394:
            dJydsigmay[394] = 1.0/sigmay394 - 1.0*pow(-my394 + y394, 2)/pow(sigmay394, 3);
            break;
        case 395:
            dJydsigmay[395] = 1.0/sigmay395 - 1.0*pow(-my395 + y395, 2)/pow(sigmay395, 3);
            break;
        case 396:
            dJydsigmay[396] = 1.0/sigmay396 - 1.0*pow(-my396 + y396, 2)/pow(sigmay396, 3);
            break;
        case 397:
            dJydsigmay[397] = 1.0/sigmay397 - 1.0*pow(-my397 + y397, 2)/pow(sigmay397, 3);
            break;
        case 398:
            dJydsigmay[398] = 1.0/sigmay398 - 1.0*pow(-my398 + y398, 2)/pow(sigmay398, 3);
            break;
        case 399:
            dJydsigmay[399] = 1.0/sigmay399 - 1.0*pow(-my399 + y399, 2)/pow(sigmay399, 3);
            break;
        case 400:
            dJydsigmay[400] = 1.0/sigmay400 - 1.0*pow(-my400 + y400, 2)/pow(sigmay400, 3);
            break;
        case 401:
            dJydsigmay[401] = 1.0/sigmay401 - 1.0*pow(-my401 + y401, 2)/pow(sigmay401, 3);
            break;
        case 402:
            dJydsigmay[402] = 1.0/sigmay402 - 1.0*pow(-my402 + y402, 2)/pow(sigmay402, 3);
            break;
        case 403:
            dJydsigmay[403] = 1.0/sigmay403 - 1.0*pow(-my403 + y403, 2)/pow(sigmay403, 3);
            break;
        case 404:
            dJydsigmay[404] = 1.0/sigmay404 - 1.0*pow(-my404 + y404, 2)/pow(sigmay404, 3);
            break;
        case 405:
            dJydsigmay[405] = 1.0/sigmay405 - 1.0*pow(-my405 + y405, 2)/pow(sigmay405, 3);
            break;
        case 406:
            dJydsigmay[406] = 1.0/sigmay406 - 1.0*pow(-my406 + y406, 2)/pow(sigmay406, 3);
            break;
        case 407:
            dJydsigmay[407] = 1.0/sigmay407 - 1.0*pow(-my407 + y407, 2)/pow(sigmay407, 3);
            break;
        case 408:
            dJydsigmay[408] = 1.0/sigmay408 - 1.0*pow(-my408 + y408, 2)/pow(sigmay408, 3);
            break;
        case 409:
            dJydsigmay[409] = 1.0/sigmay409 - 1.0*pow(-my409 + y409, 2)/pow(sigmay409, 3);
            break;
        case 410:
            dJydsigmay[410] = 1.0/sigmay410 - 1.0*pow(-my410 + y410, 2)/pow(sigmay410, 3);
            break;
        case 411:
            dJydsigmay[411] = 1.0/sigmay411 - 1.0*pow(-my411 + y411, 2)/pow(sigmay411, 3);
            break;
        case 412:
            dJydsigmay[412] = 1.0/sigmay412 - 1.0*pow(-my412 + y412, 2)/pow(sigmay412, 3);
            break;
        case 413:
            dJydsigmay[413] = 1.0/sigmay413 - 1.0*pow(-my413 + y413, 2)/pow(sigmay413, 3);
            break;
        case 414:
            dJydsigmay[414] = 1.0/sigmay414 - 1.0*pow(-my414 + y414, 2)/pow(sigmay414, 3);
            break;
        case 415:
            dJydsigmay[415] = 1.0/sigmay415 - 1.0*pow(-my415 + y415, 2)/pow(sigmay415, 3);
            break;
        case 416:
            dJydsigmay[416] = 1.0/sigmay416 - 1.0*pow(-my416 + y416, 2)/pow(sigmay416, 3);
            break;
        case 417:
            dJydsigmay[417] = 1.0/sigmay417 - 1.0*pow(-my417 + y417, 2)/pow(sigmay417, 3);
            break;
        case 418:
            dJydsigmay[418] = 1.0/sigmay418 - 1.0*pow(-my418 + y418, 2)/pow(sigmay418, 3);
            break;
        case 419:
            dJydsigmay[419] = 1.0/sigmay419 - 1.0*pow(-my419 + y419, 2)/pow(sigmay419, 3);
            break;
        case 420:
            dJydsigmay[420] = 1.0/sigmay420 - 1.0*pow(-my420 + y420, 2)/pow(sigmay420, 3);
            break;
        case 421:
            dJydsigmay[421] = 1.0/sigmay421 - 1.0*pow(-my421 + y421, 2)/pow(sigmay421, 3);
            break;
        case 422:
            dJydsigmay[422] = 1.0/sigmay422 - 1.0*pow(-my422 + y422, 2)/pow(sigmay422, 3);
            break;
        case 423:
            dJydsigmay[423] = 1.0/sigmay423 - 1.0*pow(-my423 + y423, 2)/pow(sigmay423, 3);
            break;
        case 424:
            dJydsigmay[424] = 1.0/sigmay424 - 1.0*pow(-my424 + y424, 2)/pow(sigmay424, 3);
            break;
        case 425:
            dJydsigmay[425] = 1.0/sigmay425 - 1.0*pow(-my425 + y425, 2)/pow(sigmay425, 3);
            break;
        case 426:
            dJydsigmay[426] = 1.0/sigmay426 - 1.0*pow(-my426 + y426, 2)/pow(sigmay426, 3);
            break;
        case 427:
            dJydsigmay[427] = 1.0/sigmay427 - 1.0*pow(-my427 + y427, 2)/pow(sigmay427, 3);
            break;
        case 428:
            dJydsigmay[428] = 1.0/sigmay428 - 1.0*pow(-my428 + y428, 2)/pow(sigmay428, 3);
            break;
        case 429:
            dJydsigmay[429] = 1.0/sigmay429 - 1.0*pow(-my429 + y429, 2)/pow(sigmay429, 3);
            break;
        case 430:
            dJydsigmay[430] = 1.0/sigmay430 - 1.0*pow(-my430 + y430, 2)/pow(sigmay430, 3);
            break;
        case 431:
            dJydsigmay[431] = 1.0/sigmay431 - 1.0*pow(-my431 + y431, 2)/pow(sigmay431, 3);
            break;
        case 432:
            dJydsigmay[432] = 1.0/sigmay432 - 1.0*pow(-my432 + y432, 2)/pow(sigmay432, 3);
            break;
        case 433:
            dJydsigmay[433] = 1.0/sigmay433 - 1.0*pow(-my433 + y433, 2)/pow(sigmay433, 3);
            break;
        case 434:
            dJydsigmay[434] = 1.0/sigmay434 - 1.0*pow(-my434 + y434, 2)/pow(sigmay434, 3);
            break;
        case 435:
            dJydsigmay[435] = 1.0/sigmay435 - 1.0*pow(-my435 + y435, 2)/pow(sigmay435, 3);
            break;
        case 436:
            dJydsigmay[436] = 1.0/sigmay436 - 1.0*pow(-my436 + y436, 2)/pow(sigmay436, 3);
            break;
        case 437:
            dJydsigmay[437] = 1.0/sigmay437 - 1.0*pow(-my437 + y437, 2)/pow(sigmay437, 3);
            break;
        case 438:
            dJydsigmay[438] = 1.0/sigmay438 - 1.0*pow(-my438 + y438, 2)/pow(sigmay438, 3);
            break;
        case 439:
            dJydsigmay[439] = 1.0/sigmay439 - 1.0*pow(-my439 + y439, 2)/pow(sigmay439, 3);
            break;
        case 440:
            dJydsigmay[440] = 1.0/sigmay440 - 1.0*pow(-my440 + y440, 2)/pow(sigmay440, 3);
            break;
        case 441:
            dJydsigmay[441] = 1.0/sigmay441 - 1.0*pow(-my441 + y441, 2)/pow(sigmay441, 3);
            break;
        case 442:
            dJydsigmay[442] = 1.0/sigmay442 - 1.0*pow(-my442 + y442, 2)/pow(sigmay442, 3);
            break;
        case 443:
            dJydsigmay[443] = 1.0/sigmay443 - 1.0*pow(-my443 + y443, 2)/pow(sigmay443, 3);
            break;
        case 444:
            dJydsigmay[444] = 1.0/sigmay444 - 1.0*pow(-my444 + y444, 2)/pow(sigmay444, 3);
            break;
        case 445:
            dJydsigmay[445] = 1.0/sigmay445 - 1.0*pow(-my445 + y445, 2)/pow(sigmay445, 3);
            break;
        case 446:
            dJydsigmay[446] = 1.0/sigmay446 - 1.0*pow(-my446 + y446, 2)/pow(sigmay446, 3);
            break;
        case 447:
            dJydsigmay[447] = 1.0/sigmay447 - 1.0*pow(-my447 + y447, 2)/pow(sigmay447, 3);
            break;
        case 448:
            dJydsigmay[448] = 1.0/sigmay448 - 1.0*pow(-my448 + y448, 2)/pow(sigmay448, 3);
            break;
        case 449:
            dJydsigmay[449] = 1.0/sigmay449 - 1.0*pow(-my449 + y449, 2)/pow(sigmay449, 3);
            break;
        case 450:
            dJydsigmay[450] = 1.0/sigmay450 - 1.0*pow(-my450 + y450, 2)/pow(sigmay450, 3);
            break;
        case 451:
            dJydsigmay[451] = 1.0/sigmay451 - 1.0*pow(-my451 + y451, 2)/pow(sigmay451, 3);
            break;
        case 452:
            dJydsigmay[452] = 1.0/sigmay452 - 1.0*pow(-my452 + y452, 2)/pow(sigmay452, 3);
            break;
        case 453:
            dJydsigmay[453] = 1.0/sigmay453 - 1.0*pow(-my453 + y453, 2)/pow(sigmay453, 3);
            break;
        case 454:
            dJydsigmay[454] = 1.0/sigmay454 - 1.0*pow(-my454 + y454, 2)/pow(sigmay454, 3);
            break;
        case 455:
            dJydsigmay[455] = 1.0/sigmay455 - 1.0*pow(-my455 + y455, 2)/pow(sigmay455, 3);
            break;
        case 456:
            dJydsigmay[456] = 1.0/sigmay456 - 1.0*pow(-my456 + y456, 2)/pow(sigmay456, 3);
            break;
        case 457:
            dJydsigmay[457] = 1.0/sigmay457 - 1.0*pow(-my457 + y457, 2)/pow(sigmay457, 3);
            break;
        case 458:
            dJydsigmay[458] = 1.0/sigmay458 - 1.0*pow(-my458 + y458, 2)/pow(sigmay458, 3);
            break;
        case 459:
            dJydsigmay[459] = 1.0/sigmay459 - 1.0*pow(-my459 + y459, 2)/pow(sigmay459, 3);
            break;
        case 460:
            dJydsigmay[460] = 1.0/sigmay460 - 1.0*pow(-my460 + y460, 2)/pow(sigmay460, 3);
            break;
        case 461:
            dJydsigmay[461] = 1.0/sigmay461 - 1.0*pow(-my461 + y461, 2)/pow(sigmay461, 3);
            break;
        case 462:
            dJydsigmay[462] = 1.0/sigmay462 - 1.0*pow(-my462 + y462, 2)/pow(sigmay462, 3);
            break;
        case 463:
            dJydsigmay[463] = 1.0/sigmay463 - 1.0*pow(-my463 + y463, 2)/pow(sigmay463, 3);
            break;
        case 464:
            dJydsigmay[464] = 1.0/sigmay464 - 1.0*pow(-my464 + y464, 2)/pow(sigmay464, 3);
            break;
        case 465:
            dJydsigmay[465] = 1.0/sigmay465 - 1.0*pow(-my465 + y465, 2)/pow(sigmay465, 3);
            break;
        case 466:
            dJydsigmay[466] = 1.0/sigmay466 - 1.0*pow(-my466 + y466, 2)/pow(sigmay466, 3);
            break;
        case 467:
            dJydsigmay[467] = 1.0/sigmay467 - 1.0*pow(-my467 + y467, 2)/pow(sigmay467, 3);
            break;
        case 468:
            dJydsigmay[468] = 1.0/sigmay468 - 1.0*pow(-my468 + y468, 2)/pow(sigmay468, 3);
            break;
        case 469:
            dJydsigmay[469] = 1.0/sigmay469 - 1.0*pow(-my469 + y469, 2)/pow(sigmay469, 3);
            break;
        case 470:
            dJydsigmay[470] = 1.0/sigmay470 - 1.0*pow(-my470 + y470, 2)/pow(sigmay470, 3);
            break;
        case 471:
            dJydsigmay[471] = 1.0/sigmay471 - 1.0*pow(-my471 + y471, 2)/pow(sigmay471, 3);
            break;
        case 472:
            dJydsigmay[472] = 1.0/sigmay472 - 1.0*pow(-my472 + y472, 2)/pow(sigmay472, 3);
            break;
        case 473:
            dJydsigmay[473] = 1.0/sigmay473 - 1.0*pow(-my473 + y473, 2)/pow(sigmay473, 3);
            break;
        case 474:
            dJydsigmay[474] = 1.0/sigmay474 - 1.0*pow(-my474 + y474, 2)/pow(sigmay474, 3);
            break;
        case 475:
            dJydsigmay[475] = 1.0/sigmay475 - 1.0*pow(-my475 + y475, 2)/pow(sigmay475, 3);
            break;
        case 476:
            dJydsigmay[476] = 1.0/sigmay476 - 1.0*pow(-my476 + y476, 2)/pow(sigmay476, 3);
            break;
        case 477:
            dJydsigmay[477] = 1.0/sigmay477 - 1.0*pow(-my477 + y477, 2)/pow(sigmay477, 3);
            break;
        case 478:
            dJydsigmay[478] = 1.0/sigmay478 - 1.0*pow(-my478 + y478, 2)/pow(sigmay478, 3);
            break;
        case 479:
            dJydsigmay[479] = 1.0/sigmay479 - 1.0*pow(-my479 + y479, 2)/pow(sigmay479, 3);
            break;
        case 480:
            dJydsigmay[480] = 1.0/sigmay480 - 1.0*pow(-my480 + y480, 2)/pow(sigmay480, 3);
            break;
        case 481:
            dJydsigmay[481] = 1.0/sigmay481 - 1.0*pow(-my481 + y481, 2)/pow(sigmay481, 3);
            break;
        case 482:
            dJydsigmay[482] = 1.0/sigmay482 - 1.0*pow(-my482 + y482, 2)/pow(sigmay482, 3);
            break;
        case 483:
            dJydsigmay[483] = 1.0/sigmay483 - 1.0*pow(-my483 + y483, 2)/pow(sigmay483, 3);
            break;
        case 484:
            dJydsigmay[484] = 1.0/sigmay484 - 1.0*pow(-my484 + y484, 2)/pow(sigmay484, 3);
            break;
        case 485:
            dJydsigmay[485] = 1.0/sigmay485 - 1.0*pow(-my485 + y485, 2)/pow(sigmay485, 3);
            break;
        case 486:
            dJydsigmay[486] = 1.0/sigmay486 - 1.0*pow(-my486 + y486, 2)/pow(sigmay486, 3);
            break;
        case 487:
            dJydsigmay[487] = 1.0/sigmay487 - 1.0*pow(-my487 + y487, 2)/pow(sigmay487, 3);
            break;
        case 488:
            dJydsigmay[488] = 1.0/sigmay488 - 1.0*pow(-my488 + y488, 2)/pow(sigmay488, 3);
            break;
        case 489:
            dJydsigmay[489] = 1.0/sigmay489 - 1.0*pow(-my489 + y489, 2)/pow(sigmay489, 3);
            break;
        case 490:
            dJydsigmay[490] = 1.0/sigmay490 - 1.0*pow(-my490 + y490, 2)/pow(sigmay490, 3);
            break;
        case 491:
            dJydsigmay[491] = 1.0/sigmay491 - 1.0*pow(-my491 + y491, 2)/pow(sigmay491, 3);
            break;
        case 492:
            dJydsigmay[492] = 1.0/sigmay492 - 1.0*pow(-my492 + y492, 2)/pow(sigmay492, 3);
            break;
        case 493:
            dJydsigmay[493] = 1.0/sigmay493 - 1.0*pow(-my493 + y493, 2)/pow(sigmay493, 3);
            break;
        case 494:
            dJydsigmay[494] = 1.0/sigmay494 - 1.0*pow(-my494 + y494, 2)/pow(sigmay494, 3);
            break;
        case 495:
            dJydsigmay[495] = 1.0/sigmay495 - 1.0*pow(-my495 + y495, 2)/pow(sigmay495, 3);
            break;
        case 496:
            dJydsigmay[496] = 1.0/sigmay496 - 1.0*pow(-my496 + y496, 2)/pow(sigmay496, 3);
            break;
        case 497:
            dJydsigmay[497] = 1.0/sigmay497 - 1.0*pow(-my497 + y497, 2)/pow(sigmay497, 3);
            break;
        case 498:
            dJydsigmay[498] = 1.0/sigmay498 - 1.0*pow(-my498 + y498, 2)/pow(sigmay498, 3);
            break;
        case 499:
            dJydsigmay[499] = 1.0/sigmay499 - 1.0*pow(-my499 + y499, 2)/pow(sigmay499, 3);
            break;
        case 500:
            dJydsigmay[500] = 1.0/sigmay500 - 1.0*pow(-my500 + y500, 2)/pow(sigmay500, 3);
            break;
        case 501:
            dJydsigmay[501] = 1.0/sigmay501 - 1.0*pow(-my501 + y501, 2)/pow(sigmay501, 3);
            break;
        case 502:
            dJydsigmay[502] = 1.0/sigmay502 - 1.0*pow(-my502 + y502, 2)/pow(sigmay502, 3);
            break;
        case 503:
            dJydsigmay[503] = 1.0/sigmay503 - 1.0*pow(-my503 + y503, 2)/pow(sigmay503, 3);
            break;
        case 504:
            dJydsigmay[504] = 1.0/sigmay504 - 1.0*pow(-my504 + y504, 2)/pow(sigmay504, 3);
            break;
        case 505:
            dJydsigmay[505] = 1.0/sigmay505 - 1.0*pow(-my505 + y505, 2)/pow(sigmay505, 3);
            break;
        case 506:
            dJydsigmay[506] = 1.0/sigmay506 - 1.0*pow(-my506 + y506, 2)/pow(sigmay506, 3);
            break;
        case 507:
            dJydsigmay[507] = 1.0/sigmay507 - 1.0*pow(-my507 + y507, 2)/pow(sigmay507, 3);
            break;
        case 508:
            dJydsigmay[508] = 1.0/sigmay508 - 1.0*pow(-my508 + y508, 2)/pow(sigmay508, 3);
            break;
        case 509:
            dJydsigmay[509] = 1.0/sigmay509 - 1.0*pow(-my509 + y509, 2)/pow(sigmay509, 3);
            break;
        case 510:
            dJydsigmay[510] = 1.0/sigmay510 - 1.0*pow(-my510 + y510, 2)/pow(sigmay510, 3);
            break;
        case 511:
            dJydsigmay[511] = 1.0/sigmay511 - 1.0*pow(-my511 + y511, 2)/pow(sigmay511, 3);
            break;
        case 512:
            dJydsigmay[512] = 1.0/sigmay512 - 1.0*pow(-my512 + y512, 2)/pow(sigmay512, 3);
            break;
        case 513:
            dJydsigmay[513] = 1.0/sigmay513 - 1.0*pow(-my513 + y513, 2)/pow(sigmay513, 3);
            break;
        case 514:
            dJydsigmay[514] = 1.0/sigmay514 - 1.0*pow(-my514 + y514, 2)/pow(sigmay514, 3);
            break;
        case 515:
            dJydsigmay[515] = 1.0/sigmay515 - 1.0*pow(-my515 + y515, 2)/pow(sigmay515, 3);
            break;
        case 516:
            dJydsigmay[516] = 1.0/sigmay516 - 1.0*pow(-my516 + y516, 2)/pow(sigmay516, 3);
            break;
        case 517:
            dJydsigmay[517] = 1.0/sigmay517 - 1.0*pow(-my517 + y517, 2)/pow(sigmay517, 3);
            break;
        case 518:
            dJydsigmay[518] = 1.0/sigmay518 - 1.0*pow(-my518 + y518, 2)/pow(sigmay518, 3);
            break;
        case 519:
            dJydsigmay[519] = 1.0/sigmay519 - 1.0*pow(-my519 + y519, 2)/pow(sigmay519, 3);
            break;
        case 520:
            dJydsigmay[520] = 1.0/sigmay520 - 1.0*pow(-my520 + y520, 2)/pow(sigmay520, 3);
            break;
        case 521:
            dJydsigmay[521] = 1.0/sigmay521 - 1.0*pow(-my521 + y521, 2)/pow(sigmay521, 3);
            break;
        case 522:
            dJydsigmay[522] = 1.0/sigmay522 - 1.0*pow(-my522 + y522, 2)/pow(sigmay522, 3);
            break;
        case 523:
            dJydsigmay[523] = 1.0/sigmay523 - 1.0*pow(-my523 + y523, 2)/pow(sigmay523, 3);
            break;
        case 524:
            dJydsigmay[524] = 1.0/sigmay524 - 1.0*pow(-my524 + y524, 2)/pow(sigmay524, 3);
            break;
        case 525:
            dJydsigmay[525] = 1.0/sigmay525 - 1.0*pow(-my525 + y525, 2)/pow(sigmay525, 3);
            break;
        case 526:
            dJydsigmay[526] = 1.0/sigmay526 - 1.0*pow(-my526 + y526, 2)/pow(sigmay526, 3);
            break;
        case 527:
            dJydsigmay[527] = 1.0/sigmay527 - 1.0*pow(-my527 + y527, 2)/pow(sigmay527, 3);
            break;
        case 528:
            dJydsigmay[528] = 1.0/sigmay528 - 1.0*pow(-my528 + y528, 2)/pow(sigmay528, 3);
            break;
        case 529:
            dJydsigmay[529] = 1.0/sigmay529 - 1.0*pow(-my529 + y529, 2)/pow(sigmay529, 3);
            break;
        case 530:
            dJydsigmay[530] = 1.0/sigmay530 - 1.0*pow(-my530 + y530, 2)/pow(sigmay530, 3);
            break;
        case 531:
            dJydsigmay[531] = 1.0/sigmay531 - 1.0*pow(-my531 + y531, 2)/pow(sigmay531, 3);
            break;
        case 532:
            dJydsigmay[532] = 1.0/sigmay532 - 1.0*pow(-my532 + y532, 2)/pow(sigmay532, 3);
            break;
        case 533:
            dJydsigmay[533] = 1.0/sigmay533 - 1.0*pow(-my533 + y533, 2)/pow(sigmay533, 3);
            break;
        case 534:
            dJydsigmay[534] = 1.0/sigmay534 - 1.0*pow(-my534 + y534, 2)/pow(sigmay534, 3);
            break;
        case 535:
            dJydsigmay[535] = 1.0/sigmay535 - 1.0*pow(-my535 + y535, 2)/pow(sigmay535, 3);
            break;
        case 536:
            dJydsigmay[536] = 1.0/sigmay536 - 1.0*pow(-my536 + y536, 2)/pow(sigmay536, 3);
            break;
        case 537:
            dJydsigmay[537] = 1.0/sigmay537 - 1.0*pow(-my537 + y537, 2)/pow(sigmay537, 3);
            break;
        case 538:
            dJydsigmay[538] = 1.0/sigmay538 - 1.0*pow(-my538 + y538, 2)/pow(sigmay538, 3);
            break;
        case 539:
            dJydsigmay[539] = 1.0/sigmay539 - 1.0*pow(-my539 + y539, 2)/pow(sigmay539, 3);
            break;
        case 540:
            dJydsigmay[540] = 1.0/sigmay540 - 1.0*pow(-my540 + y540, 2)/pow(sigmay540, 3);
            break;
        case 541:
            dJydsigmay[541] = 1.0/sigmay541 - 1.0*pow(-my541 + y541, 2)/pow(sigmay541, 3);
            break;
        case 542:
            dJydsigmay[542] = 1.0/sigmay542 - 1.0*pow(-my542 + y542, 2)/pow(sigmay542, 3);
            break;
        case 543:
            dJydsigmay[543] = 1.0/sigmay543 - 1.0*pow(-my543 + y543, 2)/pow(sigmay543, 3);
            break;
        case 544:
            dJydsigmay[544] = 1.0/sigmay544 - 1.0*pow(-my544 + y544, 2)/pow(sigmay544, 3);
            break;
        case 545:
            dJydsigmay[545] = 1.0/sigmay545 - 1.0*pow(-my545 + y545, 2)/pow(sigmay545, 3);
            break;
        case 546:
            dJydsigmay[546] = 1.0/sigmay546 - 1.0*pow(-my546 + y546, 2)/pow(sigmay546, 3);
            break;
        case 547:
            dJydsigmay[547] = 1.0/sigmay547 - 1.0*pow(-my547 + y547, 2)/pow(sigmay547, 3);
            break;
        case 548:
            dJydsigmay[548] = 1.0/sigmay548 - 1.0*pow(-my548 + y548, 2)/pow(sigmay548, 3);
            break;
        case 549:
            dJydsigmay[549] = 1.0/sigmay549 - 1.0*pow(-my549 + y549, 2)/pow(sigmay549, 3);
            break;
        case 550:
            dJydsigmay[550] = 1.0/sigmay550 - 1.0*pow(-my550 + y550, 2)/pow(sigmay550, 3);
            break;
        case 551:
            dJydsigmay[551] = 1.0/sigmay551 - 1.0*pow(-my551 + y551, 2)/pow(sigmay551, 3);
            break;
        case 552:
            dJydsigmay[552] = 1.0/sigmay552 - 1.0*pow(-my552 + y552, 2)/pow(sigmay552, 3);
            break;
        case 553:
            dJydsigmay[553] = 1.0/sigmay553 - 1.0*pow(-my553 + y553, 2)/pow(sigmay553, 3);
            break;
        case 554:
            dJydsigmay[554] = 1.0/sigmay554 - 1.0*pow(-my554 + y554, 2)/pow(sigmay554, 3);
            break;
        case 555:
            dJydsigmay[555] = 1.0/sigmay555 - 1.0*pow(-my555 + y555, 2)/pow(sigmay555, 3);
            break;
        case 556:
            dJydsigmay[556] = 1.0/sigmay556 - 1.0*pow(-my556 + y556, 2)/pow(sigmay556, 3);
            break;
        case 557:
            dJydsigmay[557] = 1.0/sigmay557 - 1.0*pow(-my557 + y557, 2)/pow(sigmay557, 3);
            break;
        case 558:
            dJydsigmay[558] = 1.0/sigmay558 - 1.0*pow(-my558 + y558, 2)/pow(sigmay558, 3);
            break;
        case 559:
            dJydsigmay[559] = 1.0/sigmay559 - 1.0*pow(-my559 + y559, 2)/pow(sigmay559, 3);
            break;
        case 560:
            dJydsigmay[560] = 1.0/sigmay560 - 1.0*pow(-my560 + y560, 2)/pow(sigmay560, 3);
            break;
        case 561:
            dJydsigmay[561] = 1.0/sigmay561 - 1.0*pow(-my561 + y561, 2)/pow(sigmay561, 3);
            break;
        case 562:
            dJydsigmay[562] = 1.0/sigmay562 - 1.0*pow(-my562 + y562, 2)/pow(sigmay562, 3);
            break;
        case 563:
            dJydsigmay[563] = 1.0/sigmay563 - 1.0*pow(-my563 + y563, 2)/pow(sigmay563, 3);
            break;
        case 564:
            dJydsigmay[564] = 1.0/sigmay564 - 1.0*pow(-my564 + y564, 2)/pow(sigmay564, 3);
            break;
        case 565:
            dJydsigmay[565] = 1.0/sigmay565 - 1.0*pow(-my565 + y565, 2)/pow(sigmay565, 3);
            break;
        case 566:
            dJydsigmay[566] = 1.0/sigmay566 - 1.0*pow(-my566 + y566, 2)/pow(sigmay566, 3);
            break;
        case 567:
            dJydsigmay[567] = 1.0/sigmay567 - 1.0*pow(-my567 + y567, 2)/pow(sigmay567, 3);
            break;
        case 568:
            dJydsigmay[568] = 1.0/sigmay568 - 1.0*pow(-my568 + y568, 2)/pow(sigmay568, 3);
            break;
        case 569:
            dJydsigmay[569] = 1.0/sigmay569 - 1.0*pow(-my569 + y569, 2)/pow(sigmay569, 3);
            break;
        case 570:
            dJydsigmay[570] = 1.0/sigmay570 - 1.0*pow(-my570 + y570, 2)/pow(sigmay570, 3);
            break;
        case 571:
            dJydsigmay[571] = 1.0/sigmay571 - 1.0*pow(-my571 + y571, 2)/pow(sigmay571, 3);
            break;
        case 572:
            dJydsigmay[572] = 1.0/sigmay572 - 1.0*pow(-my572 + y572, 2)/pow(sigmay572, 3);
            break;
        case 573:
            dJydsigmay[573] = 1.0/sigmay573 - 1.0*pow(-my573 + y573, 2)/pow(sigmay573, 3);
            break;
        case 574:
            dJydsigmay[574] = 1.0/sigmay574 - 1.0*pow(-my574 + y574, 2)/pow(sigmay574, 3);
            break;
        case 575:
            dJydsigmay[575] = 1.0/sigmay575 - 1.0*pow(-my575 + y575, 2)/pow(sigmay575, 3);
            break;
        case 576:
            dJydsigmay[576] = 1.0/sigmay576 - 1.0*pow(-my576 + y576, 2)/pow(sigmay576, 3);
            break;
        case 577:
            dJydsigmay[577] = 1.0/sigmay577 - 1.0*pow(-my577 + y577, 2)/pow(sigmay577, 3);
            break;
        case 578:
            dJydsigmay[578] = 1.0/sigmay578 - 1.0*pow(-my578 + y578, 2)/pow(sigmay578, 3);
            break;
        case 579:
            dJydsigmay[579] = 1.0/sigmay579 - 1.0*pow(-my579 + y579, 2)/pow(sigmay579, 3);
            break;
        case 580:
            dJydsigmay[580] = 1.0/sigmay580 - 1.0*pow(-my580 + y580, 2)/pow(sigmay580, 3);
            break;
        case 581:
            dJydsigmay[581] = 1.0/sigmay581 - 1.0*pow(-my581 + y581, 2)/pow(sigmay581, 3);
            break;
        case 582:
            dJydsigmay[582] = 1.0/sigmay582 - 1.0*pow(-my582 + y582, 2)/pow(sigmay582, 3);
            break;
        case 583:
            dJydsigmay[583] = 1.0/sigmay583 - 1.0*pow(-my583 + y583, 2)/pow(sigmay583, 3);
            break;
        case 584:
            dJydsigmay[584] = 1.0/sigmay584 - 1.0*pow(-my584 + y584, 2)/pow(sigmay584, 3);
            break;
        case 585:
            dJydsigmay[585] = 1.0/sigmay585 - 1.0*pow(-my585 + y585, 2)/pow(sigmay585, 3);
            break;
        case 586:
            dJydsigmay[586] = 1.0/sigmay586 - 1.0*pow(-my586 + y586, 2)/pow(sigmay586, 3);
            break;
        case 587:
            dJydsigmay[587] = 1.0/sigmay587 - 1.0*pow(-my587 + y587, 2)/pow(sigmay587, 3);
            break;
        case 588:
            dJydsigmay[588] = 1.0/sigmay588 - 1.0*pow(-my588 + y588, 2)/pow(sigmay588, 3);
            break;
        case 589:
            dJydsigmay[589] = 1.0/sigmay589 - 1.0*pow(-my589 + y589, 2)/pow(sigmay589, 3);
            break;
        case 590:
            dJydsigmay[590] = 1.0/sigmay590 - 1.0*pow(-my590 + y590, 2)/pow(sigmay590, 3);
            break;
        case 591:
            dJydsigmay[591] = 1.0/sigmay591 - 1.0*pow(-my591 + y591, 2)/pow(sigmay591, 3);
            break;
        case 592:
            dJydsigmay[592] = 1.0/sigmay592 - 1.0*pow(-my592 + y592, 2)/pow(sigmay592, 3);
            break;
        case 593:
            dJydsigmay[593] = 1.0/sigmay593 - 1.0*pow(-my593 + y593, 2)/pow(sigmay593, 3);
            break;
        case 594:
            dJydsigmay[594] = 1.0/sigmay594 - 1.0*pow(-my594 + y594, 2)/pow(sigmay594, 3);
            break;
        case 595:
            dJydsigmay[595] = 1.0/sigmay595 - 1.0*pow(-my595 + y595, 2)/pow(sigmay595, 3);
            break;
        case 596:
            dJydsigmay[596] = 1.0/sigmay596 - 1.0*pow(-my596 + y596, 2)/pow(sigmay596, 3);
            break;
        case 597:
            dJydsigmay[597] = 1.0/sigmay597 - 1.0*pow(-my597 + y597, 2)/pow(sigmay597, 3);
            break;
        case 598:
            dJydsigmay[598] = 1.0/sigmay598 - 1.0*pow(-my598 + y598, 2)/pow(sigmay598, 3);
            break;
        case 599:
            dJydsigmay[599] = 1.0/sigmay599 - 1.0*pow(-my599 + y599, 2)/pow(sigmay599, 3);
            break;
        case 600:
            dJydsigmay[600] = 1.0/sigmay600 - 1.0*pow(-my600 + y600, 2)/pow(sigmay600, 3);
            break;
        case 601:
            dJydsigmay[601] = 1.0/sigmay601 - 1.0*pow(-my601 + y601, 2)/pow(sigmay601, 3);
            break;
        case 602:
            dJydsigmay[602] = 1.0/sigmay602 - 1.0*pow(-my602 + y602, 2)/pow(sigmay602, 3);
            break;
        case 603:
            dJydsigmay[603] = 1.0/sigmay603 - 1.0*pow(-my603 + y603, 2)/pow(sigmay603, 3);
            break;
        case 604:
            dJydsigmay[604] = 1.0/sigmay604 - 1.0*pow(-my604 + y604, 2)/pow(sigmay604, 3);
            break;
        case 605:
            dJydsigmay[605] = 1.0/sigmay605 - 1.0*pow(-my605 + y605, 2)/pow(sigmay605, 3);
            break;
        case 606:
            dJydsigmay[606] = 1.0/sigmay606 - 1.0*pow(-my606 + y606, 2)/pow(sigmay606, 3);
            break;
        case 607:
            dJydsigmay[607] = 1.0/sigmay607 - 1.0*pow(-my607 + y607, 2)/pow(sigmay607, 3);
            break;
        case 608:
            dJydsigmay[608] = 1.0/sigmay608 - 1.0*pow(-my608 + y608, 2)/pow(sigmay608, 3);
            break;
        case 609:
            dJydsigmay[609] = 1.0/sigmay609 - 1.0*pow(-my609 + y609, 2)/pow(sigmay609, 3);
            break;
        case 610:
            dJydsigmay[610] = 1.0/sigmay610 - 1.0*pow(-my610 + y610, 2)/pow(sigmay610, 3);
            break;
        case 611:
            dJydsigmay[611] = 1.0/sigmay611 - 1.0*pow(-my611 + y611, 2)/pow(sigmay611, 3);
            break;
        case 612:
            dJydsigmay[612] = 1.0/sigmay612 - 1.0*pow(-my612 + y612, 2)/pow(sigmay612, 3);
            break;
        case 613:
            dJydsigmay[613] = 1.0/sigmay613 - 1.0*pow(-my613 + y613, 2)/pow(sigmay613, 3);
            break;
        case 614:
            dJydsigmay[614] = 1.0/sigmay614 - 1.0*pow(-my614 + y614, 2)/pow(sigmay614, 3);
            break;
        case 615:
            dJydsigmay[615] = 1.0/sigmay615 - 1.0*pow(-my615 + y615, 2)/pow(sigmay615, 3);
            break;
        case 616:
            dJydsigmay[616] = 1.0/sigmay616 - 1.0*pow(-my616 + y616, 2)/pow(sigmay616, 3);
            break;
        case 617:
            dJydsigmay[617] = 1.0/sigmay617 - 1.0*pow(-my617 + y617, 2)/pow(sigmay617, 3);
            break;
        case 618:
            dJydsigmay[618] = 1.0/sigmay618 - 1.0*pow(-my618 + y618, 2)/pow(sigmay618, 3);
            break;
        case 619:
            dJydsigmay[619] = 1.0/sigmay619 - 1.0*pow(-my619 + y619, 2)/pow(sigmay619, 3);
            break;
        case 620:
            dJydsigmay[620] = 1.0/sigmay620 - 1.0*pow(-my620 + y620, 2)/pow(sigmay620, 3);
            break;
        case 621:
            dJydsigmay[621] = 1.0/sigmay621 - 1.0*pow(-my621 + y621, 2)/pow(sigmay621, 3);
            break;
        case 622:
            dJydsigmay[622] = 1.0/sigmay622 - 1.0*pow(-my622 + y622, 2)/pow(sigmay622, 3);
            break;
        case 623:
            dJydsigmay[623] = 1.0/sigmay623 - 1.0*pow(-my623 + y623, 2)/pow(sigmay623, 3);
            break;
        case 624:
            dJydsigmay[624] = 1.0/sigmay624 - 1.0*pow(-my624 + y624, 2)/pow(sigmay624, 3);
            break;
        case 625:
            dJydsigmay[625] = 1.0/sigmay625 - 1.0*pow(-my625 + y625, 2)/pow(sigmay625, 3);
            break;
        case 626:
            dJydsigmay[626] = 1.0/sigmay626 - 1.0*pow(-my626 + y626, 2)/pow(sigmay626, 3);
            break;
        case 627:
            dJydsigmay[627] = 1.0/sigmay627 - 1.0*pow(-my627 + y627, 2)/pow(sigmay627, 3);
            break;
        case 628:
            dJydsigmay[628] = 1.0/sigmay628 - 1.0*pow(-my628 + y628, 2)/pow(sigmay628, 3);
            break;
        case 629:
            dJydsigmay[629] = 1.0/sigmay629 - 1.0*pow(-my629 + y629, 2)/pow(sigmay629, 3);
            break;
        case 630:
            dJydsigmay[630] = 1.0/sigmay630 - 1.0*pow(-my630 + y630, 2)/pow(sigmay630, 3);
            break;
        case 631:
            dJydsigmay[631] = 1.0/sigmay631 - 1.0*pow(-my631 + y631, 2)/pow(sigmay631, 3);
            break;
        case 632:
            dJydsigmay[632] = 1.0/sigmay632 - 1.0*pow(-my632 + y632, 2)/pow(sigmay632, 3);
            break;
        case 633:
            dJydsigmay[633] = 1.0/sigmay633 - 1.0*pow(-my633 + y633, 2)/pow(sigmay633, 3);
            break;
        case 634:
            dJydsigmay[634] = 1.0/sigmay634 - 1.0*pow(-my634 + y634, 2)/pow(sigmay634, 3);
            break;
        case 635:
            dJydsigmay[635] = 1.0/sigmay635 - 1.0*pow(-my635 + y635, 2)/pow(sigmay635, 3);
            break;
        case 636:
            dJydsigmay[636] = 1.0/sigmay636 - 1.0*pow(-my636 + y636, 2)/pow(sigmay636, 3);
            break;
        case 637:
            dJydsigmay[637] = 1.0/sigmay637 - 1.0*pow(-my637 + y637, 2)/pow(sigmay637, 3);
            break;
        case 638:
            dJydsigmay[638] = 1.0/sigmay638 - 1.0*pow(-my638 + y638, 2)/pow(sigmay638, 3);
            break;
        case 639:
            dJydsigmay[639] = 1.0/sigmay639 - 1.0*pow(-my639 + y639, 2)/pow(sigmay639, 3);
            break;
        case 640:
            dJydsigmay[640] = 1.0/sigmay640 - 1.0*pow(-my640 + y640, 2)/pow(sigmay640, 3);
            break;
        case 641:
            dJydsigmay[641] = 1.0/sigmay641 - 1.0*pow(-my641 + y641, 2)/pow(sigmay641, 3);
            break;
        case 642:
            dJydsigmay[642] = 1.0/sigmay642 - 1.0*pow(-my642 + y642, 2)/pow(sigmay642, 3);
            break;
        case 643:
            dJydsigmay[643] = 1.0/sigmay643 - 1.0*pow(-my643 + y643, 2)/pow(sigmay643, 3);
            break;
        case 644:
            dJydsigmay[644] = 1.0/sigmay644 - 1.0*pow(-my644 + y644, 2)/pow(sigmay644, 3);
            break;
        case 645:
            dJydsigmay[645] = 1.0/sigmay645 - 1.0*pow(-my645 + y645, 2)/pow(sigmay645, 3);
            break;
        case 646:
            dJydsigmay[646] = 1.0/sigmay646 - 1.0*pow(-my646 + y646, 2)/pow(sigmay646, 3);
            break;
        case 647:
            dJydsigmay[647] = 1.0/sigmay647 - 1.0*pow(-my647 + y647, 2)/pow(sigmay647, 3);
            break;
        case 648:
            dJydsigmay[648] = 1.0/sigmay648 - 1.0*pow(-my648 + y648, 2)/pow(sigmay648, 3);
            break;
        case 649:
            dJydsigmay[649] = 1.0/sigmay649 - 1.0*pow(-my649 + y649, 2)/pow(sigmay649, 3);
            break;
        case 650:
            dJydsigmay[650] = 1.0/sigmay650 - 1.0*pow(-my650 + y650, 2)/pow(sigmay650, 3);
            break;
        case 651:
            dJydsigmay[651] = 1.0/sigmay651 - 1.0*pow(-my651 + y651, 2)/pow(sigmay651, 3);
            break;
        case 652:
            dJydsigmay[652] = 1.0/sigmay652 - 1.0*pow(-my652 + y652, 2)/pow(sigmay652, 3);
            break;
        case 653:
            dJydsigmay[653] = 1.0/sigmay653 - 1.0*pow(-my653 + y653, 2)/pow(sigmay653, 3);
            break;
        case 654:
            dJydsigmay[654] = 1.0/sigmay654 - 1.0*pow(-my654 + y654, 2)/pow(sigmay654, 3);
            break;
        case 655:
            dJydsigmay[655] = 1.0/sigmay655 - 1.0*pow(-my655 + y655, 2)/pow(sigmay655, 3);
            break;
        case 656:
            dJydsigmay[656] = 1.0/sigmay656 - 1.0*pow(-my656 + y656, 2)/pow(sigmay656, 3);
            break;
        case 657:
            dJydsigmay[657] = 1.0/sigmay657 - 1.0*pow(-my657 + y657, 2)/pow(sigmay657, 3);
            break;
        case 658:
            dJydsigmay[658] = 1.0/sigmay658 - 1.0*pow(-my658 + y658, 2)/pow(sigmay658, 3);
            break;
        case 659:
            dJydsigmay[659] = 1.0/sigmay659 - 1.0*pow(-my659 + y659, 2)/pow(sigmay659, 3);
            break;
        case 660:
            dJydsigmay[660] = 1.0/sigmay660 - 1.0*pow(-my660 + y660, 2)/pow(sigmay660, 3);
            break;
        case 661:
            dJydsigmay[661] = 1.0/sigmay661 - 1.0*pow(-my661 + y661, 2)/pow(sigmay661, 3);
            break;
        case 662:
            dJydsigmay[662] = 1.0/sigmay662 - 1.0*pow(-my662 + y662, 2)/pow(sigmay662, 3);
            break;
        case 663:
            dJydsigmay[663] = 1.0/sigmay663 - 1.0*pow(-my663 + y663, 2)/pow(sigmay663, 3);
            break;
        case 664:
            dJydsigmay[664] = 1.0/sigmay664 - 1.0*pow(-my664 + y664, 2)/pow(sigmay664, 3);
            break;
        case 665:
            dJydsigmay[665] = 1.0/sigmay665 - 1.0*pow(-my665 + y665, 2)/pow(sigmay665, 3);
            break;
        case 666:
            dJydsigmay[666] = 1.0/sigmay666 - 1.0*pow(-my666 + y666, 2)/pow(sigmay666, 3);
            break;
        case 667:
            dJydsigmay[667] = 1.0/sigmay667 - 1.0*pow(-my667 + y667, 2)/pow(sigmay667, 3);
            break;
        case 668:
            dJydsigmay[668] = 1.0/sigmay668 - 1.0*pow(-my668 + y668, 2)/pow(sigmay668, 3);
            break;
        case 669:
            dJydsigmay[669] = 1.0/sigmay669 - 1.0*pow(-my669 + y669, 2)/pow(sigmay669, 3);
            break;
        case 670:
            dJydsigmay[670] = 1.0/sigmay670 - 1.0*pow(-my670 + y670, 2)/pow(sigmay670, 3);
            break;
        case 671:
            dJydsigmay[671] = 1.0/sigmay671 - 1.0*pow(-my671 + y671, 2)/pow(sigmay671, 3);
            break;
        case 672:
            dJydsigmay[672] = 1.0/sigmay672 - 1.0*pow(-my672 + y672, 2)/pow(sigmay672, 3);
            break;
        case 673:
            dJydsigmay[673] = 1.0/sigmay673 - 1.0*pow(-my673 + y673, 2)/pow(sigmay673, 3);
            break;
        case 674:
            dJydsigmay[674] = 1.0/sigmay674 - 1.0*pow(-my674 + y674, 2)/pow(sigmay674, 3);
            break;
        case 675:
            dJydsigmay[675] = 1.0/sigmay675 - 1.0*pow(-my675 + y675, 2)/pow(sigmay675, 3);
            break;
        case 676:
            dJydsigmay[676] = 1.0/sigmay676 - 1.0*pow(-my676 + y676, 2)/pow(sigmay676, 3);
            break;
        case 677:
            dJydsigmay[677] = 1.0/sigmay677 - 1.0*pow(-my677 + y677, 2)/pow(sigmay677, 3);
            break;
        case 678:
            dJydsigmay[678] = 1.0/sigmay678 - 1.0*pow(-my678 + y678, 2)/pow(sigmay678, 3);
            break;
        case 679:
            dJydsigmay[679] = 1.0/sigmay679 - 1.0*pow(-my679 + y679, 2)/pow(sigmay679, 3);
            break;
        case 680:
            dJydsigmay[680] = 1.0/sigmay680 - 1.0*pow(-my680 + y680, 2)/pow(sigmay680, 3);
            break;
        case 681:
            dJydsigmay[681] = 1.0/sigmay681 - 1.0*pow(-my681 + y681, 2)/pow(sigmay681, 3);
            break;
        case 682:
            dJydsigmay[682] = 1.0/sigmay682 - 1.0*pow(-my682 + y682, 2)/pow(sigmay682, 3);
            break;
        case 683:
            dJydsigmay[683] = 1.0/sigmay683 - 1.0*pow(-my683 + y683, 2)/pow(sigmay683, 3);
            break;
        case 684:
            dJydsigmay[684] = 1.0/sigmay684 - 1.0*pow(-my684 + y684, 2)/pow(sigmay684, 3);
            break;
        case 685:
            dJydsigmay[685] = 1.0/sigmay685 - 1.0*pow(-my685 + y685, 2)/pow(sigmay685, 3);
            break;
        case 686:
            dJydsigmay[686] = 1.0/sigmay686 - 1.0*pow(-my686 + y686, 2)/pow(sigmay686, 3);
            break;
        case 687:
            dJydsigmay[687] = 1.0/sigmay687 - 1.0*pow(-my687 + y687, 2)/pow(sigmay687, 3);
            break;
        case 688:
            dJydsigmay[688] = 1.0/sigmay688 - 1.0*pow(-my688 + y688, 2)/pow(sigmay688, 3);
            break;
        case 689:
            dJydsigmay[689] = 1.0/sigmay689 - 1.0*pow(-my689 + y689, 2)/pow(sigmay689, 3);
            break;
        case 690:
            dJydsigmay[690] = 1.0/sigmay690 - 1.0*pow(-my690 + y690, 2)/pow(sigmay690, 3);
            break;
        case 691:
            dJydsigmay[691] = 1.0/sigmay691 - 1.0*pow(-my691 + y691, 2)/pow(sigmay691, 3);
            break;
        case 692:
            dJydsigmay[692] = 1.0/sigmay692 - 1.0*pow(-my692 + y692, 2)/pow(sigmay692, 3);
            break;
        case 693:
            dJydsigmay[693] = 1.0/sigmay693 - 1.0*pow(-my693 + y693, 2)/pow(sigmay693, 3);
            break;
        case 694:
            dJydsigmay[694] = 1.0/sigmay694 - 1.0*pow(-my694 + y694, 2)/pow(sigmay694, 3);
            break;
        case 695:
            dJydsigmay[695] = 1.0/sigmay695 - 1.0*pow(-my695 + y695, 2)/pow(sigmay695, 3);
            break;
        case 696:
            dJydsigmay[696] = 1.0/sigmay696 - 1.0*pow(-my696 + y696, 2)/pow(sigmay696, 3);
            break;
        case 697:
            dJydsigmay[697] = 1.0/sigmay697 - 1.0*pow(-my697 + y697, 2)/pow(sigmay697, 3);
            break;
        case 698:
            dJydsigmay[698] = 1.0/sigmay698 - 1.0*pow(-my698 + y698, 2)/pow(sigmay698, 3);
            break;
        case 699:
            dJydsigmay[699] = 1.0/sigmay699 - 1.0*pow(-my699 + y699, 2)/pow(sigmay699, 3);
            break;
        case 700:
            dJydsigmay[700] = 1.0/sigmay700 - 1.0*pow(-my700 + y700, 2)/pow(sigmay700, 3);
            break;
        case 701:
            dJydsigmay[701] = 1.0/sigmay701 - 1.0*pow(-my701 + y701, 2)/pow(sigmay701, 3);
            break;
        case 702:
            dJydsigmay[702] = 1.0/sigmay702 - 1.0*pow(-my702 + y702, 2)/pow(sigmay702, 3);
            break;
        case 703:
            dJydsigmay[703] = 1.0/sigmay703 - 1.0*pow(-my703 + y703, 2)/pow(sigmay703, 3);
            break;
        case 704:
            dJydsigmay[704] = 1.0/sigmay704 - 1.0*pow(-my704 + y704, 2)/pow(sigmay704, 3);
            break;
        case 705:
            dJydsigmay[705] = 1.0/sigmay705 - 1.0*pow(-my705 + y705, 2)/pow(sigmay705, 3);
            break;
        case 706:
            dJydsigmay[706] = 1.0/sigmay706 - 1.0*pow(-my706 + y706, 2)/pow(sigmay706, 3);
            break;
        case 707:
            dJydsigmay[707] = 1.0/sigmay707 - 1.0*pow(-my707 + y707, 2)/pow(sigmay707, 3);
            break;
        case 708:
            dJydsigmay[708] = 1.0/sigmay708 - 1.0*pow(-my708 + y708, 2)/pow(sigmay708, 3);
            break;
        case 709:
            dJydsigmay[709] = 1.0/sigmay709 - 1.0*pow(-my709 + y709, 2)/pow(sigmay709, 3);
            break;
        case 710:
            dJydsigmay[710] = 1.0/sigmay710 - 1.0*pow(-my710 + y710, 2)/pow(sigmay710, 3);
            break;
        case 711:
            dJydsigmay[711] = 1.0/sigmay711 - 1.0*pow(-my711 + y711, 2)/pow(sigmay711, 3);
            break;
        case 712:
            dJydsigmay[712] = 1.0/sigmay712 - 1.0*pow(-my712 + y712, 2)/pow(sigmay712, 3);
            break;
        case 713:
            dJydsigmay[713] = 1.0/sigmay713 - 1.0*pow(-my713 + y713, 2)/pow(sigmay713, 3);
            break;
        case 714:
            dJydsigmay[714] = 1.0/sigmay714 - 1.0*pow(-my714 + y714, 2)/pow(sigmay714, 3);
            break;
        case 715:
            dJydsigmay[715] = 1.0/sigmay715 - 1.0*pow(-my715 + y715, 2)/pow(sigmay715, 3);
            break;
        case 716:
            dJydsigmay[716] = 1.0/sigmay716 - 1.0*pow(-my716 + y716, 2)/pow(sigmay716, 3);
            break;
        case 717:
            dJydsigmay[717] = 1.0/sigmay717 - 1.0*pow(-my717 + y717, 2)/pow(sigmay717, 3);
            break;
        case 718:
            dJydsigmay[718] = 1.0/sigmay718 - 1.0*pow(-my718 + y718, 2)/pow(sigmay718, 3);
            break;
        case 719:
            dJydsigmay[719] = 1.0/sigmay719 - 1.0*pow(-my719 + y719, 2)/pow(sigmay719, 3);
            break;
        case 720:
            dJydsigmay[720] = 1.0/sigmay720 - 1.0*pow(-my720 + y720, 2)/pow(sigmay720, 3);
            break;
        case 721:
            dJydsigmay[721] = 1.0/sigmay721 - 1.0*pow(-my721 + y721, 2)/pow(sigmay721, 3);
            break;
        case 722:
            dJydsigmay[722] = 1.0/sigmay722 - 1.0*pow(-my722 + y722, 2)/pow(sigmay722, 3);
            break;
        case 723:
            dJydsigmay[723] = 1.0/sigmay723 - 1.0*pow(-my723 + y723, 2)/pow(sigmay723, 3);
            break;
        case 724:
            dJydsigmay[724] = 1.0/sigmay724 - 1.0*pow(-my724 + y724, 2)/pow(sigmay724, 3);
            break;
        case 725:
            dJydsigmay[725] = 1.0/sigmay725 - 1.0*pow(-my725 + y725, 2)/pow(sigmay725, 3);
            break;
        case 726:
            dJydsigmay[726] = 1.0/sigmay726 - 1.0*pow(-my726 + y726, 2)/pow(sigmay726, 3);
            break;
        case 727:
            dJydsigmay[727] = 1.0/sigmay727 - 1.0*pow(-my727 + y727, 2)/pow(sigmay727, 3);
            break;
        case 728:
            dJydsigmay[728] = 1.0/sigmay728 - 1.0*pow(-my728 + y728, 2)/pow(sigmay728, 3);
            break;
        case 729:
            dJydsigmay[729] = 1.0/sigmay729 - 1.0*pow(-my729 + y729, 2)/pow(sigmay729, 3);
            break;
        case 730:
            dJydsigmay[730] = 1.0/sigmay730 - 1.0*pow(-my730 + y730, 2)/pow(sigmay730, 3);
            break;
        case 731:
            dJydsigmay[731] = 1.0/sigmay731 - 1.0*pow(-my731 + y731, 2)/pow(sigmay731, 3);
            break;
        case 732:
            dJydsigmay[732] = 1.0/sigmay732 - 1.0*pow(-my732 + y732, 2)/pow(sigmay732, 3);
            break;
        case 733:
            dJydsigmay[733] = 1.0/sigmay733 - 1.0*pow(-my733 + y733, 2)/pow(sigmay733, 3);
            break;
        case 734:
            dJydsigmay[734] = 1.0/sigmay734 - 1.0*pow(-my734 + y734, 2)/pow(sigmay734, 3);
            break;
        case 735:
            dJydsigmay[735] = 1.0/sigmay735 - 1.0*pow(-my735 + y735, 2)/pow(sigmay735, 3);
            break;
        case 736:
            dJydsigmay[736] = 1.0/sigmay736 - 1.0*pow(-my736 + y736, 2)/pow(sigmay736, 3);
            break;
        case 737:
            dJydsigmay[737] = 1.0/sigmay737 - 1.0*pow(-my737 + y737, 2)/pow(sigmay737, 3);
            break;
        case 738:
            dJydsigmay[738] = 1.0/sigmay738 - 1.0*pow(-my738 + y738, 2)/pow(sigmay738, 3);
            break;
        case 739:
            dJydsigmay[739] = 1.0/sigmay739 - 1.0*pow(-my739 + y739, 2)/pow(sigmay739, 3);
            break;
        case 740:
            dJydsigmay[740] = 1.0/sigmay740 - 1.0*pow(-my740 + y740, 2)/pow(sigmay740, 3);
            break;
        case 741:
            dJydsigmay[741] = 1.0/sigmay741 - 1.0*pow(-my741 + y741, 2)/pow(sigmay741, 3);
            break;
        case 742:
            dJydsigmay[742] = 1.0/sigmay742 - 1.0*pow(-my742 + y742, 2)/pow(sigmay742, 3);
            break;
        case 743:
            dJydsigmay[743] = 1.0/sigmay743 - 1.0*pow(-my743 + y743, 2)/pow(sigmay743, 3);
            break;
        case 744:
            dJydsigmay[744] = 1.0/sigmay744 - 1.0*pow(-my744 + y744, 2)/pow(sigmay744, 3);
            break;
        case 745:
            dJydsigmay[745] = 1.0/sigmay745 - 1.0*pow(-my745 + y745, 2)/pow(sigmay745, 3);
            break;
        case 746:
            dJydsigmay[746] = 1.0/sigmay746 - 1.0*pow(-my746 + y746, 2)/pow(sigmay746, 3);
            break;
        case 747:
            dJydsigmay[747] = 1.0/sigmay747 - 1.0*pow(-my747 + y747, 2)/pow(sigmay747, 3);
            break;
        case 748:
            dJydsigmay[748] = 1.0/sigmay748 - 1.0*pow(-my748 + y748, 2)/pow(sigmay748, 3);
            break;
        case 749:
            dJydsigmay[749] = 1.0/sigmay749 - 1.0*pow(-my749 + y749, 2)/pow(sigmay749, 3);
            break;
        case 750:
            dJydsigmay[750] = 1.0/sigmay750 - 1.0*pow(-my750 + y750, 2)/pow(sigmay750, 3);
            break;
        case 751:
            dJydsigmay[751] = 1.0/sigmay751 - 1.0*pow(-my751 + y751, 2)/pow(sigmay751, 3);
            break;
        case 752:
            dJydsigmay[752] = 1.0/sigmay752 - 1.0*pow(-my752 + y752, 2)/pow(sigmay752, 3);
            break;
        case 753:
            dJydsigmay[753] = 1.0/sigmay753 - 1.0*pow(-my753 + y753, 2)/pow(sigmay753, 3);
            break;
        case 754:
            dJydsigmay[754] = 1.0/sigmay754 - 1.0*pow(-my754 + y754, 2)/pow(sigmay754, 3);
            break;
        case 755:
            dJydsigmay[755] = 1.0/sigmay755 - 1.0*pow(-my755 + y755, 2)/pow(sigmay755, 3);
            break;
        case 756:
            dJydsigmay[756] = 1.0/sigmay756 - 1.0*pow(-my756 + y756, 2)/pow(sigmay756, 3);
            break;
        case 757:
            dJydsigmay[757] = 1.0/sigmay757 - 1.0*pow(-my757 + y757, 2)/pow(sigmay757, 3);
            break;
        case 758:
            dJydsigmay[758] = 1.0/sigmay758 - 1.0*pow(-my758 + y758, 2)/pow(sigmay758, 3);
            break;
        case 759:
            dJydsigmay[759] = 1.0/sigmay759 - 1.0*pow(-my759 + y759, 2)/pow(sigmay759, 3);
            break;
        case 760:
            dJydsigmay[760] = 1.0/sigmay760 - 1.0*pow(-my760 + y760, 2)/pow(sigmay760, 3);
            break;
        case 761:
            dJydsigmay[761] = 1.0/sigmay761 - 1.0*pow(-my761 + y761, 2)/pow(sigmay761, 3);
            break;
        case 762:
            dJydsigmay[762] = 1.0/sigmay762 - 1.0*pow(-my762 + y762, 2)/pow(sigmay762, 3);
            break;
        case 763:
            dJydsigmay[763] = 1.0/sigmay763 - 1.0*pow(-my763 + y763, 2)/pow(sigmay763, 3);
            break;
        case 764:
            dJydsigmay[764] = 1.0/sigmay764 - 1.0*pow(-my764 + y764, 2)/pow(sigmay764, 3);
            break;
        case 765:
            dJydsigmay[765] = 1.0/sigmay765 - 1.0*pow(-my765 + y765, 2)/pow(sigmay765, 3);
            break;
        case 766:
            dJydsigmay[766] = 1.0/sigmay766 - 1.0*pow(-my766 + y766, 2)/pow(sigmay766, 3);
            break;
        case 767:
            dJydsigmay[767] = 1.0/sigmay767 - 1.0*pow(-my767 + y767, 2)/pow(sigmay767, 3);
            break;
        case 768:
            dJydsigmay[768] = 1.0/sigmay768 - 1.0*pow(-my768 + y768, 2)/pow(sigmay768, 3);
            break;
        case 769:
            dJydsigmay[769] = 1.0/sigmay769 - 1.0*pow(-my769 + y769, 2)/pow(sigmay769, 3);
            break;
        case 770:
            dJydsigmay[770] = 1.0/sigmay770 - 1.0*pow(-my770 + y770, 2)/pow(sigmay770, 3);
            break;
        case 771:
            dJydsigmay[771] = 1.0/sigmay771 - 1.0*pow(-my771 + y771, 2)/pow(sigmay771, 3);
            break;
        case 772:
            dJydsigmay[772] = 1.0/sigmay772 - 1.0*pow(-my772 + y772, 2)/pow(sigmay772, 3);
            break;
        case 773:
            dJydsigmay[773] = 1.0/sigmay773 - 1.0*pow(-my773 + y773, 2)/pow(sigmay773, 3);
            break;
        case 774:
            dJydsigmay[774] = 1.0/sigmay774 - 1.0*pow(-my774 + y774, 2)/pow(sigmay774, 3);
            break;
        case 775:
            dJydsigmay[775] = 1.0/sigmay775 - 1.0*pow(-my775 + y775, 2)/pow(sigmay775, 3);
            break;
        case 776:
            dJydsigmay[776] = 1.0/sigmay776 - 1.0*pow(-my776 + y776, 2)/pow(sigmay776, 3);
            break;
        case 777:
            dJydsigmay[777] = 1.0/sigmay777 - 1.0*pow(-my777 + y777, 2)/pow(sigmay777, 3);
            break;
        case 778:
            dJydsigmay[778] = 1.0/sigmay778 - 1.0*pow(-my778 + y778, 2)/pow(sigmay778, 3);
            break;
        case 779:
            dJydsigmay[779] = 1.0/sigmay779 - 1.0*pow(-my779 + y779, 2)/pow(sigmay779, 3);
            break;
        case 780:
            dJydsigmay[780] = 1.0/sigmay780 - 1.0*pow(-my780 + y780, 2)/pow(sigmay780, 3);
            break;
        case 781:
            dJydsigmay[781] = 1.0/sigmay781 - 1.0*pow(-my781 + y781, 2)/pow(sigmay781, 3);
            break;
        case 782:
            dJydsigmay[782] = 1.0/sigmay782 - 1.0*pow(-my782 + y782, 2)/pow(sigmay782, 3);
            break;
        case 783:
            dJydsigmay[783] = 1.0/sigmay783 - 1.0*pow(-my783 + y783, 2)/pow(sigmay783, 3);
            break;
        case 784:
            dJydsigmay[784] = 1.0/sigmay784 - 1.0*pow(-my784 + y784, 2)/pow(sigmay784, 3);
            break;
        case 785:
            dJydsigmay[785] = 1.0/sigmay785 - 1.0*pow(-my785 + y785, 2)/pow(sigmay785, 3);
            break;
        case 786:
            dJydsigmay[786] = 1.0/sigmay786 - 1.0*pow(-my786 + y786, 2)/pow(sigmay786, 3);
            break;
        case 787:
            dJydsigmay[787] = 1.0/sigmay787 - 1.0*pow(-my787 + y787, 2)/pow(sigmay787, 3);
            break;
        case 788:
            dJydsigmay[788] = 1.0/sigmay788 - 1.0*pow(-my788 + y788, 2)/pow(sigmay788, 3);
            break;
        case 789:
            dJydsigmay[789] = 1.0/sigmay789 - 1.0*pow(-my789 + y789, 2)/pow(sigmay789, 3);
            break;
        case 790:
            dJydsigmay[790] = 1.0/sigmay790 - 1.0*pow(-my790 + y790, 2)/pow(sigmay790, 3);
            break;
        case 791:
            dJydsigmay[791] = 1.0/sigmay791 - 1.0*pow(-my791 + y791, 2)/pow(sigmay791, 3);
            break;
        case 792:
            dJydsigmay[792] = 1.0/sigmay792 - 1.0*pow(-my792 + y792, 2)/pow(sigmay792, 3);
            break;
        case 793:
            dJydsigmay[793] = 1.0/sigmay793 - 1.0*pow(-my793 + y793, 2)/pow(sigmay793, 3);
            break;
        case 794:
            dJydsigmay[794] = 1.0/sigmay794 - 1.0*pow(-my794 + y794, 2)/pow(sigmay794, 3);
            break;
        case 795:
            dJydsigmay[795] = 1.0/sigmay795 - 1.0*pow(-my795 + y795, 2)/pow(sigmay795, 3);
            break;
        case 796:
            dJydsigmay[796] = 1.0/sigmay796 - 1.0*pow(-my796 + y796, 2)/pow(sigmay796, 3);
            break;
        case 797:
            dJydsigmay[797] = 1.0/sigmay797 - 1.0*pow(-my797 + y797, 2)/pow(sigmay797, 3);
            break;
        case 798:
            dJydsigmay[798] = 1.0/sigmay798 - 1.0*pow(-my798 + y798, 2)/pow(sigmay798, 3);
            break;
        case 799:
            dJydsigmay[799] = 1.0/sigmay799 - 1.0*pow(-my799 + y799, 2)/pow(sigmay799, 3);
            break;
        case 800:
            dJydsigmay[800] = 1.0/sigmay800 - 1.0*pow(-my800 + y800, 2)/pow(sigmay800, 3);
            break;
        case 801:
            dJydsigmay[801] = 1.0/sigmay801 - 1.0*pow(-my801 + y801, 2)/pow(sigmay801, 3);
            break;
        case 802:
            dJydsigmay[802] = 1.0/sigmay802 - 1.0*pow(-my802 + y802, 2)/pow(sigmay802, 3);
            break;
        case 803:
            dJydsigmay[803] = 1.0/sigmay803 - 1.0*pow(-my803 + y803, 2)/pow(sigmay803, 3);
            break;
        case 804:
            dJydsigmay[804] = 1.0/sigmay804 - 1.0*pow(-my804 + y804, 2)/pow(sigmay804, 3);
            break;
        case 805:
            dJydsigmay[805] = 1.0/sigmay805 - 1.0*pow(-my805 + y805, 2)/pow(sigmay805, 3);
            break;
        case 806:
            dJydsigmay[806] = 1.0/sigmay806 - 1.0*pow(-my806 + y806, 2)/pow(sigmay806, 3);
            break;
        case 807:
            dJydsigmay[807] = 1.0/sigmay807 - 1.0*pow(-my807 + y807, 2)/pow(sigmay807, 3);
            break;
        case 808:
            dJydsigmay[808] = 1.0/sigmay808 - 1.0*pow(-my808 + y808, 2)/pow(sigmay808, 3);
            break;
        case 809:
            dJydsigmay[809] = 1.0/sigmay809 - 1.0*pow(-my809 + y809, 2)/pow(sigmay809, 3);
            break;
        case 810:
            dJydsigmay[810] = 1.0/sigmay810 - 1.0*pow(-my810 + y810, 2)/pow(sigmay810, 3);
            break;
        case 811:
            dJydsigmay[811] = 1.0/sigmay811 - 1.0*pow(-my811 + y811, 2)/pow(sigmay811, 3);
            break;
        case 812:
            dJydsigmay[812] = 1.0/sigmay812 - 1.0*pow(-my812 + y812, 2)/pow(sigmay812, 3);
            break;
        case 813:
            dJydsigmay[813] = 1.0/sigmay813 - 1.0*pow(-my813 + y813, 2)/pow(sigmay813, 3);
            break;
        case 814:
            dJydsigmay[814] = 1.0/sigmay814 - 1.0*pow(-my814 + y814, 2)/pow(sigmay814, 3);
            break;
        case 815:
            dJydsigmay[815] = 1.0/sigmay815 - 1.0*pow(-my815 + y815, 2)/pow(sigmay815, 3);
            break;
        case 816:
            dJydsigmay[816] = 1.0/sigmay816 - 1.0*pow(-my816 + y816, 2)/pow(sigmay816, 3);
            break;
        case 817:
            dJydsigmay[817] = 1.0/sigmay817 - 1.0*pow(-my817 + y817, 2)/pow(sigmay817, 3);
            break;
        case 818:
            dJydsigmay[818] = 1.0/sigmay818 - 1.0*pow(-my818 + y818, 2)/pow(sigmay818, 3);
            break;
        case 819:
            dJydsigmay[819] = 1.0/sigmay819 - 1.0*pow(-my819 + y819, 2)/pow(sigmay819, 3);
            break;
        case 820:
            dJydsigmay[820] = 1.0/sigmay820 - 1.0*pow(-my820 + y820, 2)/pow(sigmay820, 3);
            break;
        case 821:
            dJydsigmay[821] = 1.0/sigmay821 - 1.0*pow(-my821 + y821, 2)/pow(sigmay821, 3);
            break;
        case 822:
            dJydsigmay[822] = 1.0/sigmay822 - 1.0*pow(-my822 + y822, 2)/pow(sigmay822, 3);
            break;
        case 823:
            dJydsigmay[823] = 1.0/sigmay823 - 1.0*pow(-my823 + y823, 2)/pow(sigmay823, 3);
            break;
        case 824:
            dJydsigmay[824] = 1.0/sigmay824 - 1.0*pow(-my824 + y824, 2)/pow(sigmay824, 3);
            break;
        case 825:
            dJydsigmay[825] = 1.0/sigmay825 - 1.0*pow(-my825 + y825, 2)/pow(sigmay825, 3);
            break;
        case 826:
            dJydsigmay[826] = 1.0/sigmay826 - 1.0*pow(-my826 + y826, 2)/pow(sigmay826, 3);
            break;
        case 827:
            dJydsigmay[827] = 1.0/sigmay827 - 1.0*pow(-my827 + y827, 2)/pow(sigmay827, 3);
            break;
        case 828:
            dJydsigmay[828] = 1.0/sigmay828 - 1.0*pow(-my828 + y828, 2)/pow(sigmay828, 3);
            break;
        case 829:
            dJydsigmay[829] = 1.0/sigmay829 - 1.0*pow(-my829 + y829, 2)/pow(sigmay829, 3);
            break;
        case 830:
            dJydsigmay[830] = 1.0/sigmay830 - 1.0*pow(-my830 + y830, 2)/pow(sigmay830, 3);
            break;
        case 831:
            dJydsigmay[831] = 1.0/sigmay831 - 1.0*pow(-my831 + y831, 2)/pow(sigmay831, 3);
            break;
        case 832:
            dJydsigmay[832] = 1.0/sigmay832 - 1.0*pow(-my832 + y832, 2)/pow(sigmay832, 3);
            break;
        case 833:
            dJydsigmay[833] = 1.0/sigmay833 - 1.0*pow(-my833 + y833, 2)/pow(sigmay833, 3);
            break;
        case 834:
            dJydsigmay[834] = 1.0/sigmay834 - 1.0*pow(-my834 + y834, 2)/pow(sigmay834, 3);
            break;
        case 835:
            dJydsigmay[835] = 1.0/sigmay835 - 1.0*pow(-my835 + y835, 2)/pow(sigmay835, 3);
            break;
        case 836:
            dJydsigmay[836] = 1.0/sigmay836 - 1.0*pow(-my836 + y836, 2)/pow(sigmay836, 3);
            break;
        case 837:
            dJydsigmay[837] = 1.0/sigmay837 - 1.0*pow(-my837 + y837, 2)/pow(sigmay837, 3);
            break;
        case 838:
            dJydsigmay[838] = 1.0/sigmay838 - 1.0*pow(-my838 + y838, 2)/pow(sigmay838, 3);
            break;
        case 839:
            dJydsigmay[839] = 1.0/sigmay839 - 1.0*pow(-my839 + y839, 2)/pow(sigmay839, 3);
            break;
        case 840:
            dJydsigmay[840] = 1.0/sigmay840 - 1.0*pow(-my840 + y840, 2)/pow(sigmay840, 3);
            break;
        case 841:
            dJydsigmay[841] = 1.0/sigmay841 - 1.0*pow(-my841 + y841, 2)/pow(sigmay841, 3);
            break;
        case 842:
            dJydsigmay[842] = 1.0/sigmay842 - 1.0*pow(-my842 + y842, 2)/pow(sigmay842, 3);
            break;
        case 843:
            dJydsigmay[843] = 1.0/sigmay843 - 1.0*pow(-my843 + y843, 2)/pow(sigmay843, 3);
            break;
        case 844:
            dJydsigmay[844] = 1.0/sigmay844 - 1.0*pow(-my844 + y844, 2)/pow(sigmay844, 3);
            break;
        case 845:
            dJydsigmay[845] = 1.0/sigmay845 - 1.0*pow(-my845 + y845, 2)/pow(sigmay845, 3);
            break;
        case 846:
            dJydsigmay[846] = 1.0/sigmay846 - 1.0*pow(-my846 + y846, 2)/pow(sigmay846, 3);
            break;
        case 847:
            dJydsigmay[847] = 1.0/sigmay847 - 1.0*pow(-my847 + y847, 2)/pow(sigmay847, 3);
            break;
        case 848:
            dJydsigmay[848] = 1.0/sigmay848 - 1.0*pow(-my848 + y848, 2)/pow(sigmay848, 3);
            break;
        case 849:
            dJydsigmay[849] = 1.0/sigmay849 - 1.0*pow(-my849 + y849, 2)/pow(sigmay849, 3);
            break;
        case 850:
            dJydsigmay[850] = 1.0/sigmay850 - 1.0*pow(-my850 + y850, 2)/pow(sigmay850, 3);
            break;
        case 851:
            dJydsigmay[851] = 1.0/sigmay851 - 1.0*pow(-my851 + y851, 2)/pow(sigmay851, 3);
            break;
        case 852:
            dJydsigmay[852] = 1.0/sigmay852 - 1.0*pow(-my852 + y852, 2)/pow(sigmay852, 3);
            break;
        case 853:
            dJydsigmay[853] = 1.0/sigmay853 - 1.0*pow(-my853 + y853, 2)/pow(sigmay853, 3);
            break;
        case 854:
            dJydsigmay[854] = 1.0/sigmay854 - 1.0*pow(-my854 + y854, 2)/pow(sigmay854, 3);
            break;
        case 855:
            dJydsigmay[855] = 1.0/sigmay855 - 1.0*pow(-my855 + y855, 2)/pow(sigmay855, 3);
            break;
        case 856:
            dJydsigmay[856] = 1.0/sigmay856 - 1.0*pow(-my856 + y856, 2)/pow(sigmay856, 3);
            break;
        case 857:
            dJydsigmay[857] = 1.0/sigmay857 - 1.0*pow(-my857 + y857, 2)/pow(sigmay857, 3);
            break;
        case 858:
            dJydsigmay[858] = 1.0/sigmay858 - 1.0*pow(-my858 + y858, 2)/pow(sigmay858, 3);
            break;
        case 859:
            dJydsigmay[859] = 1.0/sigmay859 - 1.0*pow(-my859 + y859, 2)/pow(sigmay859, 3);
            break;
        case 860:
            dJydsigmay[860] = 1.0/sigmay860 - 1.0*pow(-my860 + y860, 2)/pow(sigmay860, 3);
            break;
        case 861:
            dJydsigmay[861] = 1.0/sigmay861 - 1.0*pow(-my861 + y861, 2)/pow(sigmay861, 3);
            break;
        case 862:
            dJydsigmay[862] = 1.0/sigmay862 - 1.0*pow(-my862 + y862, 2)/pow(sigmay862, 3);
            break;
        case 863:
            dJydsigmay[863] = 1.0/sigmay863 - 1.0*pow(-my863 + y863, 2)/pow(sigmay863, 3);
            break;
        case 864:
            dJydsigmay[864] = 1.0/sigmay864 - 1.0*pow(-my864 + y864, 2)/pow(sigmay864, 3);
            break;
        case 865:
            dJydsigmay[865] = 1.0/sigmay865 - 1.0*pow(-my865 + y865, 2)/pow(sigmay865, 3);
            break;
        case 866:
            dJydsigmay[866] = 1.0/sigmay866 - 1.0*pow(-my866 + y866, 2)/pow(sigmay866, 3);
            break;
        case 867:
            dJydsigmay[867] = 1.0/sigmay867 - 1.0*pow(-my867 + y867, 2)/pow(sigmay867, 3);
            break;
        case 868:
            dJydsigmay[868] = 1.0/sigmay868 - 1.0*pow(-my868 + y868, 2)/pow(sigmay868, 3);
            break;
        case 869:
            dJydsigmay[869] = 1.0/sigmay869 - 1.0*pow(-my869 + y869, 2)/pow(sigmay869, 3);
            break;
        case 870:
            dJydsigmay[870] = 1.0/sigmay870 - 1.0*pow(-my870 + y870, 2)/pow(sigmay870, 3);
            break;
        case 871:
            dJydsigmay[871] = 1.0/sigmay871 - 1.0*pow(-my871 + y871, 2)/pow(sigmay871, 3);
            break;
        case 872:
            dJydsigmay[872] = 1.0/sigmay872 - 1.0*pow(-my872 + y872, 2)/pow(sigmay872, 3);
            break;
        case 873:
            dJydsigmay[873] = 1.0/sigmay873 - 1.0*pow(-my873 + y873, 2)/pow(sigmay873, 3);
            break;
        case 874:
            dJydsigmay[874] = 1.0/sigmay874 - 1.0*pow(-my874 + y874, 2)/pow(sigmay874, 3);
            break;
        case 875:
            dJydsigmay[875] = 1.0/sigmay875 - 1.0*pow(-my875 + y875, 2)/pow(sigmay875, 3);
            break;
        case 876:
            dJydsigmay[876] = 1.0/sigmay876 - 1.0*pow(-my876 + y876, 2)/pow(sigmay876, 3);
            break;
        case 877:
            dJydsigmay[877] = 1.0/sigmay877 - 1.0*pow(-my877 + y877, 2)/pow(sigmay877, 3);
            break;
        case 878:
            dJydsigmay[878] = 1.0/sigmay878 - 1.0*pow(-my878 + y878, 2)/pow(sigmay878, 3);
            break;
        case 879:
            dJydsigmay[879] = 1.0/sigmay879 - 1.0*pow(-my879 + y879, 2)/pow(sigmay879, 3);
            break;
        case 880:
            dJydsigmay[880] = 1.0/sigmay880 - 1.0*pow(-my880 + y880, 2)/pow(sigmay880, 3);
            break;
        case 881:
            dJydsigmay[881] = 1.0/sigmay881 - 1.0*pow(-my881 + y881, 2)/pow(sigmay881, 3);
            break;
        case 882:
            dJydsigmay[882] = 1.0/sigmay882 - 1.0*pow(-my882 + y882, 2)/pow(sigmay882, 3);
            break;
        case 883:
            dJydsigmay[883] = 1.0/sigmay883 - 1.0*pow(-my883 + y883, 2)/pow(sigmay883, 3);
            break;
        case 884:
            dJydsigmay[884] = 1.0/sigmay884 - 1.0*pow(-my884 + y884, 2)/pow(sigmay884, 3);
            break;
        case 885:
            dJydsigmay[885] = 1.0/sigmay885 - 1.0*pow(-my885 + y885, 2)/pow(sigmay885, 3);
            break;
        case 886:
            dJydsigmay[886] = 1.0/sigmay886 - 1.0*pow(-my886 + y886, 2)/pow(sigmay886, 3);
            break;
        case 887:
            dJydsigmay[887] = 1.0/sigmay887 - 1.0*pow(-my887 + y887, 2)/pow(sigmay887, 3);
            break;
        case 888:
            dJydsigmay[888] = 1.0/sigmay888 - 1.0*pow(-my888 + y888, 2)/pow(sigmay888, 3);
            break;
        case 889:
            dJydsigmay[889] = 1.0/sigmay889 - 1.0*pow(-my889 + y889, 2)/pow(sigmay889, 3);
            break;
        case 890:
            dJydsigmay[890] = 1.0/sigmay890 - 1.0*pow(-my890 + y890, 2)/pow(sigmay890, 3);
            break;
        case 891:
            dJydsigmay[891] = 1.0/sigmay891 - 1.0*pow(-my891 + y891, 2)/pow(sigmay891, 3);
            break;
        case 892:
            dJydsigmay[892] = 1.0/sigmay892 - 1.0*pow(-my892 + y892, 2)/pow(sigmay892, 3);
            break;
        case 893:
            dJydsigmay[893] = 1.0/sigmay893 - 1.0*pow(-my893 + y893, 2)/pow(sigmay893, 3);
            break;
        case 894:
            dJydsigmay[894] = 1.0/sigmay894 - 1.0*pow(-my894 + y894, 2)/pow(sigmay894, 3);
            break;
        case 895:
            dJydsigmay[895] = 1.0/sigmay895 - 1.0*pow(-my895 + y895, 2)/pow(sigmay895, 3);
            break;
        case 896:
            dJydsigmay[896] = 1.0/sigmay896 - 1.0*pow(-my896 + y896, 2)/pow(sigmay896, 3);
            break;
        case 897:
            dJydsigmay[897] = 1.0/sigmay897 - 1.0*pow(-my897 + y897, 2)/pow(sigmay897, 3);
            break;
        case 898:
            dJydsigmay[898] = 1.0/sigmay898 - 1.0*pow(-my898 + y898, 2)/pow(sigmay898, 3);
            break;
        case 899:
            dJydsigmay[899] = 1.0/sigmay899 - 1.0*pow(-my899 + y899, 2)/pow(sigmay899, 3);
            break;
        case 900:
            dJydsigmay[900] = 1.0/sigmay900 - 1.0*pow(-my900 + y900, 2)/pow(sigmay900, 3);
            break;
        case 901:
            dJydsigmay[901] = 1.0/sigmay901 - 1.0*pow(-my901 + y901, 2)/pow(sigmay901, 3);
            break;
        case 902:
            dJydsigmay[902] = 1.0/sigmay902 - 1.0*pow(-my902 + y902, 2)/pow(sigmay902, 3);
            break;
        case 903:
            dJydsigmay[903] = 1.0/sigmay903 - 1.0*pow(-my903 + y903, 2)/pow(sigmay903, 3);
            break;
        case 904:
            dJydsigmay[904] = 1.0/sigmay904 - 1.0*pow(-my904 + y904, 2)/pow(sigmay904, 3);
            break;
        case 905:
            dJydsigmay[905] = 1.0/sigmay905 - 1.0*pow(-my905 + y905, 2)/pow(sigmay905, 3);
            break;
        case 906:
            dJydsigmay[906] = 1.0/sigmay906 - 1.0*pow(-my906 + y906, 2)/pow(sigmay906, 3);
            break;
        case 907:
            dJydsigmay[907] = 1.0/sigmay907 - 1.0*pow(-my907 + y907, 2)/pow(sigmay907, 3);
            break;
        case 908:
            dJydsigmay[908] = 1.0/sigmay908 - 1.0*pow(-my908 + y908, 2)/pow(sigmay908, 3);
            break;
        case 909:
            dJydsigmay[909] = 1.0/sigmay909 - 1.0*pow(-my909 + y909, 2)/pow(sigmay909, 3);
            break;
        case 910:
            dJydsigmay[910] = 1.0/sigmay910 - 1.0*pow(-my910 + y910, 2)/pow(sigmay910, 3);
            break;
        case 911:
            dJydsigmay[911] = 1.0/sigmay911 - 1.0*pow(-my911 + y911, 2)/pow(sigmay911, 3);
            break;
        case 912:
            dJydsigmay[912] = 1.0/sigmay912 - 1.0*pow(-my912 + y912, 2)/pow(sigmay912, 3);
            break;
        case 913:
            dJydsigmay[913] = 1.0/sigmay913 - 1.0*pow(-my913 + y913, 2)/pow(sigmay913, 3);
            break;
        case 914:
            dJydsigmay[914] = 1.0/sigmay914 - 1.0*pow(-my914 + y914, 2)/pow(sigmay914, 3);
            break;
        case 915:
            dJydsigmay[915] = 1.0/sigmay915 - 1.0*pow(-my915 + y915, 2)/pow(sigmay915, 3);
            break;
        case 916:
            dJydsigmay[916] = 1.0/sigmay916 - 1.0*pow(-my916 + y916, 2)/pow(sigmay916, 3);
            break;
        case 917:
            dJydsigmay[917] = 1.0/sigmay917 - 1.0*pow(-my917 + y917, 2)/pow(sigmay917, 3);
            break;
        case 918:
            dJydsigmay[918] = 1.0/sigmay918 - 1.0*pow(-my918 + y918, 2)/pow(sigmay918, 3);
            break;
        case 919:
            dJydsigmay[919] = 1.0/sigmay919 - 1.0*pow(-my919 + y919, 2)/pow(sigmay919, 3);
            break;
        case 920:
            dJydsigmay[920] = 1.0/sigmay920 - 1.0*pow(-my920 + y920, 2)/pow(sigmay920, 3);
            break;
        case 921:
            dJydsigmay[921] = 1.0/sigmay921 - 1.0*pow(-my921 + y921, 2)/pow(sigmay921, 3);
            break;
        case 922:
            dJydsigmay[922] = 1.0/sigmay922 - 1.0*pow(-my922 + y922, 2)/pow(sigmay922, 3);
            break;
        case 923:
            dJydsigmay[923] = 1.0/sigmay923 - 1.0*pow(-my923 + y923, 2)/pow(sigmay923, 3);
            break;
        case 924:
            dJydsigmay[924] = 1.0/sigmay924 - 1.0*pow(-my924 + y924, 2)/pow(sigmay924, 3);
            break;
        case 925:
            dJydsigmay[925] = 1.0/sigmay925 - 1.0*pow(-my925 + y925, 2)/pow(sigmay925, 3);
            break;
        case 926:
            dJydsigmay[926] = 1.0/sigmay926 - 1.0*pow(-my926 + y926, 2)/pow(sigmay926, 3);
            break;
        case 927:
            dJydsigmay[927] = 1.0/sigmay927 - 1.0*pow(-my927 + y927, 2)/pow(sigmay927, 3);
            break;
        case 928:
            dJydsigmay[928] = 1.0/sigmay928 - 1.0*pow(-my928 + y928, 2)/pow(sigmay928, 3);
            break;
        case 929:
            dJydsigmay[929] = 1.0/sigmay929 - 1.0*pow(-my929 + y929, 2)/pow(sigmay929, 3);
            break;
        case 930:
            dJydsigmay[930] = 1.0/sigmay930 - 1.0*pow(-my930 + y930, 2)/pow(sigmay930, 3);
            break;
        case 931:
            dJydsigmay[931] = 1.0/sigmay931 - 1.0*pow(-my931 + y931, 2)/pow(sigmay931, 3);
            break;
        case 932:
            dJydsigmay[932] = 1.0/sigmay932 - 1.0*pow(-my932 + y932, 2)/pow(sigmay932, 3);
            break;
        case 933:
            dJydsigmay[933] = 1.0/sigmay933 - 1.0*pow(-my933 + y933, 2)/pow(sigmay933, 3);
            break;
        case 934:
            dJydsigmay[934] = 1.0/sigmay934 - 1.0*pow(-my934 + y934, 2)/pow(sigmay934, 3);
            break;
        case 935:
            dJydsigmay[935] = 1.0/sigmay935 - 1.0*pow(-my935 + y935, 2)/pow(sigmay935, 3);
            break;
        case 936:
            dJydsigmay[936] = 1.0/sigmay936 - 1.0*pow(-my936 + y936, 2)/pow(sigmay936, 3);
            break;
        case 937:
            dJydsigmay[937] = 1.0/sigmay937 - 1.0*pow(-my937 + y937, 2)/pow(sigmay937, 3);
            break;
        case 938:
            dJydsigmay[938] = 1.0/sigmay938 - 1.0*pow(-my938 + y938, 2)/pow(sigmay938, 3);
            break;
        case 939:
            dJydsigmay[939] = 1.0/sigmay939 - 1.0*pow(-my939 + y939, 2)/pow(sigmay939, 3);
            break;
        case 940:
            dJydsigmay[940] = 1.0/sigmay940 - 1.0*pow(-my940 + y940, 2)/pow(sigmay940, 3);
            break;
        case 941:
            dJydsigmay[941] = 1.0/sigmay941 - 1.0*pow(-my941 + y941, 2)/pow(sigmay941, 3);
            break;
        case 942:
            dJydsigmay[942] = 1.0/sigmay942 - 1.0*pow(-my942 + y942, 2)/pow(sigmay942, 3);
            break;
        case 943:
            dJydsigmay[943] = 1.0/sigmay943 - 1.0*pow(-my943 + y943, 2)/pow(sigmay943, 3);
            break;
        case 944:
            dJydsigmay[944] = 1.0/sigmay944 - 1.0*pow(-my944 + y944, 2)/pow(sigmay944, 3);
            break;
        case 945:
            dJydsigmay[945] = 1.0/sigmay945 - 1.0*pow(-my945 + y945, 2)/pow(sigmay945, 3);
            break;
        case 946:
            dJydsigmay[946] = 1.0/sigmay946 - 1.0*pow(-my946 + y946, 2)/pow(sigmay946, 3);
            break;
        case 947:
            dJydsigmay[947] = 1.0/sigmay947 - 1.0*pow(-my947 + y947, 2)/pow(sigmay947, 3);
            break;
        case 948:
            dJydsigmay[948] = 1.0/sigmay948 - 1.0*pow(-my948 + y948, 2)/pow(sigmay948, 3);
            break;
        case 949:
            dJydsigmay[949] = 1.0/sigmay949 - 1.0*pow(-my949 + y949, 2)/pow(sigmay949, 3);
            break;
        case 950:
            dJydsigmay[950] = 1.0/sigmay950 - 1.0*pow(-my950 + y950, 2)/pow(sigmay950, 3);
            break;
        case 951:
            dJydsigmay[951] = 1.0/sigmay951 - 1.0*pow(-my951 + y951, 2)/pow(sigmay951, 3);
            break;
        case 952:
            dJydsigmay[952] = 1.0/sigmay952 - 1.0*pow(-my952 + y952, 2)/pow(sigmay952, 3);
            break;
        case 953:
            dJydsigmay[953] = 1.0/sigmay953 - 1.0*pow(-my953 + y953, 2)/pow(sigmay953, 3);
            break;
        case 954:
            dJydsigmay[954] = 1.0/sigmay954 - 1.0*pow(-my954 + y954, 2)/pow(sigmay954, 3);
            break;
        case 955:
            dJydsigmay[955] = 1.0/sigmay955 - 1.0*pow(-my955 + y955, 2)/pow(sigmay955, 3);
            break;
        case 956:
            dJydsigmay[956] = 1.0/sigmay956 - 1.0*pow(-my956 + y956, 2)/pow(sigmay956, 3);
            break;
        case 957:
            dJydsigmay[957] = 1.0/sigmay957 - 1.0*pow(-my957 + y957, 2)/pow(sigmay957, 3);
            break;
        case 958:
            dJydsigmay[958] = 1.0/sigmay958 - 1.0*pow(-my958 + y958, 2)/pow(sigmay958, 3);
            break;
        case 959:
            dJydsigmay[959] = 1.0/sigmay959 - 1.0*pow(-my959 + y959, 2)/pow(sigmay959, 3);
            break;
        case 960:
            dJydsigmay[960] = 1.0/sigmay960 - 1.0*pow(-my960 + y960, 2)/pow(sigmay960, 3);
            break;
        case 961:
            dJydsigmay[961] = 1.0/sigmay961 - 1.0*pow(-my961 + y961, 2)/pow(sigmay961, 3);
            break;
        case 962:
            dJydsigmay[962] = 1.0/sigmay962 - 1.0*pow(-my962 + y962, 2)/pow(sigmay962, 3);
            break;
        case 963:
            dJydsigmay[963] = 1.0/sigmay963 - 1.0*pow(-my963 + y963, 2)/pow(sigmay963, 3);
            break;
        case 964:
            dJydsigmay[964] = 1.0/sigmay964 - 1.0*pow(-my964 + y964, 2)/pow(sigmay964, 3);
            break;
        case 965:
            dJydsigmay[965] = 1.0/sigmay965 - 1.0*pow(-my965 + y965, 2)/pow(sigmay965, 3);
            break;
        case 966:
            dJydsigmay[966] = 1.0/sigmay966 - 1.0*pow(-my966 + y966, 2)/pow(sigmay966, 3);
            break;
        case 967:
            dJydsigmay[967] = 1.0/sigmay967 - 1.0*pow(-my967 + y967, 2)/pow(sigmay967, 3);
            break;
        case 968:
            dJydsigmay[968] = 1.0/sigmay968 - 1.0*pow(-my968 + y968, 2)/pow(sigmay968, 3);
            break;
        case 969:
            dJydsigmay[969] = 1.0/sigmay969 - 1.0*pow(-my969 + y969, 2)/pow(sigmay969, 3);
            break;
        case 970:
            dJydsigmay[970] = 1.0/sigmay970 - 1.0*pow(-my970 + y970, 2)/pow(sigmay970, 3);
            break;
        case 971:
            dJydsigmay[971] = 1.0/sigmay971 - 1.0*pow(-my971 + y971, 2)/pow(sigmay971, 3);
            break;
        case 972:
            dJydsigmay[972] = 1.0/sigmay972 - 1.0*pow(-my972 + y972, 2)/pow(sigmay972, 3);
            break;
        case 973:
            dJydsigmay[973] = 1.0/sigmay973 - 1.0*pow(-my973 + y973, 2)/pow(sigmay973, 3);
            break;
        case 974:
            dJydsigmay[974] = 1.0/sigmay974 - 1.0*pow(-my974 + y974, 2)/pow(sigmay974, 3);
            break;
        case 975:
            dJydsigmay[975] = 1.0/sigmay975 - 1.0*pow(-my975 + y975, 2)/pow(sigmay975, 3);
            break;
        case 976:
            dJydsigmay[976] = 1.0/sigmay976 - 1.0*pow(-my976 + y976, 2)/pow(sigmay976, 3);
            break;
        case 977:
            dJydsigmay[977] = 1.0/sigmay977 - 1.0*pow(-my977 + y977, 2)/pow(sigmay977, 3);
            break;
        case 978:
            dJydsigmay[978] = 1.0/sigmay978 - 1.0*pow(-my978 + y978, 2)/pow(sigmay978, 3);
            break;
        case 979:
            dJydsigmay[979] = 1.0/sigmay979 - 1.0*pow(-my979 + y979, 2)/pow(sigmay979, 3);
            break;
        case 980:
            dJydsigmay[980] = 1.0/sigmay980 - 1.0*pow(-my980 + y980, 2)/pow(sigmay980, 3);
            break;
        case 981:
            dJydsigmay[981] = 1.0/sigmay981 - 1.0*pow(-my981 + y981, 2)/pow(sigmay981, 3);
            break;
        case 982:
            dJydsigmay[982] = 1.0/sigmay982 - 1.0*pow(-my982 + y982, 2)/pow(sigmay982, 3);
            break;
        case 983:
            dJydsigmay[983] = 1.0/sigmay983 - 1.0*pow(-my983 + y983, 2)/pow(sigmay983, 3);
            break;
        case 984:
            dJydsigmay[984] = 1.0/sigmay984 - 1.0*pow(-my984 + y984, 2)/pow(sigmay984, 3);
            break;
        case 985:
            dJydsigmay[985] = 1.0/sigmay985 - 1.0*pow(-my985 + y985, 2)/pow(sigmay985, 3);
            break;
        case 986:
            dJydsigmay[986] = 1.0/sigmay986 - 1.0*pow(-my986 + y986, 2)/pow(sigmay986, 3);
            break;
        case 987:
            dJydsigmay[987] = 1.0/sigmay987 - 1.0*pow(-my987 + y987, 2)/pow(sigmay987, 3);
            break;
        case 988:
            dJydsigmay[988] = 1.0/sigmay988 - 1.0*pow(-my988 + y988, 2)/pow(sigmay988, 3);
            break;
        case 989:
            dJydsigmay[989] = 1.0/sigmay989 - 1.0*pow(-my989 + y989, 2)/pow(sigmay989, 3);
            break;
        case 990:
            dJydsigmay[990] = 1.0/sigmay990 - 1.0*pow(-my990 + y990, 2)/pow(sigmay990, 3);
            break;
        case 991:
            dJydsigmay[991] = 1.0/sigmay991 - 1.0*pow(-my991 + y991, 2)/pow(sigmay991, 3);
            break;
        case 992:
            dJydsigmay[992] = 1.0/sigmay992 - 1.0*pow(-my992 + y992, 2)/pow(sigmay992, 3);
            break;
        case 993:
            dJydsigmay[993] = 1.0/sigmay993 - 1.0*pow(-my993 + y993, 2)/pow(sigmay993, 3);
            break;
        case 994:
            dJydsigmay[994] = 1.0/sigmay994 - 1.0*pow(-my994 + y994, 2)/pow(sigmay994, 3);
            break;
        case 995:
            dJydsigmay[995] = 1.0/sigmay995 - 1.0*pow(-my995 + y995, 2)/pow(sigmay995, 3);
            break;
        case 996:
            dJydsigmay[996] = 1.0/sigmay996 - 1.0*pow(-my996 + y996, 2)/pow(sigmay996, 3);
            break;
        case 997:
            dJydsigmay[997] = 1.0/sigmay997 - 1.0*pow(-my997 + y997, 2)/pow(sigmay997, 3);
            break;
        case 998:
            dJydsigmay[998] = 1.0/sigmay998 - 1.0*pow(-my998 + y998, 2)/pow(sigmay998, 3);
            break;
        case 999:
            dJydsigmay[999] = 1.0/sigmay999 - 1.0*pow(-my999 + y999, 2)/pow(sigmay999, 3);
            break;
        case 1000:
            dJydsigmay[1000] = 1.0/sigmay1000 - 1.0*pow(-my1000 + y1000, 2)/pow(sigmay1000, 3);
            break;
        case 1001:
            dJydsigmay[1001] = 1.0/sigmay1001 - 1.0*pow(-my1001 + y1001, 2)/pow(sigmay1001, 3);
            break;
        case 1002:
            dJydsigmay[1002] = 1.0/sigmay1002 - 1.0*pow(-my1002 + y1002, 2)/pow(sigmay1002, 3);
            break;
        case 1003:
            dJydsigmay[1003] = 1.0/sigmay1003 - 1.0*pow(-my1003 + y1003, 2)/pow(sigmay1003, 3);
            break;
        case 1004:
            dJydsigmay[1004] = 1.0/sigmay1004 - 1.0*pow(-my1004 + y1004, 2)/pow(sigmay1004, 3);
            break;
        case 1005:
            dJydsigmay[1005] = 1.0/sigmay1005 - 1.0*pow(-my1005 + y1005, 2)/pow(sigmay1005, 3);
            break;
        case 1006:
            dJydsigmay[1006] = 1.0/sigmay1006 - 1.0*pow(-my1006 + y1006, 2)/pow(sigmay1006, 3);
            break;
        case 1007:
            dJydsigmay[1007] = 1.0/sigmay1007 - 1.0*pow(-my1007 + y1007, 2)/pow(sigmay1007, 3);
            break;
        case 1008:
            dJydsigmay[1008] = 1.0/sigmay1008 - 1.0*pow(-my1008 + y1008, 2)/pow(sigmay1008, 3);
            break;
        case 1009:
            dJydsigmay[1009] = 1.0/sigmay1009 - 1.0*pow(-my1009 + y1009, 2)/pow(sigmay1009, 3);
            break;
        case 1010:
            dJydsigmay[1010] = 1.0/sigmay1010 - 1.0*pow(-my1010 + y1010, 2)/pow(sigmay1010, 3);
            break;
        case 1011:
            dJydsigmay[1011] = 1.0/sigmay1011 - 1.0*pow(-my1011 + y1011, 2)/pow(sigmay1011, 3);
            break;
        case 1012:
            dJydsigmay[1012] = 1.0/sigmay1012 - 1.0*pow(-my1012 + y1012, 2)/pow(sigmay1012, 3);
            break;
        case 1013:
            dJydsigmay[1013] = 1.0/sigmay1013 - 1.0*pow(-my1013 + y1013, 2)/pow(sigmay1013, 3);
            break;
        case 1014:
            dJydsigmay[1014] = 1.0/sigmay1014 - 1.0*pow(-my1014 + y1014, 2)/pow(sigmay1014, 3);
            break;
        case 1015:
            dJydsigmay[1015] = 1.0/sigmay1015 - 1.0*pow(-my1015 + y1015, 2)/pow(sigmay1015, 3);
            break;
        case 1016:
            dJydsigmay[1016] = 1.0/sigmay1016 - 1.0*pow(-my1016 + y1016, 2)/pow(sigmay1016, 3);
            break;
        case 1017:
            dJydsigmay[1017] = 1.0/sigmay1017 - 1.0*pow(-my1017 + y1017, 2)/pow(sigmay1017, 3);
            break;
        case 1018:
            dJydsigmay[1018] = 1.0/sigmay1018 - 1.0*pow(-my1018 + y1018, 2)/pow(sigmay1018, 3);
            break;
        case 1019:
            dJydsigmay[1019] = 1.0/sigmay1019 - 1.0*pow(-my1019 + y1019, 2)/pow(sigmay1019, 3);
            break;
        case 1020:
            dJydsigmay[1020] = 1.0/sigmay1020 - 1.0*pow(-my1020 + y1020, 2)/pow(sigmay1020, 3);
            break;
        case 1021:
            dJydsigmay[1021] = 1.0/sigmay1021 - 1.0*pow(-my1021 + y1021, 2)/pow(sigmay1021, 3);
            break;
        case 1022:
            dJydsigmay[1022] = 1.0/sigmay1022 - 1.0*pow(-my1022 + y1022, 2)/pow(sigmay1022, 3);
            break;
        case 1023:
            dJydsigmay[1023] = 1.0/sigmay1023 - 1.0*pow(-my1023 + y1023, 2)/pow(sigmay1023, 3);
            break;
        case 1024:
            dJydsigmay[1024] = 1.0/sigmay1024 - 1.0*pow(-my1024 + y1024, 2)/pow(sigmay1024, 3);
            break;
        case 1025:
            dJydsigmay[1025] = 1.0/sigmay1025 - 1.0*pow(-my1025 + y1025, 2)/pow(sigmay1025, 3);
            break;
        case 1026:
            dJydsigmay[1026] = 1.0/sigmay1026 - 1.0*pow(-my1026 + y1026, 2)/pow(sigmay1026, 3);
            break;
        case 1027:
            dJydsigmay[1027] = 1.0/sigmay1027 - 1.0*pow(-my1027 + y1027, 2)/pow(sigmay1027, 3);
            break;
        case 1028:
            dJydsigmay[1028] = 1.0/sigmay1028 - 1.0*pow(-my1028 + y1028, 2)/pow(sigmay1028, 3);
            break;
        case 1029:
            dJydsigmay[1029] = 1.0/sigmay1029 - 1.0*pow(-my1029 + y1029, 2)/pow(sigmay1029, 3);
            break;
        case 1030:
            dJydsigmay[1030] = 1.0/sigmay1030 - 1.0*pow(-my1030 + y1030, 2)/pow(sigmay1030, 3);
            break;
        case 1031:
            dJydsigmay[1031] = 1.0/sigmay1031 - 1.0*pow(-my1031 + y1031, 2)/pow(sigmay1031, 3);
            break;
        case 1032:
            dJydsigmay[1032] = 1.0/sigmay1032 - 1.0*pow(-my1032 + y1032, 2)/pow(sigmay1032, 3);
            break;
        case 1033:
            dJydsigmay[1033] = 1.0/sigmay1033 - 1.0*pow(-my1033 + y1033, 2)/pow(sigmay1033, 3);
            break;
        case 1034:
            dJydsigmay[1034] = 1.0/sigmay1034 - 1.0*pow(-my1034 + y1034, 2)/pow(sigmay1034, 3);
            break;
        case 1035:
            dJydsigmay[1035] = 1.0/sigmay1035 - 1.0*pow(-my1035 + y1035, 2)/pow(sigmay1035, 3);
            break;
        case 1036:
            dJydsigmay[1036] = 1.0/sigmay1036 - 1.0*pow(-my1036 + y1036, 2)/pow(sigmay1036, 3);
            break;
        case 1037:
            dJydsigmay[1037] = 1.0/sigmay1037 - 1.0*pow(-my1037 + y1037, 2)/pow(sigmay1037, 3);
            break;
        case 1038:
            dJydsigmay[1038] = 1.0/sigmay1038 - 1.0*pow(-my1038 + y1038, 2)/pow(sigmay1038, 3);
            break;
        case 1039:
            dJydsigmay[1039] = 1.0/sigmay1039 - 1.0*pow(-my1039 + y1039, 2)/pow(sigmay1039, 3);
            break;
        case 1040:
            dJydsigmay[1040] = 1.0/sigmay1040 - 1.0*pow(-my1040 + y1040, 2)/pow(sigmay1040, 3);
            break;
        case 1041:
            dJydsigmay[1041] = 1.0/sigmay1041 - 1.0*pow(-my1041 + y1041, 2)/pow(sigmay1041, 3);
            break;
        case 1042:
            dJydsigmay[1042] = 1.0/sigmay1042 - 1.0*pow(-my1042 + y1042, 2)/pow(sigmay1042, 3);
            break;
        case 1043:
            dJydsigmay[1043] = 1.0/sigmay1043 - 1.0*pow(-my1043 + y1043, 2)/pow(sigmay1043, 3);
            break;
        case 1044:
            dJydsigmay[1044] = 1.0/sigmay1044 - 1.0*pow(-my1044 + y1044, 2)/pow(sigmay1044, 3);
            break;
        case 1045:
            dJydsigmay[1045] = 1.0/sigmay1045 - 1.0*pow(-my1045 + y1045, 2)/pow(sigmay1045, 3);
            break;
        case 1046:
            dJydsigmay[1046] = 1.0/sigmay1046 - 1.0*pow(-my1046 + y1046, 2)/pow(sigmay1046, 3);
            break;
        case 1047:
            dJydsigmay[1047] = 1.0/sigmay1047 - 1.0*pow(-my1047 + y1047, 2)/pow(sigmay1047, 3);
            break;
        case 1048:
            dJydsigmay[1048] = 1.0/sigmay1048 - 1.0*pow(-my1048 + y1048, 2)/pow(sigmay1048, 3);
            break;
        case 1049:
            dJydsigmay[1049] = 1.0/sigmay1049 - 1.0*pow(-my1049 + y1049, 2)/pow(sigmay1049, 3);
            break;
        case 1050:
            dJydsigmay[1050] = 1.0/sigmay1050 - 1.0*pow(-my1050 + y1050, 2)/pow(sigmay1050, 3);
            break;
        case 1051:
            dJydsigmay[1051] = 1.0/sigmay1051 - 1.0*pow(-my1051 + y1051, 2)/pow(sigmay1051, 3);
            break;
        case 1052:
            dJydsigmay[1052] = 1.0/sigmay1052 - 1.0*pow(-my1052 + y1052, 2)/pow(sigmay1052, 3);
            break;
        case 1053:
            dJydsigmay[1053] = 1.0/sigmay1053 - 1.0*pow(-my1053 + y1053, 2)/pow(sigmay1053, 3);
            break;
        case 1054:
            dJydsigmay[1054] = 1.0/sigmay1054 - 1.0*pow(-my1054 + y1054, 2)/pow(sigmay1054, 3);
            break;
        case 1055:
            dJydsigmay[1055] = 1.0/sigmay1055 - 1.0*pow(-my1055 + y1055, 2)/pow(sigmay1055, 3);
            break;
        case 1056:
            dJydsigmay[1056] = 1.0/sigmay1056 - 1.0*pow(-my1056 + y1056, 2)/pow(sigmay1056, 3);
            break;
        case 1057:
            dJydsigmay[1057] = 1.0/sigmay1057 - 1.0*pow(-my1057 + y1057, 2)/pow(sigmay1057, 3);
            break;
        case 1058:
            dJydsigmay[1058] = 1.0/sigmay1058 - 1.0*pow(-my1058 + y1058, 2)/pow(sigmay1058, 3);
            break;
        case 1059:
            dJydsigmay[1059] = 1.0/sigmay1059 - 1.0*pow(-my1059 + y1059, 2)/pow(sigmay1059, 3);
            break;
        case 1060:
            dJydsigmay[1060] = 1.0/sigmay1060 - 1.0*pow(-my1060 + y1060, 2)/pow(sigmay1060, 3);
            break;
        case 1061:
            dJydsigmay[1061] = 1.0/sigmay1061 - 1.0*pow(-my1061 + y1061, 2)/pow(sigmay1061, 3);
            break;
        case 1062:
            dJydsigmay[1062] = 1.0/sigmay1062 - 1.0*pow(-my1062 + y1062, 2)/pow(sigmay1062, 3);
            break;
        case 1063:
            dJydsigmay[1063] = 1.0/sigmay1063 - 1.0*pow(-my1063 + y1063, 2)/pow(sigmay1063, 3);
            break;
        case 1064:
            dJydsigmay[1064] = 1.0/sigmay1064 - 1.0*pow(-my1064 + y1064, 2)/pow(sigmay1064, 3);
            break;
        case 1065:
            dJydsigmay[1065] = 1.0/sigmay1065 - 1.0*pow(-my1065 + y1065, 2)/pow(sigmay1065, 3);
            break;
        case 1066:
            dJydsigmay[1066] = 1.0/sigmay1066 - 1.0*pow(-my1066 + y1066, 2)/pow(sigmay1066, 3);
            break;
        case 1067:
            dJydsigmay[1067] = 1.0/sigmay1067 - 1.0*pow(-my1067 + y1067, 2)/pow(sigmay1067, 3);
            break;
        case 1068:
            dJydsigmay[1068] = 1.0/sigmay1068 - 1.0*pow(-my1068 + y1068, 2)/pow(sigmay1068, 3);
            break;
        case 1069:
            dJydsigmay[1069] = 1.0/sigmay1069 - 1.0*pow(-my1069 + y1069, 2)/pow(sigmay1069, 3);
            break;
        case 1070:
            dJydsigmay[1070] = 1.0/sigmay1070 - 1.0*pow(-my1070 + y1070, 2)/pow(sigmay1070, 3);
            break;
        case 1071:
            dJydsigmay[1071] = 1.0/sigmay1071 - 1.0*pow(-my1071 + y1071, 2)/pow(sigmay1071, 3);
            break;
        case 1072:
            dJydsigmay[1072] = 1.0/sigmay1072 - 1.0*pow(-my1072 + y1072, 2)/pow(sigmay1072, 3);
            break;
        case 1073:
            dJydsigmay[1073] = 1.0/sigmay1073 - 1.0*pow(-my1073 + y1073, 2)/pow(sigmay1073, 3);
            break;
        case 1074:
            dJydsigmay[1074] = 1.0/sigmay1074 - 1.0*pow(-my1074 + y1074, 2)/pow(sigmay1074, 3);
            break;
        case 1075:
            dJydsigmay[1075] = 1.0/sigmay1075 - 1.0*pow(-my1075 + y1075, 2)/pow(sigmay1075, 3);
            break;
        case 1076:
            dJydsigmay[1076] = 1.0/sigmay1076 - 1.0*pow(-my1076 + y1076, 2)/pow(sigmay1076, 3);
            break;
        case 1077:
            dJydsigmay[1077] = 1.0/sigmay1077 - 1.0*pow(-my1077 + y1077, 2)/pow(sigmay1077, 3);
            break;
        case 1078:
            dJydsigmay[1078] = 1.0/sigmay1078 - 1.0*pow(-my1078 + y1078, 2)/pow(sigmay1078, 3);
            break;
        case 1079:
            dJydsigmay[1079] = 1.0/sigmay1079 - 1.0*pow(-my1079 + y1079, 2)/pow(sigmay1079, 3);
            break;
        case 1080:
            dJydsigmay[1080] = 1.0/sigmay1080 - 1.0*pow(-my1080 + y1080, 2)/pow(sigmay1080, 3);
            break;
        case 1081:
            dJydsigmay[1081] = 1.0/sigmay1081 - 1.0*pow(-my1081 + y1081, 2)/pow(sigmay1081, 3);
            break;
        case 1082:
            dJydsigmay[1082] = 1.0/sigmay1082 - 1.0*pow(-my1082 + y1082, 2)/pow(sigmay1082, 3);
            break;
        case 1083:
            dJydsigmay[1083] = 1.0/sigmay1083 - 1.0*pow(-my1083 + y1083, 2)/pow(sigmay1083, 3);
            break;
        case 1084:
            dJydsigmay[1084] = 1.0/sigmay1084 - 1.0*pow(-my1084 + y1084, 2)/pow(sigmay1084, 3);
            break;
        case 1085:
            dJydsigmay[1085] = 1.0/sigmay1085 - 1.0*pow(-my1085 + y1085, 2)/pow(sigmay1085, 3);
            break;
        case 1086:
            dJydsigmay[1086] = 1.0/sigmay1086 - 1.0*pow(-my1086 + y1086, 2)/pow(sigmay1086, 3);
            break;
        case 1087:
            dJydsigmay[1087] = 1.0/sigmay1087 - 1.0*pow(-my1087 + y1087, 2)/pow(sigmay1087, 3);
            break;
        case 1088:
            dJydsigmay[1088] = 1.0/sigmay1088 - 1.0*pow(-my1088 + y1088, 2)/pow(sigmay1088, 3);
            break;
        case 1089:
            dJydsigmay[1089] = 1.0/sigmay1089 - 1.0*pow(-my1089 + y1089, 2)/pow(sigmay1089, 3);
            break;
        case 1090:
            dJydsigmay[1090] = 1.0/sigmay1090 - 1.0*pow(-my1090 + y1090, 2)/pow(sigmay1090, 3);
            break;
        case 1091:
            dJydsigmay[1091] = 1.0/sigmay1091 - 1.0*pow(-my1091 + y1091, 2)/pow(sigmay1091, 3);
            break;
        case 1092:
            dJydsigmay[1092] = 1.0/sigmay1092 - 1.0*pow(-my1092 + y1092, 2)/pow(sigmay1092, 3);
            break;
        case 1093:
            dJydsigmay[1093] = 1.0/sigmay1093 - 1.0*pow(-my1093 + y1093, 2)/pow(sigmay1093, 3);
            break;
        case 1094:
            dJydsigmay[1094] = 1.0/sigmay1094 - 1.0*pow(-my1094 + y1094, 2)/pow(sigmay1094, 3);
            break;
        case 1095:
            dJydsigmay[1095] = 1.0/sigmay1095 - 1.0*pow(-my1095 + y1095, 2)/pow(sigmay1095, 3);
            break;
        case 1096:
            dJydsigmay[1096] = 1.0/sigmay1096 - 1.0*pow(-my1096 + y1096, 2)/pow(sigmay1096, 3);
            break;
        case 1097:
            dJydsigmay[1097] = 1.0/sigmay1097 - 1.0*pow(-my1097 + y1097, 2)/pow(sigmay1097, 3);
            break;
        case 1098:
            dJydsigmay[1098] = 1.0/sigmay1098 - 1.0*pow(-my1098 + y1098, 2)/pow(sigmay1098, 3);
            break;
        case 1099:
            dJydsigmay[1099] = 1.0/sigmay1099 - 1.0*pow(-my1099 + y1099, 2)/pow(sigmay1099, 3);
            break;
        case 1100:
            dJydsigmay[1100] = 1.0/sigmay1100 - 1.0*pow(-my1100 + y1100, 2)/pow(sigmay1100, 3);
            break;
        case 1101:
            dJydsigmay[1101] = 1.0/sigmay1101 - 1.0*pow(-my1101 + y1101, 2)/pow(sigmay1101, 3);
            break;
        case 1102:
            dJydsigmay[1102] = 1.0/sigmay1102 - 1.0*pow(-my1102 + y1102, 2)/pow(sigmay1102, 3);
            break;
        case 1103:
            dJydsigmay[1103] = 1.0/sigmay1103 - 1.0*pow(-my1103 + y1103, 2)/pow(sigmay1103, 3);
            break;
        case 1104:
            dJydsigmay[1104] = 1.0/sigmay1104 - 1.0*pow(-my1104 + y1104, 2)/pow(sigmay1104, 3);
            break;
        case 1105:
            dJydsigmay[1105] = 1.0/sigmay1105 - 1.0*pow(-my1105 + y1105, 2)/pow(sigmay1105, 3);
            break;
        case 1106:
            dJydsigmay[1106] = 1.0/sigmay1106 - 1.0*pow(-my1106 + y1106, 2)/pow(sigmay1106, 3);
            break;
        case 1107:
            dJydsigmay[1107] = 1.0/sigmay1107 - 1.0*pow(-my1107 + y1107, 2)/pow(sigmay1107, 3);
            break;
        case 1108:
            dJydsigmay[1108] = 1.0/sigmay1108 - 1.0*pow(-my1108 + y1108, 2)/pow(sigmay1108, 3);
            break;
        case 1109:
            dJydsigmay[1109] = 1.0/sigmay1109 - 1.0*pow(-my1109 + y1109, 2)/pow(sigmay1109, 3);
            break;
        case 1110:
            dJydsigmay[1110] = 1.0/sigmay1110 - 1.0*pow(-my1110 + y1110, 2)/pow(sigmay1110, 3);
            break;
        case 1111:
            dJydsigmay[1111] = 1.0/sigmay1111 - 1.0*pow(-my1111 + y1111, 2)/pow(sigmay1111, 3);
            break;
        case 1112:
            dJydsigmay[1112] = 1.0/sigmay1112 - 1.0*pow(-my1112 + y1112, 2)/pow(sigmay1112, 3);
            break;
        case 1113:
            dJydsigmay[1113] = 1.0/sigmay1113 - 1.0*pow(-my1113 + y1113, 2)/pow(sigmay1113, 3);
            break;
        case 1114:
            dJydsigmay[1114] = 1.0/sigmay1114 - 1.0*pow(-my1114 + y1114, 2)/pow(sigmay1114, 3);
            break;
        case 1115:
            dJydsigmay[1115] = 1.0/sigmay1115 - 1.0*pow(-my1115 + y1115, 2)/pow(sigmay1115, 3);
            break;
        case 1116:
            dJydsigmay[1116] = 1.0/sigmay1116 - 1.0*pow(-my1116 + y1116, 2)/pow(sigmay1116, 3);
            break;
        case 1117:
            dJydsigmay[1117] = 1.0/sigmay1117 - 1.0*pow(-my1117 + y1117, 2)/pow(sigmay1117, 3);
            break;
        case 1118:
            dJydsigmay[1118] = 1.0/sigmay1118 - 1.0*pow(-my1118 + y1118, 2)/pow(sigmay1118, 3);
            break;
        case 1119:
            dJydsigmay[1119] = 1.0/sigmay1119 - 1.0*pow(-my1119 + y1119, 2)/pow(sigmay1119, 3);
            break;
        case 1120:
            dJydsigmay[1120] = 1.0/sigmay1120 - 1.0*pow(-my1120 + y1120, 2)/pow(sigmay1120, 3);
            break;
        case 1121:
            dJydsigmay[1121] = 1.0/sigmay1121 - 1.0*pow(-my1121 + y1121, 2)/pow(sigmay1121, 3);
            break;
        case 1122:
            dJydsigmay[1122] = 1.0/sigmay1122 - 1.0*pow(-my1122 + y1122, 2)/pow(sigmay1122, 3);
            break;
        case 1123:
            dJydsigmay[1123] = 1.0/sigmay1123 - 1.0*pow(-my1123 + y1123, 2)/pow(sigmay1123, 3);
            break;
        case 1124:
            dJydsigmay[1124] = 1.0/sigmay1124 - 1.0*pow(-my1124 + y1124, 2)/pow(sigmay1124, 3);
            break;
        case 1125:
            dJydsigmay[1125] = 1.0/sigmay1125 - 1.0*pow(-my1125 + y1125, 2)/pow(sigmay1125, 3);
            break;
        case 1126:
            dJydsigmay[1126] = 1.0/sigmay1126 - 1.0*pow(-my1126 + y1126, 2)/pow(sigmay1126, 3);
            break;
        case 1127:
            dJydsigmay[1127] = 1.0/sigmay1127 - 1.0*pow(-my1127 + y1127, 2)/pow(sigmay1127, 3);
            break;
        case 1128:
            dJydsigmay[1128] = 1.0/sigmay1128 - 1.0*pow(-my1128 + y1128, 2)/pow(sigmay1128, 3);
            break;
        case 1129:
            dJydsigmay[1129] = 1.0/sigmay1129 - 1.0*pow(-my1129 + y1129, 2)/pow(sigmay1129, 3);
            break;
        case 1130:
            dJydsigmay[1130] = 1.0/sigmay1130 - 1.0*pow(-my1130 + y1130, 2)/pow(sigmay1130, 3);
            break;
        case 1131:
            dJydsigmay[1131] = 1.0/sigmay1131 - 1.0*pow(-my1131 + y1131, 2)/pow(sigmay1131, 3);
            break;
        case 1132:
            dJydsigmay[1132] = 1.0/sigmay1132 - 1.0*pow(-my1132 + y1132, 2)/pow(sigmay1132, 3);
            break;
        case 1133:
            dJydsigmay[1133] = 1.0/sigmay1133 - 1.0*pow(-my1133 + y1133, 2)/pow(sigmay1133, 3);
            break;
        case 1134:
            dJydsigmay[1134] = 1.0/sigmay1134 - 1.0*pow(-my1134 + y1134, 2)/pow(sigmay1134, 3);
            break;
        case 1135:
            dJydsigmay[1135] = 1.0/sigmay1135 - 1.0*pow(-my1135 + y1135, 2)/pow(sigmay1135, 3);
            break;
        case 1136:
            dJydsigmay[1136] = 1.0/sigmay1136 - 1.0*pow(-my1136 + y1136, 2)/pow(sigmay1136, 3);
            break;
        case 1137:
            dJydsigmay[1137] = 1.0/sigmay1137 - 1.0*pow(-my1137 + y1137, 2)/pow(sigmay1137, 3);
            break;
        case 1138:
            dJydsigmay[1138] = 1.0/sigmay1138 - 1.0*pow(-my1138 + y1138, 2)/pow(sigmay1138, 3);
            break;
        case 1139:
            dJydsigmay[1139] = 1.0/sigmay1139 - 1.0*pow(-my1139 + y1139, 2)/pow(sigmay1139, 3);
            break;
        case 1140:
            dJydsigmay[1140] = 1.0/sigmay1140 - 1.0*pow(-my1140 + y1140, 2)/pow(sigmay1140, 3);
            break;
        case 1141:
            dJydsigmay[1141] = 1.0/sigmay1141 - 1.0*pow(-my1141 + y1141, 2)/pow(sigmay1141, 3);
            break;
        case 1142:
            dJydsigmay[1142] = 1.0/sigmay1142 - 1.0*pow(-my1142 + y1142, 2)/pow(sigmay1142, 3);
            break;
        case 1143:
            dJydsigmay[1143] = 1.0/sigmay1143 - 1.0*pow(-my1143 + y1143, 2)/pow(sigmay1143, 3);
            break;
        case 1144:
            dJydsigmay[1144] = 1.0/sigmay1144 - 1.0*pow(-my1144 + y1144, 2)/pow(sigmay1144, 3);
            break;
        case 1145:
            dJydsigmay[1145] = 1.0/sigmay1145 - 1.0*pow(-my1145 + y1145, 2)/pow(sigmay1145, 3);
            break;
        case 1146:
            dJydsigmay[1146] = 1.0/sigmay1146 - 1.0*pow(-my1146 + y1146, 2)/pow(sigmay1146, 3);
            break;
        case 1147:
            dJydsigmay[1147] = 1.0/sigmay1147 - 1.0*pow(-my1147 + y1147, 2)/pow(sigmay1147, 3);
            break;
        case 1148:
            dJydsigmay[1148] = 1.0/sigmay1148 - 1.0*pow(-my1148 + y1148, 2)/pow(sigmay1148, 3);
            break;
        case 1149:
            dJydsigmay[1149] = 1.0/sigmay1149 - 1.0*pow(-my1149 + y1149, 2)/pow(sigmay1149, 3);
            break;
        case 1150:
            dJydsigmay[1150] = 1.0/sigmay1150 - 1.0*pow(-my1150 + y1150, 2)/pow(sigmay1150, 3);
            break;
        case 1151:
            dJydsigmay[1151] = 1.0/sigmay1151 - 1.0*pow(-my1151 + y1151, 2)/pow(sigmay1151, 3);
            break;
        case 1152:
            dJydsigmay[1152] = 1.0/sigmay1152 - 1.0*pow(-my1152 + y1152, 2)/pow(sigmay1152, 3);
            break;
        case 1153:
            dJydsigmay[1153] = 1.0/sigmay1153 - 1.0*pow(-my1153 + y1153, 2)/pow(sigmay1153, 3);
            break;
        case 1154:
            dJydsigmay[1154] = 1.0/sigmay1154 - 1.0*pow(-my1154 + y1154, 2)/pow(sigmay1154, 3);
            break;
        case 1155:
            dJydsigmay[1155] = 1.0/sigmay1155 - 1.0*pow(-my1155 + y1155, 2)/pow(sigmay1155, 3);
            break;
        case 1156:
            dJydsigmay[1156] = 1.0/sigmay1156 - 1.0*pow(-my1156 + y1156, 2)/pow(sigmay1156, 3);
            break;
        case 1157:
            dJydsigmay[1157] = 1.0/sigmay1157 - 1.0*pow(-my1157 + y1157, 2)/pow(sigmay1157, 3);
            break;
        case 1158:
            dJydsigmay[1158] = 1.0/sigmay1158 - 1.0*pow(-my1158 + y1158, 2)/pow(sigmay1158, 3);
            break;
        case 1159:
            dJydsigmay[1159] = 1.0/sigmay1159 - 1.0*pow(-my1159 + y1159, 2)/pow(sigmay1159, 3);
            break;
        case 1160:
            dJydsigmay[1160] = 1.0/sigmay1160 - 1.0*pow(-my1160 + y1160, 2)/pow(sigmay1160, 3);
            break;
        case 1161:
            dJydsigmay[1161] = 1.0/sigmay1161 - 1.0*pow(-my1161 + y1161, 2)/pow(sigmay1161, 3);
            break;
        case 1162:
            dJydsigmay[1162] = 1.0/sigmay1162 - 1.0*pow(-my1162 + y1162, 2)/pow(sigmay1162, 3);
            break;
        case 1163:
            dJydsigmay[1163] = 1.0/sigmay1163 - 1.0*pow(-my1163 + y1163, 2)/pow(sigmay1163, 3);
            break;
        case 1164:
            dJydsigmay[1164] = 1.0/sigmay1164 - 1.0*pow(-my1164 + y1164, 2)/pow(sigmay1164, 3);
            break;
        case 1165:
            dJydsigmay[1165] = 1.0/sigmay1165 - 1.0*pow(-my1165 + y1165, 2)/pow(sigmay1165, 3);
            break;
        case 1166:
            dJydsigmay[1166] = 1.0/sigmay1166 - 1.0*pow(-my1166 + y1166, 2)/pow(sigmay1166, 3);
            break;
        case 1167:
            dJydsigmay[1167] = 1.0/sigmay1167 - 1.0*pow(-my1167 + y1167, 2)/pow(sigmay1167, 3);
            break;
        case 1168:
            dJydsigmay[1168] = 1.0/sigmay1168 - 1.0*pow(-my1168 + y1168, 2)/pow(sigmay1168, 3);
            break;
        case 1169:
            dJydsigmay[1169] = 1.0/sigmay1169 - 1.0*pow(-my1169 + y1169, 2)/pow(sigmay1169, 3);
            break;
        case 1170:
            dJydsigmay[1170] = 1.0/sigmay1170 - 1.0*pow(-my1170 + y1170, 2)/pow(sigmay1170, 3);
            break;
        case 1171:
            dJydsigmay[1171] = 1.0/sigmay1171 - 1.0*pow(-my1171 + y1171, 2)/pow(sigmay1171, 3);
            break;
        case 1172:
            dJydsigmay[1172] = 1.0/sigmay1172 - 1.0*pow(-my1172 + y1172, 2)/pow(sigmay1172, 3);
            break;
        case 1173:
            dJydsigmay[1173] = 1.0/sigmay1173 - 1.0*pow(-my1173 + y1173, 2)/pow(sigmay1173, 3);
            break;
        case 1174:
            dJydsigmay[1174] = 1.0/sigmay1174 - 1.0*pow(-my1174 + y1174, 2)/pow(sigmay1174, 3);
            break;
        case 1175:
            dJydsigmay[1175] = 1.0/sigmay1175 - 1.0*pow(-my1175 + y1175, 2)/pow(sigmay1175, 3);
            break;
        case 1176:
            dJydsigmay[1176] = 1.0/sigmay1176 - 1.0*pow(-my1176 + y1176, 2)/pow(sigmay1176, 3);
            break;
        case 1177:
            dJydsigmay[1177] = 1.0/sigmay1177 - 1.0*pow(-my1177 + y1177, 2)/pow(sigmay1177, 3);
            break;
        case 1178:
            dJydsigmay[1178] = 1.0/sigmay1178 - 1.0*pow(-my1178 + y1178, 2)/pow(sigmay1178, 3);
            break;
        case 1179:
            dJydsigmay[1179] = 1.0/sigmay1179 - 1.0*pow(-my1179 + y1179, 2)/pow(sigmay1179, 3);
            break;
        case 1180:
            dJydsigmay[1180] = 1.0/sigmay1180 - 1.0*pow(-my1180 + y1180, 2)/pow(sigmay1180, 3);
            break;
        case 1181:
            dJydsigmay[1181] = 1.0/sigmay1181 - 1.0*pow(-my1181 + y1181, 2)/pow(sigmay1181, 3);
            break;
        case 1182:
            dJydsigmay[1182] = 1.0/sigmay1182 - 1.0*pow(-my1182 + y1182, 2)/pow(sigmay1182, 3);
            break;
        case 1183:
            dJydsigmay[1183] = 1.0/sigmay1183 - 1.0*pow(-my1183 + y1183, 2)/pow(sigmay1183, 3);
            break;
        case 1184:
            dJydsigmay[1184] = 1.0/sigmay1184 - 1.0*pow(-my1184 + y1184, 2)/pow(sigmay1184, 3);
            break;
        case 1185:
            dJydsigmay[1185] = 1.0/sigmay1185 - 1.0*pow(-my1185 + y1185, 2)/pow(sigmay1185, 3);
            break;
        case 1186:
            dJydsigmay[1186] = 1.0/sigmay1186 - 1.0*pow(-my1186 + y1186, 2)/pow(sigmay1186, 3);
            break;
        case 1187:
            dJydsigmay[1187] = 1.0/sigmay1187 - 1.0*pow(-my1187 + y1187, 2)/pow(sigmay1187, 3);
            break;
        case 1188:
            dJydsigmay[1188] = 1.0/sigmay1188 - 1.0*pow(-my1188 + y1188, 2)/pow(sigmay1188, 3);
            break;
        case 1189:
            dJydsigmay[1189] = 1.0/sigmay1189 - 1.0*pow(-my1189 + y1189, 2)/pow(sigmay1189, 3);
            break;
        case 1190:
            dJydsigmay[1190] = 1.0/sigmay1190 - 1.0*pow(-my1190 + y1190, 2)/pow(sigmay1190, 3);
            break;
        case 1191:
            dJydsigmay[1191] = 1.0/sigmay1191 - 1.0*pow(-my1191 + y1191, 2)/pow(sigmay1191, 3);
            break;
        case 1192:
            dJydsigmay[1192] = 1.0/sigmay1192 - 1.0*pow(-my1192 + y1192, 2)/pow(sigmay1192, 3);
            break;
        case 1193:
            dJydsigmay[1193] = 1.0/sigmay1193 - 1.0*pow(-my1193 + y1193, 2)/pow(sigmay1193, 3);
            break;
        case 1194:
            dJydsigmay[1194] = 1.0/sigmay1194 - 1.0*pow(-my1194 + y1194, 2)/pow(sigmay1194, 3);
            break;
        case 1195:
            dJydsigmay[1195] = 1.0/sigmay1195 - 1.0*pow(-my1195 + y1195, 2)/pow(sigmay1195, 3);
            break;
        case 1196:
            dJydsigmay[1196] = 1.0/sigmay1196 - 1.0*pow(-my1196 + y1196, 2)/pow(sigmay1196, 3);
            break;
        case 1197:
            dJydsigmay[1197] = 1.0/sigmay1197 - 1.0*pow(-my1197 + y1197, 2)/pow(sigmay1197, 3);
            break;
        case 1198:
            dJydsigmay[1198] = 1.0/sigmay1198 - 1.0*pow(-my1198 + y1198, 2)/pow(sigmay1198, 3);
            break;
        case 1199:
            dJydsigmay[1199] = 1.0/sigmay1199 - 1.0*pow(-my1199 + y1199, 2)/pow(sigmay1199, 3);
            break;
        case 1200:
            dJydsigmay[1200] = 1.0/sigmay1200 - 1.0*pow(-my1200 + y1200, 2)/pow(sigmay1200, 3);
            break;
        case 1201:
            dJydsigmay[1201] = 1.0/sigmay1201 - 1.0*pow(-my1201 + y1201, 2)/pow(sigmay1201, 3);
            break;
        case 1202:
            dJydsigmay[1202] = 1.0/sigmay1202 - 1.0*pow(-my1202 + y1202, 2)/pow(sigmay1202, 3);
            break;
        case 1203:
            dJydsigmay[1203] = 1.0/sigmay1203 - 1.0*pow(-my1203 + y1203, 2)/pow(sigmay1203, 3);
            break;
        case 1204:
            dJydsigmay[1204] = 1.0/sigmay1204 - 1.0*pow(-my1204 + y1204, 2)/pow(sigmay1204, 3);
            break;
        case 1205:
            dJydsigmay[1205] = 1.0/sigmay1205 - 1.0*pow(-my1205 + y1205, 2)/pow(sigmay1205, 3);
            break;
        case 1206:
            dJydsigmay[1206] = 1.0/sigmay1206 - 1.0*pow(-my1206 + y1206, 2)/pow(sigmay1206, 3);
            break;
        case 1207:
            dJydsigmay[1207] = 1.0/sigmay1207 - 1.0*pow(-my1207 + y1207, 2)/pow(sigmay1207, 3);
            break;
        case 1208:
            dJydsigmay[1208] = 1.0/sigmay1208 - 1.0*pow(-my1208 + y1208, 2)/pow(sigmay1208, 3);
            break;
        case 1209:
            dJydsigmay[1209] = 1.0/sigmay1209 - 1.0*pow(-my1209 + y1209, 2)/pow(sigmay1209, 3);
            break;
        case 1210:
            dJydsigmay[1210] = 1.0/sigmay1210 - 1.0*pow(-my1210 + y1210, 2)/pow(sigmay1210, 3);
            break;
        case 1211:
            dJydsigmay[1211] = 1.0/sigmay1211 - 1.0*pow(-my1211 + y1211, 2)/pow(sigmay1211, 3);
            break;
        case 1212:
            dJydsigmay[1212] = 1.0/sigmay1212 - 1.0*pow(-my1212 + y1212, 2)/pow(sigmay1212, 3);
            break;
        case 1213:
            dJydsigmay[1213] = 1.0/sigmay1213 - 1.0*pow(-my1213 + y1213, 2)/pow(sigmay1213, 3);
            break;
        case 1214:
            dJydsigmay[1214] = 1.0/sigmay1214 - 1.0*pow(-my1214 + y1214, 2)/pow(sigmay1214, 3);
            break;
        case 1215:
            dJydsigmay[1215] = 1.0/sigmay1215 - 1.0*pow(-my1215 + y1215, 2)/pow(sigmay1215, 3);
            break;
        case 1216:
            dJydsigmay[1216] = 1.0/sigmay1216 - 1.0*pow(-my1216 + y1216, 2)/pow(sigmay1216, 3);
            break;
        case 1217:
            dJydsigmay[1217] = 1.0/sigmay1217 - 1.0*pow(-my1217 + y1217, 2)/pow(sigmay1217, 3);
            break;
        case 1218:
            dJydsigmay[1218] = 1.0/sigmay1218 - 1.0*pow(-my1218 + y1218, 2)/pow(sigmay1218, 3);
            break;
        case 1219:
            dJydsigmay[1219] = 1.0/sigmay1219 - 1.0*pow(-my1219 + y1219, 2)/pow(sigmay1219, 3);
            break;
        case 1220:
            dJydsigmay[1220] = 1.0/sigmay1220 - 1.0*pow(-my1220 + y1220, 2)/pow(sigmay1220, 3);
            break;
        case 1221:
            dJydsigmay[1221] = 1.0/sigmay1221 - 1.0*pow(-my1221 + y1221, 2)/pow(sigmay1221, 3);
            break;
        case 1222:
            dJydsigmay[1222] = 1.0/sigmay1222 - 1.0*pow(-my1222 + y1222, 2)/pow(sigmay1222, 3);
            break;
        case 1223:
            dJydsigmay[1223] = 1.0/sigmay1223 - 1.0*pow(-my1223 + y1223, 2)/pow(sigmay1223, 3);
            break;
        case 1224:
            dJydsigmay[1224] = 1.0/sigmay1224 - 1.0*pow(-my1224 + y1224, 2)/pow(sigmay1224, 3);
            break;
        case 1225:
            dJydsigmay[1225] = 1.0/sigmay1225 - 1.0*pow(-my1225 + y1225, 2)/pow(sigmay1225, 3);
            break;
        case 1226:
            dJydsigmay[1226] = 1.0/sigmay1226 - 1.0*pow(-my1226 + y1226, 2)/pow(sigmay1226, 3);
            break;
        case 1227:
            dJydsigmay[1227] = 1.0/sigmay1227 - 1.0*pow(-my1227 + y1227, 2)/pow(sigmay1227, 3);
            break;
        case 1228:
            dJydsigmay[1228] = 1.0/sigmay1228 - 1.0*pow(-my1228 + y1228, 2)/pow(sigmay1228, 3);
            break;
        case 1229:
            dJydsigmay[1229] = 1.0/sigmay1229 - 1.0*pow(-my1229 + y1229, 2)/pow(sigmay1229, 3);
            break;
        case 1230:
            dJydsigmay[1230] = 1.0/sigmay1230 - 1.0*pow(-my1230 + y1230, 2)/pow(sigmay1230, 3);
            break;
        case 1231:
            dJydsigmay[1231] = 1.0/sigmay1231 - 1.0*pow(-my1231 + y1231, 2)/pow(sigmay1231, 3);
            break;
        case 1232:
            dJydsigmay[1232] = 1.0/sigmay1232 - 1.0*pow(-my1232 + y1232, 2)/pow(sigmay1232, 3);
            break;
        case 1233:
            dJydsigmay[1233] = 1.0/sigmay1233 - 1.0*pow(-my1233 + y1233, 2)/pow(sigmay1233, 3);
            break;
        case 1234:
            dJydsigmay[1234] = 1.0/sigmay1234 - 1.0*pow(-my1234 + y1234, 2)/pow(sigmay1234, 3);
            break;
        case 1235:
            dJydsigmay[1235] = 1.0/sigmay1235 - 1.0*pow(-my1235 + y1235, 2)/pow(sigmay1235, 3);
            break;
        case 1236:
            dJydsigmay[1236] = 1.0/sigmay1236 - 1.0*pow(-my1236 + y1236, 2)/pow(sigmay1236, 3);
            break;
        case 1237:
            dJydsigmay[1237] = 1.0/sigmay1237 - 1.0*pow(-my1237 + y1237, 2)/pow(sigmay1237, 3);
            break;
        case 1238:
            dJydsigmay[1238] = 1.0/sigmay1238 - 1.0*pow(-my1238 + y1238, 2)/pow(sigmay1238, 3);
            break;
        case 1239:
            dJydsigmay[1239] = 1.0/sigmay1239 - 1.0*pow(-my1239 + y1239, 2)/pow(sigmay1239, 3);
            break;
        case 1240:
            dJydsigmay[1240] = 1.0/sigmay1240 - 1.0*pow(-my1240 + y1240, 2)/pow(sigmay1240, 3);
            break;
        case 1241:
            dJydsigmay[1241] = 1.0/sigmay1241 - 1.0*pow(-my1241 + y1241, 2)/pow(sigmay1241, 3);
            break;
        case 1242:
            dJydsigmay[1242] = 1.0/sigmay1242 - 1.0*pow(-my1242 + y1242, 2)/pow(sigmay1242, 3);
            break;
        case 1243:
            dJydsigmay[1243] = 1.0/sigmay1243 - 1.0*pow(-my1243 + y1243, 2)/pow(sigmay1243, 3);
            break;
        case 1244:
            dJydsigmay[1244] = 1.0/sigmay1244 - 1.0*pow(-my1244 + y1244, 2)/pow(sigmay1244, 3);
            break;
        case 1245:
            dJydsigmay[1245] = 1.0/sigmay1245 - 1.0*pow(-my1245 + y1245, 2)/pow(sigmay1245, 3);
            break;
        case 1246:
            dJydsigmay[1246] = 1.0/sigmay1246 - 1.0*pow(-my1246 + y1246, 2)/pow(sigmay1246, 3);
            break;
        case 1247:
            dJydsigmay[1247] = 1.0/sigmay1247 - 1.0*pow(-my1247 + y1247, 2)/pow(sigmay1247, 3);
            break;
        case 1248:
            dJydsigmay[1248] = 1.0/sigmay1248 - 1.0*pow(-my1248 + y1248, 2)/pow(sigmay1248, 3);
            break;
        case 1249:
            dJydsigmay[1249] = 1.0/sigmay1249 - 1.0*pow(-my1249 + y1249, 2)/pow(sigmay1249, 3);
            break;
        case 1250:
            dJydsigmay[1250] = 1.0/sigmay1250 - 1.0*pow(-my1250 + y1250, 2)/pow(sigmay1250, 3);
            break;
        case 1251:
            dJydsigmay[1251] = 1.0/sigmay1251 - 1.0*pow(-my1251 + y1251, 2)/pow(sigmay1251, 3);
            break;
        case 1252:
            dJydsigmay[1252] = 1.0/sigmay1252 - 1.0*pow(-my1252 + y1252, 2)/pow(sigmay1252, 3);
            break;
        case 1253:
            dJydsigmay[1253] = 1.0/sigmay1253 - 1.0*pow(-my1253 + y1253, 2)/pow(sigmay1253, 3);
            break;
        case 1254:
            dJydsigmay[1254] = 1.0/sigmay1254 - 1.0*pow(-my1254 + y1254, 2)/pow(sigmay1254, 3);
            break;
        case 1255:
            dJydsigmay[1255] = 1.0/sigmay1255 - 1.0*pow(-my1255 + y1255, 2)/pow(sigmay1255, 3);
            break;
        case 1256:
            dJydsigmay[1256] = 1.0/sigmay1256 - 1.0*pow(-my1256 + y1256, 2)/pow(sigmay1256, 3);
            break;
        case 1257:
            dJydsigmay[1257] = 1.0/sigmay1257 - 1.0*pow(-my1257 + y1257, 2)/pow(sigmay1257, 3);
            break;
        case 1258:
            dJydsigmay[1258] = 1.0/sigmay1258 - 1.0*pow(-my1258 + y1258, 2)/pow(sigmay1258, 3);
            break;
        case 1259:
            dJydsigmay[1259] = 1.0/sigmay1259 - 1.0*pow(-my1259 + y1259, 2)/pow(sigmay1259, 3);
            break;
        case 1260:
            dJydsigmay[1260] = 1.0/sigmay1260 - 1.0*pow(-my1260 + y1260, 2)/pow(sigmay1260, 3);
            break;
        case 1261:
            dJydsigmay[1261] = 1.0/sigmay1261 - 1.0*pow(-my1261 + y1261, 2)/pow(sigmay1261, 3);
            break;
        case 1262:
            dJydsigmay[1262] = 1.0/sigmay1262 - 1.0*pow(-my1262 + y1262, 2)/pow(sigmay1262, 3);
            break;
        case 1263:
            dJydsigmay[1263] = 1.0/sigmay1263 - 1.0*pow(-my1263 + y1263, 2)/pow(sigmay1263, 3);
            break;
        case 1264:
            dJydsigmay[1264] = 1.0/sigmay1264 - 1.0*pow(-my1264 + y1264, 2)/pow(sigmay1264, 3);
            break;
        case 1265:
            dJydsigmay[1265] = 1.0/sigmay1265 - 1.0*pow(-my1265 + y1265, 2)/pow(sigmay1265, 3);
            break;
        case 1266:
            dJydsigmay[1266] = 1.0/sigmay1266 - 1.0*pow(-my1266 + y1266, 2)/pow(sigmay1266, 3);
            break;
        case 1267:
            dJydsigmay[1267] = 1.0/sigmay1267 - 1.0*pow(-my1267 + y1267, 2)/pow(sigmay1267, 3);
            break;
        case 1268:
            dJydsigmay[1268] = 1.0/sigmay1268 - 1.0*pow(-my1268 + y1268, 2)/pow(sigmay1268, 3);
            break;
        case 1269:
            dJydsigmay[1269] = 1.0/sigmay1269 - 1.0*pow(-my1269 + y1269, 2)/pow(sigmay1269, 3);
            break;
        case 1270:
            dJydsigmay[1270] = 1.0/sigmay1270 - 1.0*pow(-my1270 + y1270, 2)/pow(sigmay1270, 3);
            break;
        case 1271:
            dJydsigmay[1271] = 1.0/sigmay1271 - 1.0*pow(-my1271 + y1271, 2)/pow(sigmay1271, 3);
            break;
        case 1272:
            dJydsigmay[1272] = 1.0/sigmay1272 - 1.0*pow(-my1272 + y1272, 2)/pow(sigmay1272, 3);
            break;
        case 1273:
            dJydsigmay[1273] = 1.0/sigmay1273 - 1.0*pow(-my1273 + y1273, 2)/pow(sigmay1273, 3);
            break;
        case 1274:
            dJydsigmay[1274] = 1.0/sigmay1274 - 1.0*pow(-my1274 + y1274, 2)/pow(sigmay1274, 3);
            break;
        case 1275:
            dJydsigmay[1275] = 1.0/sigmay1275 - 1.0*pow(-my1275 + y1275, 2)/pow(sigmay1275, 3);
            break;
        case 1276:
            dJydsigmay[1276] = 1.0/sigmay1276 - 1.0*pow(-my1276 + y1276, 2)/pow(sigmay1276, 3);
            break;
        case 1277:
            dJydsigmay[1277] = 1.0/sigmay1277 - 1.0*pow(-my1277 + y1277, 2)/pow(sigmay1277, 3);
            break;
        case 1278:
            dJydsigmay[1278] = 1.0/sigmay1278 - 1.0*pow(-my1278 + y1278, 2)/pow(sigmay1278, 3);
            break;
        case 1279:
            dJydsigmay[1279] = 1.0/sigmay1279 - 1.0*pow(-my1279 + y1279, 2)/pow(sigmay1279, 3);
            break;
        case 1280:
            dJydsigmay[1280] = 1.0/sigmay1280 - 1.0*pow(-my1280 + y1280, 2)/pow(sigmay1280, 3);
            break;
        case 1281:
            dJydsigmay[1281] = 1.0/sigmay1281 - 1.0*pow(-my1281 + y1281, 2)/pow(sigmay1281, 3);
            break;
        case 1282:
            dJydsigmay[1282] = 1.0/sigmay1282 - 1.0*pow(-my1282 + y1282, 2)/pow(sigmay1282, 3);
            break;
        case 1283:
            dJydsigmay[1283] = 1.0/sigmay1283 - 1.0*pow(-my1283 + y1283, 2)/pow(sigmay1283, 3);
            break;
        case 1284:
            dJydsigmay[1284] = 1.0/sigmay1284 - 1.0*pow(-my1284 + y1284, 2)/pow(sigmay1284, 3);
            break;
        case 1285:
            dJydsigmay[1285] = 1.0/sigmay1285 - 1.0*pow(-my1285 + y1285, 2)/pow(sigmay1285, 3);
            break;
        case 1286:
            dJydsigmay[1286] = 1.0/sigmay1286 - 1.0*pow(-my1286 + y1286, 2)/pow(sigmay1286, 3);
            break;
        case 1287:
            dJydsigmay[1287] = 1.0/sigmay1287 - 1.0*pow(-my1287 + y1287, 2)/pow(sigmay1287, 3);
            break;
        case 1288:
            dJydsigmay[1288] = 1.0/sigmay1288 - 1.0*pow(-my1288 + y1288, 2)/pow(sigmay1288, 3);
            break;
        case 1289:
            dJydsigmay[1289] = 1.0/sigmay1289 - 1.0*pow(-my1289 + y1289, 2)/pow(sigmay1289, 3);
            break;
        case 1290:
            dJydsigmay[1290] = 1.0/sigmay1290 - 1.0*pow(-my1290 + y1290, 2)/pow(sigmay1290, 3);
            break;
        case 1291:
            dJydsigmay[1291] = 1.0/sigmay1291 - 1.0*pow(-my1291 + y1291, 2)/pow(sigmay1291, 3);
            break;
        case 1292:
            dJydsigmay[1292] = 1.0/sigmay1292 - 1.0*pow(-my1292 + y1292, 2)/pow(sigmay1292, 3);
            break;
        case 1293:
            dJydsigmay[1293] = 1.0/sigmay1293 - 1.0*pow(-my1293 + y1293, 2)/pow(sigmay1293, 3);
            break;
        case 1294:
            dJydsigmay[1294] = 1.0/sigmay1294 - 1.0*pow(-my1294 + y1294, 2)/pow(sigmay1294, 3);
            break;
        case 1295:
            dJydsigmay[1295] = 1.0/sigmay1295 - 1.0*pow(-my1295 + y1295, 2)/pow(sigmay1295, 3);
            break;
        case 1296:
            dJydsigmay[1296] = 1.0/sigmay1296 - 1.0*pow(-my1296 + y1296, 2)/pow(sigmay1296, 3);
            break;
        case 1297:
            dJydsigmay[1297] = 1.0/sigmay1297 - 1.0*pow(-my1297 + y1297, 2)/pow(sigmay1297, 3);
            break;
        case 1298:
            dJydsigmay[1298] = 1.0/sigmay1298 - 1.0*pow(-my1298 + y1298, 2)/pow(sigmay1298, 3);
            break;
        case 1299:
            dJydsigmay[1299] = 1.0/sigmay1299 - 1.0*pow(-my1299 + y1299, 2)/pow(sigmay1299, 3);
            break;
        case 1300:
            dJydsigmay[1300] = 1.0/sigmay1300 - 1.0*pow(-my1300 + y1300, 2)/pow(sigmay1300, 3);
            break;
        case 1301:
            dJydsigmay[1301] = 1.0/sigmay1301 - 1.0*pow(-my1301 + y1301, 2)/pow(sigmay1301, 3);
            break;
        case 1302:
            dJydsigmay[1302] = 1.0/sigmay1302 - 1.0*pow(-my1302 + y1302, 2)/pow(sigmay1302, 3);
            break;
        case 1303:
            dJydsigmay[1303] = 1.0/sigmay1303 - 1.0*pow(-my1303 + y1303, 2)/pow(sigmay1303, 3);
            break;
        case 1304:
            dJydsigmay[1304] = 1.0/sigmay1304 - 1.0*pow(-my1304 + y1304, 2)/pow(sigmay1304, 3);
            break;
        case 1305:
            dJydsigmay[1305] = 1.0/sigmay1305 - 1.0*pow(-my1305 + y1305, 2)/pow(sigmay1305, 3);
            break;
        case 1306:
            dJydsigmay[1306] = 1.0/sigmay1306 - 1.0*pow(-my1306 + y1306, 2)/pow(sigmay1306, 3);
            break;
        case 1307:
            dJydsigmay[1307] = 1.0/sigmay1307 - 1.0*pow(-my1307 + y1307, 2)/pow(sigmay1307, 3);
            break;
        case 1308:
            dJydsigmay[1308] = 1.0/sigmay1308 - 1.0*pow(-my1308 + y1308, 2)/pow(sigmay1308, 3);
            break;
        case 1309:
            dJydsigmay[1309] = 1.0/sigmay1309 - 1.0*pow(-my1309 + y1309, 2)/pow(sigmay1309, 3);
            break;
        case 1310:
            dJydsigmay[1310] = 1.0/sigmay1310 - 1.0*pow(-my1310 + y1310, 2)/pow(sigmay1310, 3);
            break;
        case 1311:
            dJydsigmay[1311] = 1.0/sigmay1311 - 1.0*pow(-my1311 + y1311, 2)/pow(sigmay1311, 3);
            break;
        case 1312:
            dJydsigmay[1312] = 1.0/sigmay1312 - 1.0*pow(-my1312 + y1312, 2)/pow(sigmay1312, 3);
            break;
        case 1313:
            dJydsigmay[1313] = 1.0/sigmay1313 - 1.0*pow(-my1313 + y1313, 2)/pow(sigmay1313, 3);
            break;
        case 1314:
            dJydsigmay[1314] = 1.0/sigmay1314 - 1.0*pow(-my1314 + y1314, 2)/pow(sigmay1314, 3);
            break;
        case 1315:
            dJydsigmay[1315] = 1.0/sigmay1315 - 1.0*pow(-my1315 + y1315, 2)/pow(sigmay1315, 3);
            break;
        case 1316:
            dJydsigmay[1316] = 1.0/sigmay1316 - 1.0*pow(-my1316 + y1316, 2)/pow(sigmay1316, 3);
            break;
        case 1317:
            dJydsigmay[1317] = 1.0/sigmay1317 - 1.0*pow(-my1317 + y1317, 2)/pow(sigmay1317, 3);
            break;
        case 1318:
            dJydsigmay[1318] = 1.0/sigmay1318 - 1.0*pow(-my1318 + y1318, 2)/pow(sigmay1318, 3);
            break;
        case 1319:
            dJydsigmay[1319] = 1.0/sigmay1319 - 1.0*pow(-my1319 + y1319, 2)/pow(sigmay1319, 3);
            break;
        case 1320:
            dJydsigmay[1320] = 1.0/sigmay1320 - 1.0*pow(-my1320 + y1320, 2)/pow(sigmay1320, 3);
            break;
        case 1321:
            dJydsigmay[1321] = 1.0/sigmay1321 - 1.0*pow(-my1321 + y1321, 2)/pow(sigmay1321, 3);
            break;
        case 1322:
            dJydsigmay[1322] = 1.0/sigmay1322 - 1.0*pow(-my1322 + y1322, 2)/pow(sigmay1322, 3);
            break;
        case 1323:
            dJydsigmay[1323] = 1.0/sigmay1323 - 1.0*pow(-my1323 + y1323, 2)/pow(sigmay1323, 3);
            break;
        case 1324:
            dJydsigmay[1324] = 1.0/sigmay1324 - 1.0*pow(-my1324 + y1324, 2)/pow(sigmay1324, 3);
            break;
        case 1325:
            dJydsigmay[1325] = 1.0/sigmay1325 - 1.0*pow(-my1325 + y1325, 2)/pow(sigmay1325, 3);
            break;
        case 1326:
            dJydsigmay[1326] = 1.0/sigmay1326 - 1.0*pow(-my1326 + y1326, 2)/pow(sigmay1326, 3);
            break;
        case 1327:
            dJydsigmay[1327] = 1.0/sigmay1327 - 1.0*pow(-my1327 + y1327, 2)/pow(sigmay1327, 3);
            break;
        case 1328:
            dJydsigmay[1328] = 1.0/sigmay1328 - 1.0*pow(-my1328 + y1328, 2)/pow(sigmay1328, 3);
            break;
        case 1329:
            dJydsigmay[1329] = 1.0/sigmay1329 - 1.0*pow(-my1329 + y1329, 2)/pow(sigmay1329, 3);
            break;
        case 1330:
            dJydsigmay[1330] = 1.0/sigmay1330 - 1.0*pow(-my1330 + y1330, 2)/pow(sigmay1330, 3);
            break;
        case 1331:
            dJydsigmay[1331] = 1.0/sigmay1331 - 1.0*pow(-my1331 + y1331, 2)/pow(sigmay1331, 3);
            break;
        case 1332:
            dJydsigmay[1332] = 1.0/sigmay1332 - 1.0*pow(-my1332 + y1332, 2)/pow(sigmay1332, 3);
            break;
        case 1333:
            dJydsigmay[1333] = 1.0/sigmay1333 - 1.0*pow(-my1333 + y1333, 2)/pow(sigmay1333, 3);
            break;
        case 1334:
            dJydsigmay[1334] = 1.0/sigmay1334 - 1.0*pow(-my1334 + y1334, 2)/pow(sigmay1334, 3);
            break;
        case 1335:
            dJydsigmay[1335] = 1.0/sigmay1335 - 1.0*pow(-my1335 + y1335, 2)/pow(sigmay1335, 3);
            break;
        case 1336:
            dJydsigmay[1336] = 1.0/sigmay1336 - 1.0*pow(-my1336 + y1336, 2)/pow(sigmay1336, 3);
            break;
        case 1337:
            dJydsigmay[1337] = 1.0/sigmay1337 - 1.0*pow(-my1337 + y1337, 2)/pow(sigmay1337, 3);
            break;
        case 1338:
            dJydsigmay[1338] = 1.0/sigmay1338 - 1.0*pow(-my1338 + y1338, 2)/pow(sigmay1338, 3);
            break;
        case 1339:
            dJydsigmay[1339] = 1.0/sigmay1339 - 1.0*pow(-my1339 + y1339, 2)/pow(sigmay1339, 3);
            break;
        case 1340:
            dJydsigmay[1340] = 1.0/sigmay1340 - 1.0*pow(-my1340 + y1340, 2)/pow(sigmay1340, 3);
            break;
        case 1341:
            dJydsigmay[1341] = 1.0/sigmay1341 - 1.0*pow(-my1341 + y1341, 2)/pow(sigmay1341, 3);
            break;
        case 1342:
            dJydsigmay[1342] = 1.0/sigmay1342 - 1.0*pow(-my1342 + y1342, 2)/pow(sigmay1342, 3);
            break;
        case 1343:
            dJydsigmay[1343] = 1.0/sigmay1343 - 1.0*pow(-my1343 + y1343, 2)/pow(sigmay1343, 3);
            break;
        case 1344:
            dJydsigmay[1344] = 1.0/sigmay1344 - 1.0*pow(-my1344 + y1344, 2)/pow(sigmay1344, 3);
            break;
        case 1345:
            dJydsigmay[1345] = 1.0/sigmay1345 - 1.0*pow(-my1345 + y1345, 2)/pow(sigmay1345, 3);
            break;
        case 1346:
            dJydsigmay[1346] = 1.0/sigmay1346 - 1.0*pow(-my1346 + y1346, 2)/pow(sigmay1346, 3);
            break;
        case 1347:
            dJydsigmay[1347] = 1.0/sigmay1347 - 1.0*pow(-my1347 + y1347, 2)/pow(sigmay1347, 3);
            break;
        case 1348:
            dJydsigmay[1348] = 1.0/sigmay1348 - 1.0*pow(-my1348 + y1348, 2)/pow(sigmay1348, 3);
            break;
        case 1349:
            dJydsigmay[1349] = 1.0/sigmay1349 - 1.0*pow(-my1349 + y1349, 2)/pow(sigmay1349, 3);
            break;
        case 1350:
            dJydsigmay[1350] = 1.0/sigmay1350 - 1.0*pow(-my1350 + y1350, 2)/pow(sigmay1350, 3);
            break;
        case 1351:
            dJydsigmay[1351] = 1.0/sigmay1351 - 1.0*pow(-my1351 + y1351, 2)/pow(sigmay1351, 3);
            break;
        case 1352:
            dJydsigmay[1352] = 1.0/sigmay1352 - 1.0*pow(-my1352 + y1352, 2)/pow(sigmay1352, 3);
            break;
        case 1353:
            dJydsigmay[1353] = 1.0/sigmay1353 - 1.0*pow(-my1353 + y1353, 2)/pow(sigmay1353, 3);
            break;
        case 1354:
            dJydsigmay[1354] = 1.0/sigmay1354 - 1.0*pow(-my1354 + y1354, 2)/pow(sigmay1354, 3);
            break;
        case 1355:
            dJydsigmay[1355] = 1.0/sigmay1355 - 1.0*pow(-my1355 + y1355, 2)/pow(sigmay1355, 3);
            break;
        case 1356:
            dJydsigmay[1356] = 1.0/sigmay1356 - 1.0*pow(-my1356 + y1356, 2)/pow(sigmay1356, 3);
            break;
        case 1357:
            dJydsigmay[1357] = 1.0/sigmay1357 - 1.0*pow(-my1357 + y1357, 2)/pow(sigmay1357, 3);
            break;
        case 1358:
            dJydsigmay[1358] = 1.0/sigmay1358 - 1.0*pow(-my1358 + y1358, 2)/pow(sigmay1358, 3);
            break;
        case 1359:
            dJydsigmay[1359] = 1.0/sigmay1359 - 1.0*pow(-my1359 + y1359, 2)/pow(sigmay1359, 3);
            break;
        case 1360:
            dJydsigmay[1360] = 1.0/sigmay1360 - 1.0*pow(-my1360 + y1360, 2)/pow(sigmay1360, 3);
            break;
        case 1361:
            dJydsigmay[1361] = 1.0/sigmay1361 - 1.0*pow(-my1361 + y1361, 2)/pow(sigmay1361, 3);
            break;
        case 1362:
            dJydsigmay[1362] = 1.0/sigmay1362 - 1.0*pow(-my1362 + y1362, 2)/pow(sigmay1362, 3);
            break;
        case 1363:
            dJydsigmay[1363] = 1.0/sigmay1363 - 1.0*pow(-my1363 + y1363, 2)/pow(sigmay1363, 3);
            break;
        case 1364:
            dJydsigmay[1364] = 1.0/sigmay1364 - 1.0*pow(-my1364 + y1364, 2)/pow(sigmay1364, 3);
            break;
        case 1365:
            dJydsigmay[1365] = 1.0/sigmay1365 - 1.0*pow(-my1365 + y1365, 2)/pow(sigmay1365, 3);
            break;
        case 1366:
            dJydsigmay[1366] = 1.0/sigmay1366 - 1.0*pow(-my1366 + y1366, 2)/pow(sigmay1366, 3);
            break;
        case 1367:
            dJydsigmay[1367] = 1.0/sigmay1367 - 1.0*pow(-my1367 + y1367, 2)/pow(sigmay1367, 3);
            break;
        case 1368:
            dJydsigmay[1368] = 1.0/sigmay1368 - 1.0*pow(-my1368 + y1368, 2)/pow(sigmay1368, 3);
            break;
        case 1369:
            dJydsigmay[1369] = 1.0/sigmay1369 - 1.0*pow(-my1369 + y1369, 2)/pow(sigmay1369, 3);
            break;
        case 1370:
            dJydsigmay[1370] = 1.0/sigmay1370 - 1.0*pow(-my1370 + y1370, 2)/pow(sigmay1370, 3);
            break;
        case 1371:
            dJydsigmay[1371] = 1.0/sigmay1371 - 1.0*pow(-my1371 + y1371, 2)/pow(sigmay1371, 3);
            break;
        case 1372:
            dJydsigmay[1372] = 1.0/sigmay1372 - 1.0*pow(-my1372 + y1372, 2)/pow(sigmay1372, 3);
            break;
        case 1373:
            dJydsigmay[1373] = 1.0/sigmay1373 - 1.0*pow(-my1373 + y1373, 2)/pow(sigmay1373, 3);
            break;
        case 1374:
            dJydsigmay[1374] = 1.0/sigmay1374 - 1.0*pow(-my1374 + y1374, 2)/pow(sigmay1374, 3);
            break;
        case 1375:
            dJydsigmay[1375] = 1.0/sigmay1375 - 1.0*pow(-my1375 + y1375, 2)/pow(sigmay1375, 3);
            break;
        case 1376:
            dJydsigmay[1376] = 1.0/sigmay1376 - 1.0*pow(-my1376 + y1376, 2)/pow(sigmay1376, 3);
            break;
        case 1377:
            dJydsigmay[1377] = 1.0/sigmay1377 - 1.0*pow(-my1377 + y1377, 2)/pow(sigmay1377, 3);
            break;
        case 1378:
            dJydsigmay[1378] = 1.0/sigmay1378 - 1.0*pow(-my1378 + y1378, 2)/pow(sigmay1378, 3);
            break;
        case 1379:
            dJydsigmay[1379] = 1.0/sigmay1379 - 1.0*pow(-my1379 + y1379, 2)/pow(sigmay1379, 3);
            break;
        case 1380:
            dJydsigmay[1380] = 1.0/sigmay1380 - 1.0*pow(-my1380 + y1380, 2)/pow(sigmay1380, 3);
            break;
        case 1381:
            dJydsigmay[1381] = 1.0/sigmay1381 - 1.0*pow(-my1381 + y1381, 2)/pow(sigmay1381, 3);
            break;
        case 1382:
            dJydsigmay[1382] = 1.0/sigmay1382 - 1.0*pow(-my1382 + y1382, 2)/pow(sigmay1382, 3);
            break;
        case 1383:
            dJydsigmay[1383] = 1.0/sigmay1383 - 1.0*pow(-my1383 + y1383, 2)/pow(sigmay1383, 3);
            break;
        case 1384:
            dJydsigmay[1384] = 1.0/sigmay1384 - 1.0*pow(-my1384 + y1384, 2)/pow(sigmay1384, 3);
            break;
        case 1385:
            dJydsigmay[1385] = 1.0/sigmay1385 - 1.0*pow(-my1385 + y1385, 2)/pow(sigmay1385, 3);
            break;
        case 1386:
            dJydsigmay[1386] = 1.0/sigmay1386 - 1.0*pow(-my1386 + y1386, 2)/pow(sigmay1386, 3);
            break;
        case 1387:
            dJydsigmay[1387] = 1.0/sigmay1387 - 1.0*pow(-my1387 + y1387, 2)/pow(sigmay1387, 3);
            break;
        case 1388:
            dJydsigmay[1388] = 1.0/sigmay1388 - 1.0*pow(-my1388 + y1388, 2)/pow(sigmay1388, 3);
            break;
        case 1389:
            dJydsigmay[1389] = 1.0/sigmay1389 - 1.0*pow(-my1389 + y1389, 2)/pow(sigmay1389, 3);
            break;
        case 1390:
            dJydsigmay[1390] = 1.0/sigmay1390 - 1.0*pow(-my1390 + y1390, 2)/pow(sigmay1390, 3);
            break;
        case 1391:
            dJydsigmay[1391] = 1.0/sigmay1391 - 1.0*pow(-my1391 + y1391, 2)/pow(sigmay1391, 3);
            break;
        case 1392:
            dJydsigmay[1392] = 1.0/sigmay1392 - 1.0*pow(-my1392 + y1392, 2)/pow(sigmay1392, 3);
            break;
        case 1393:
            dJydsigmay[1393] = 1.0/sigmay1393 - 1.0*pow(-my1393 + y1393, 2)/pow(sigmay1393, 3);
            break;
        case 1394:
            dJydsigmay[1394] = 1.0/sigmay1394 - 1.0*pow(-my1394 + y1394, 2)/pow(sigmay1394, 3);
            break;
        case 1395:
            dJydsigmay[1395] = 1.0/sigmay1395 - 1.0*pow(-my1395 + y1395, 2)/pow(sigmay1395, 3);
            break;
    }
}