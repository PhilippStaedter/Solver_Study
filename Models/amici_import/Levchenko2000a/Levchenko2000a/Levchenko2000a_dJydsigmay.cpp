#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydsigmay_Levchenko2000a(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
    }
}