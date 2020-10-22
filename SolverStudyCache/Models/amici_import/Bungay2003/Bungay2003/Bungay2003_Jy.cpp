#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void Jy_Bungay2003(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay0, 2)) + 0.5*pow(-my0 + y0, 2)/pow(sigmay0, 2);
            break;
        case 1:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1, 2)) + 0.5*pow(-my1 + y1, 2)/pow(sigmay1, 2);
            break;
        case 2:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay2, 2)) + 0.5*pow(-my2 + y2, 2)/pow(sigmay2, 2);
            break;
        case 3:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay3, 2)) + 0.5*pow(-my3 + y3, 2)/pow(sigmay3, 2);
            break;
        case 4:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay4, 2)) + 0.5*pow(-my4 + y4, 2)/pow(sigmay4, 2);
            break;
        case 5:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay5, 2)) + 0.5*pow(-my5 + y5, 2)/pow(sigmay5, 2);
            break;
        case 6:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay6, 2)) + 0.5*pow(-my6 + y6, 2)/pow(sigmay6, 2);
            break;
        case 7:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay7, 2)) + 0.5*pow(-my7 + y7, 2)/pow(sigmay7, 2);
            break;
        case 8:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay8, 2)) + 0.5*pow(-my8 + y8, 2)/pow(sigmay8, 2);
            break;
        case 9:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay9, 2)) + 0.5*pow(-my9 + y9, 2)/pow(sigmay9, 2);
            break;
        case 10:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay10, 2)) + 0.5*pow(-my10 + y10, 2)/pow(sigmay10, 2);
            break;
        case 11:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay11, 2)) + 0.5*pow(-my11 + y11, 2)/pow(sigmay11, 2);
            break;
        case 12:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay12, 2)) + 0.5*pow(-my12 + y12, 2)/pow(sigmay12, 2);
            break;
        case 13:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay13, 2)) + 0.5*pow(-my13 + y13, 2)/pow(sigmay13, 2);
            break;
        case 14:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay14, 2)) + 0.5*pow(-my14 + y14, 2)/pow(sigmay14, 2);
            break;
        case 15:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay15, 2)) + 0.5*pow(-my15 + y15, 2)/pow(sigmay15, 2);
            break;
        case 16:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay16, 2)) + 0.5*pow(-my16 + y16, 2)/pow(sigmay16, 2);
            break;
        case 17:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay17, 2)) + 0.5*pow(-my17 + y17, 2)/pow(sigmay17, 2);
            break;
        case 18:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay18, 2)) + 0.5*pow(-my18 + y18, 2)/pow(sigmay18, 2);
            break;
        case 19:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay19, 2)) + 0.5*pow(-my19 + y19, 2)/pow(sigmay19, 2);
            break;
        case 20:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay20, 2)) + 0.5*pow(-my20 + y20, 2)/pow(sigmay20, 2);
            break;
        case 21:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay21, 2)) + 0.5*pow(-my21 + y21, 2)/pow(sigmay21, 2);
            break;
        case 22:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay22, 2)) + 0.5*pow(-my22 + y22, 2)/pow(sigmay22, 2);
            break;
        case 23:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay23, 2)) + 0.5*pow(-my23 + y23, 2)/pow(sigmay23, 2);
            break;
        case 24:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay24, 2)) + 0.5*pow(-my24 + y24, 2)/pow(sigmay24, 2);
            break;
        case 25:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay25, 2)) + 0.5*pow(-my25 + y25, 2)/pow(sigmay25, 2);
            break;
        case 26:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay26, 2)) + 0.5*pow(-my26 + y26, 2)/pow(sigmay26, 2);
            break;
        case 27:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay27, 2)) + 0.5*pow(-my27 + y27, 2)/pow(sigmay27, 2);
            break;
        case 28:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay28, 2)) + 0.5*pow(-my28 + y28, 2)/pow(sigmay28, 2);
            break;
        case 29:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay29, 2)) + 0.5*pow(-my29 + y29, 2)/pow(sigmay29, 2);
            break;
        case 30:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay30, 2)) + 0.5*pow(-my30 + y30, 2)/pow(sigmay30, 2);
            break;
        case 31:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay31, 2)) + 0.5*pow(-my31 + y31, 2)/pow(sigmay31, 2);
            break;
        case 32:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay32, 2)) + 0.5*pow(-my32 + y32, 2)/pow(sigmay32, 2);
            break;
        case 33:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay33, 2)) + 0.5*pow(-my33 + y33, 2)/pow(sigmay33, 2);
            break;
        case 34:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay34, 2)) + 0.5*pow(-my34 + y34, 2)/pow(sigmay34, 2);
            break;
        case 35:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay35, 2)) + 0.5*pow(-my35 + y35, 2)/pow(sigmay35, 2);
            break;
        case 36:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay36, 2)) + 0.5*pow(-my36 + y36, 2)/pow(sigmay36, 2);
            break;
        case 37:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay37, 2)) + 0.5*pow(-my37 + y37, 2)/pow(sigmay37, 2);
            break;
        case 38:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay38, 2)) + 0.5*pow(-my38 + y38, 2)/pow(sigmay38, 2);
            break;
        case 39:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay39, 2)) + 0.5*pow(-my39 + y39, 2)/pow(sigmay39, 2);
            break;
        case 40:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay40, 2)) + 0.5*pow(-my40 + y40, 2)/pow(sigmay40, 2);
            break;
        case 41:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay41, 2)) + 0.5*pow(-my41 + y41, 2)/pow(sigmay41, 2);
            break;
        case 42:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay42, 2)) + 0.5*pow(-my42 + y42, 2)/pow(sigmay42, 2);
            break;
        case 43:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay43, 2)) + 0.5*pow(-my43 + y43, 2)/pow(sigmay43, 2);
            break;
        case 44:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay44, 2)) + 0.5*pow(-my44 + y44, 2)/pow(sigmay44, 2);
            break;
        case 45:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay45, 2)) + 0.5*pow(-my45 + y45, 2)/pow(sigmay45, 2);
            break;
        case 46:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay46, 2)) + 0.5*pow(-my46 + y46, 2)/pow(sigmay46, 2);
            break;
        case 47:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay47, 2)) + 0.5*pow(-my47 + y47, 2)/pow(sigmay47, 2);
            break;
        case 48:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay48, 2)) + 0.5*pow(-my48 + y48, 2)/pow(sigmay48, 2);
            break;
        case 49:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay49, 2)) + 0.5*pow(-my49 + y49, 2)/pow(sigmay49, 2);
            break;
        case 50:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay50, 2)) + 0.5*pow(-my50 + y50, 2)/pow(sigmay50, 2);
            break;
        case 51:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay51, 2)) + 0.5*pow(-my51 + y51, 2)/pow(sigmay51, 2);
            break;
        case 52:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay52, 2)) + 0.5*pow(-my52 + y52, 2)/pow(sigmay52, 2);
            break;
        case 53:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay53, 2)) + 0.5*pow(-my53 + y53, 2)/pow(sigmay53, 2);
            break;
        case 54:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay54, 2)) + 0.5*pow(-my54 + y54, 2)/pow(sigmay54, 2);
            break;
        case 55:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay55, 2)) + 0.5*pow(-my55 + y55, 2)/pow(sigmay55, 2);
            break;
        case 56:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay56, 2)) + 0.5*pow(-my56 + y56, 2)/pow(sigmay56, 2);
            break;
        case 57:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay57, 2)) + 0.5*pow(-my57 + y57, 2)/pow(sigmay57, 2);
            break;
        case 58:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay58, 2)) + 0.5*pow(-my58 + y58, 2)/pow(sigmay58, 2);
            break;
        case 59:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay59, 2)) + 0.5*pow(-my59 + y59, 2)/pow(sigmay59, 2);
            break;
        case 60:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay60, 2)) + 0.5*pow(-my60 + y60, 2)/pow(sigmay60, 2);
            break;
        case 61:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay61, 2)) + 0.5*pow(-my61 + y61, 2)/pow(sigmay61, 2);
            break;
        case 62:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay62, 2)) + 0.5*pow(-my62 + y62, 2)/pow(sigmay62, 2);
            break;
        case 63:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay63, 2)) + 0.5*pow(-my63 + y63, 2)/pow(sigmay63, 2);
            break;
        case 64:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay64, 2)) + 0.5*pow(-my64 + y64, 2)/pow(sigmay64, 2);
            break;
        case 65:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay65, 2)) + 0.5*pow(-my65 + y65, 2)/pow(sigmay65, 2);
            break;
        case 66:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay66, 2)) + 0.5*pow(-my66 + y66, 2)/pow(sigmay66, 2);
            break;
        case 67:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay67, 2)) + 0.5*pow(-my67 + y67, 2)/pow(sigmay67, 2);
            break;
        case 68:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay68, 2)) + 0.5*pow(-my68 + y68, 2)/pow(sigmay68, 2);
            break;
        case 69:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay69, 2)) + 0.5*pow(-my69 + y69, 2)/pow(sigmay69, 2);
            break;
        case 70:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay70, 2)) + 0.5*pow(-my70 + y70, 2)/pow(sigmay70, 2);
            break;
        case 71:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay71, 2)) + 0.5*pow(-my71 + y71, 2)/pow(sigmay71, 2);
            break;
        case 72:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay72, 2)) + 0.5*pow(-my72 + y72, 2)/pow(sigmay72, 2);
            break;
        case 73:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay73, 2)) + 0.5*pow(-my73 + y73, 2)/pow(sigmay73, 2);
            break;
    }
}