#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void Jy_Ung2008(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
        case 74:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay74, 2)) + 0.5*pow(-my74 + y74, 2)/pow(sigmay74, 2);
            break;
        case 75:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay75, 2)) + 0.5*pow(-my75 + y75, 2)/pow(sigmay75, 2);
            break;
        case 76:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay76, 2)) + 0.5*pow(-my76 + y76, 2)/pow(sigmay76, 2);
            break;
        case 77:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay77, 2)) + 0.5*pow(-my77 + y77, 2)/pow(sigmay77, 2);
            break;
        case 78:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay78, 2)) + 0.5*pow(-my78 + y78, 2)/pow(sigmay78, 2);
            break;
        case 79:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay79, 2)) + 0.5*pow(-my79 + y79, 2)/pow(sigmay79, 2);
            break;
        case 80:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay80, 2)) + 0.5*pow(-my80 + y80, 2)/pow(sigmay80, 2);
            break;
        case 81:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay81, 2)) + 0.5*pow(-my81 + y81, 2)/pow(sigmay81, 2);
            break;
        case 82:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay82, 2)) + 0.5*pow(-my82 + y82, 2)/pow(sigmay82, 2);
            break;
        case 83:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay83, 2)) + 0.5*pow(-my83 + y83, 2)/pow(sigmay83, 2);
            break;
        case 84:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay84, 2)) + 0.5*pow(-my84 + y84, 2)/pow(sigmay84, 2);
            break;
        case 85:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay85, 2)) + 0.5*pow(-my85 + y85, 2)/pow(sigmay85, 2);
            break;
        case 86:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay86, 2)) + 0.5*pow(-my86 + y86, 2)/pow(sigmay86, 2);
            break;
        case 87:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay87, 2)) + 0.5*pow(-my87 + y87, 2)/pow(sigmay87, 2);
            break;
        case 88:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay88, 2)) + 0.5*pow(-my88 + y88, 2)/pow(sigmay88, 2);
            break;
        case 89:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay89, 2)) + 0.5*pow(-my89 + y89, 2)/pow(sigmay89, 2);
            break;
        case 90:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay90, 2)) + 0.5*pow(-my90 + y90, 2)/pow(sigmay90, 2);
            break;
        case 91:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay91, 2)) + 0.5*pow(-my91 + y91, 2)/pow(sigmay91, 2);
            break;
        case 92:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay92, 2)) + 0.5*pow(-my92 + y92, 2)/pow(sigmay92, 2);
            break;
        case 93:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay93, 2)) + 0.5*pow(-my93 + y93, 2)/pow(sigmay93, 2);
            break;
        case 94:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay94, 2)) + 0.5*pow(-my94 + y94, 2)/pow(sigmay94, 2);
            break;
        case 95:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay95, 2)) + 0.5*pow(-my95 + y95, 2)/pow(sigmay95, 2);
            break;
        case 96:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay96, 2)) + 0.5*pow(-my96 + y96, 2)/pow(sigmay96, 2);
            break;
        case 97:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay97, 2)) + 0.5*pow(-my97 + y97, 2)/pow(sigmay97, 2);
            break;
        case 98:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay98, 2)) + 0.5*pow(-my98 + y98, 2)/pow(sigmay98, 2);
            break;
        case 99:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay99, 2)) + 0.5*pow(-my99 + y99, 2)/pow(sigmay99, 2);
            break;
        case 100:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay100, 2)) + 0.5*pow(-my100 + y100, 2)/pow(sigmay100, 2);
            break;
        case 101:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay101, 2)) + 0.5*pow(-my101 + y101, 2)/pow(sigmay101, 2);
            break;
        case 102:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay102, 2)) + 0.5*pow(-my102 + y102, 2)/pow(sigmay102, 2);
            break;
        case 103:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay103, 2)) + 0.5*pow(-my103 + y103, 2)/pow(sigmay103, 2);
            break;
        case 104:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay104, 2)) + 0.5*pow(-my104 + y104, 2)/pow(sigmay104, 2);
            break;
        case 105:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay105, 2)) + 0.5*pow(-my105 + y105, 2)/pow(sigmay105, 2);
            break;
        case 106:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay106, 2)) + 0.5*pow(-my106 + y106, 2)/pow(sigmay106, 2);
            break;
        case 107:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay107, 2)) + 0.5*pow(-my107 + y107, 2)/pow(sigmay107, 2);
            break;
        case 108:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay108, 2)) + 0.5*pow(-my108 + y108, 2)/pow(sigmay108, 2);
            break;
        case 109:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay109, 2)) + 0.5*pow(-my109 + y109, 2)/pow(sigmay109, 2);
            break;
        case 110:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay110, 2)) + 0.5*pow(-my110 + y110, 2)/pow(sigmay110, 2);
            break;
        case 111:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay111, 2)) + 0.5*pow(-my111 + y111, 2)/pow(sigmay111, 2);
            break;
        case 112:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay112, 2)) + 0.5*pow(-my112 + y112, 2)/pow(sigmay112, 2);
            break;
        case 113:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay113, 2)) + 0.5*pow(-my113 + y113, 2)/pow(sigmay113, 2);
            break;
        case 114:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay114, 2)) + 0.5*pow(-my114 + y114, 2)/pow(sigmay114, 2);
            break;
        case 115:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay115, 2)) + 0.5*pow(-my115 + y115, 2)/pow(sigmay115, 2);
            break;
        case 116:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay116, 2)) + 0.5*pow(-my116 + y116, 2)/pow(sigmay116, 2);
            break;
        case 117:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay117, 2)) + 0.5*pow(-my117 + y117, 2)/pow(sigmay117, 2);
            break;
        case 118:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay118, 2)) + 0.5*pow(-my118 + y118, 2)/pow(sigmay118, 2);
            break;
        case 119:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay119, 2)) + 0.5*pow(-my119 + y119, 2)/pow(sigmay119, 2);
            break;
        case 120:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay120, 2)) + 0.5*pow(-my120 + y120, 2)/pow(sigmay120, 2);
            break;
        case 121:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay121, 2)) + 0.5*pow(-my121 + y121, 2)/pow(sigmay121, 2);
            break;
        case 122:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay122, 2)) + 0.5*pow(-my122 + y122, 2)/pow(sigmay122, 2);
            break;
        case 123:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay123, 2)) + 0.5*pow(-my123 + y123, 2)/pow(sigmay123, 2);
            break;
        case 124:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay124, 2)) + 0.5*pow(-my124 + y124, 2)/pow(sigmay124, 2);
            break;
        case 125:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay125, 2)) + 0.5*pow(-my125 + y125, 2)/pow(sigmay125, 2);
            break;
        case 126:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay126, 2)) + 0.5*pow(-my126 + y126, 2)/pow(sigmay126, 2);
            break;
        case 127:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay127, 2)) + 0.5*pow(-my127 + y127, 2)/pow(sigmay127, 2);
            break;
        case 128:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay128, 2)) + 0.5*pow(-my128 + y128, 2)/pow(sigmay128, 2);
            break;
        case 129:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay129, 2)) + 0.5*pow(-my129 + y129, 2)/pow(sigmay129, 2);
            break;
        case 130:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay130, 2)) + 0.5*pow(-my130 + y130, 2)/pow(sigmay130, 2);
            break;
        case 131:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay131, 2)) + 0.5*pow(-my131 + y131, 2)/pow(sigmay131, 2);
            break;
        case 132:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay132, 2)) + 0.5*pow(-my132 + y132, 2)/pow(sigmay132, 2);
            break;
        case 133:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay133, 2)) + 0.5*pow(-my133 + y133, 2)/pow(sigmay133, 2);
            break;
        case 134:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay134, 2)) + 0.5*pow(-my134 + y134, 2)/pow(sigmay134, 2);
            break;
        case 135:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay135, 2)) + 0.5*pow(-my135 + y135, 2)/pow(sigmay135, 2);
            break;
        case 136:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay136, 2)) + 0.5*pow(-my136 + y136, 2)/pow(sigmay136, 2);
            break;
        case 137:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay137, 2)) + 0.5*pow(-my137 + y137, 2)/pow(sigmay137, 2);
            break;
        case 138:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay138, 2)) + 0.5*pow(-my138 + y138, 2)/pow(sigmay138, 2);
            break;
        case 139:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay139, 2)) + 0.5*pow(-my139 + y139, 2)/pow(sigmay139, 2);
            break;
        case 140:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay140, 2)) + 0.5*pow(-my140 + y140, 2)/pow(sigmay140, 2);
            break;
        case 141:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay141, 2)) + 0.5*pow(-my141 + y141, 2)/pow(sigmay141, 2);
            break;
        case 142:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay142, 2)) + 0.5*pow(-my142 + y142, 2)/pow(sigmay142, 2);
            break;
        case 143:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay143, 2)) + 0.5*pow(-my143 + y143, 2)/pow(sigmay143, 2);
            break;
        case 144:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay144, 2)) + 0.5*pow(-my144 + y144, 2)/pow(sigmay144, 2);
            break;
        case 145:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay145, 2)) + 0.5*pow(-my145 + y145, 2)/pow(sigmay145, 2);
            break;
        case 146:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay146, 2)) + 0.5*pow(-my146 + y146, 2)/pow(sigmay146, 2);
            break;
        case 147:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay147, 2)) + 0.5*pow(-my147 + y147, 2)/pow(sigmay147, 2);
            break;
        case 148:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay148, 2)) + 0.5*pow(-my148 + y148, 2)/pow(sigmay148, 2);
            break;
        case 149:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay149, 2)) + 0.5*pow(-my149 + y149, 2)/pow(sigmay149, 2);
            break;
        case 150:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay150, 2)) + 0.5*pow(-my150 + y150, 2)/pow(sigmay150, 2);
            break;
        case 151:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay151, 2)) + 0.5*pow(-my151 + y151, 2)/pow(sigmay151, 2);
            break;
        case 152:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay152, 2)) + 0.5*pow(-my152 + y152, 2)/pow(sigmay152, 2);
            break;
        case 153:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay153, 2)) + 0.5*pow(-my153 + y153, 2)/pow(sigmay153, 2);
            break;
        case 154:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay154, 2)) + 0.5*pow(-my154 + y154, 2)/pow(sigmay154, 2);
            break;
        case 155:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay155, 2)) + 0.5*pow(-my155 + y155, 2)/pow(sigmay155, 2);
            break;
        case 156:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay156, 2)) + 0.5*pow(-my156 + y156, 2)/pow(sigmay156, 2);
            break;
        case 157:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay157, 2)) + 0.5*pow(-my157 + y157, 2)/pow(sigmay157, 2);
            break;
        case 158:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay158, 2)) + 0.5*pow(-my158 + y158, 2)/pow(sigmay158, 2);
            break;
        case 159:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay159, 2)) + 0.5*pow(-my159 + y159, 2)/pow(sigmay159, 2);
            break;
        case 160:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay160, 2)) + 0.5*pow(-my160 + y160, 2)/pow(sigmay160, 2);
            break;
        case 161:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay161, 2)) + 0.5*pow(-my161 + y161, 2)/pow(sigmay161, 2);
            break;
        case 162:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay162, 2)) + 0.5*pow(-my162 + y162, 2)/pow(sigmay162, 2);
            break;
        case 163:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay163, 2)) + 0.5*pow(-my163 + y163, 2)/pow(sigmay163, 2);
            break;
        case 164:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay164, 2)) + 0.5*pow(-my164 + y164, 2)/pow(sigmay164, 2);
            break;
        case 165:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay165, 2)) + 0.5*pow(-my165 + y165, 2)/pow(sigmay165, 2);
            break;
        case 166:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay166, 2)) + 0.5*pow(-my166 + y166, 2)/pow(sigmay166, 2);
            break;
        case 167:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay167, 2)) + 0.5*pow(-my167 + y167, 2)/pow(sigmay167, 2);
            break;
        case 168:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay168, 2)) + 0.5*pow(-my168 + y168, 2)/pow(sigmay168, 2);
            break;
        case 169:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay169, 2)) + 0.5*pow(-my169 + y169, 2)/pow(sigmay169, 2);
            break;
        case 170:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay170, 2)) + 0.5*pow(-my170 + y170, 2)/pow(sigmay170, 2);
            break;
        case 171:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay171, 2)) + 0.5*pow(-my171 + y171, 2)/pow(sigmay171, 2);
            break;
        case 172:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay172, 2)) + 0.5*pow(-my172 + y172, 2)/pow(sigmay172, 2);
            break;
        case 173:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay173, 2)) + 0.5*pow(-my173 + y173, 2)/pow(sigmay173, 2);
            break;
        case 174:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay174, 2)) + 0.5*pow(-my174 + y174, 2)/pow(sigmay174, 2);
            break;
        case 175:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay175, 2)) + 0.5*pow(-my175 + y175, 2)/pow(sigmay175, 2);
            break;
        case 176:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay176, 2)) + 0.5*pow(-my176 + y176, 2)/pow(sigmay176, 2);
            break;
        case 177:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay177, 2)) + 0.5*pow(-my177 + y177, 2)/pow(sigmay177, 2);
            break;
        case 178:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay178, 2)) + 0.5*pow(-my178 + y178, 2)/pow(sigmay178, 2);
            break;
        case 179:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay179, 2)) + 0.5*pow(-my179 + y179, 2)/pow(sigmay179, 2);
            break;
        case 180:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay180, 2)) + 0.5*pow(-my180 + y180, 2)/pow(sigmay180, 2);
            break;
        case 181:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay181, 2)) + 0.5*pow(-my181 + y181, 2)/pow(sigmay181, 2);
            break;
        case 182:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay182, 2)) + 0.5*pow(-my182 + y182, 2)/pow(sigmay182, 2);
            break;
        case 183:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay183, 2)) + 0.5*pow(-my183 + y183, 2)/pow(sigmay183, 2);
            break;
        case 184:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay184, 2)) + 0.5*pow(-my184 + y184, 2)/pow(sigmay184, 2);
            break;
        case 185:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay185, 2)) + 0.5*pow(-my185 + y185, 2)/pow(sigmay185, 2);
            break;
        case 186:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay186, 2)) + 0.5*pow(-my186 + y186, 2)/pow(sigmay186, 2);
            break;
        case 187:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay187, 2)) + 0.5*pow(-my187 + y187, 2)/pow(sigmay187, 2);
            break;
        case 188:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay188, 2)) + 0.5*pow(-my188 + y188, 2)/pow(sigmay188, 2);
            break;
        case 189:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay189, 2)) + 0.5*pow(-my189 + y189, 2)/pow(sigmay189, 2);
            break;
        case 190:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay190, 2)) + 0.5*pow(-my190 + y190, 2)/pow(sigmay190, 2);
            break;
        case 191:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay191, 2)) + 0.5*pow(-my191 + y191, 2)/pow(sigmay191, 2);
            break;
        case 192:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay192, 2)) + 0.5*pow(-my192 + y192, 2)/pow(sigmay192, 2);
            break;
        case 193:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay193, 2)) + 0.5*pow(-my193 + y193, 2)/pow(sigmay193, 2);
            break;
    }
}