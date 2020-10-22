#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void Jy_Froehlich2018(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
        case 194:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay194, 2)) + 0.5*pow(-my194 + y194, 2)/pow(sigmay194, 2);
            break;
        case 195:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay195, 2)) + 0.5*pow(-my195 + y195, 2)/pow(sigmay195, 2);
            break;
        case 196:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay196, 2)) + 0.5*pow(-my196 + y196, 2)/pow(sigmay196, 2);
            break;
        case 197:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay197, 2)) + 0.5*pow(-my197 + y197, 2)/pow(sigmay197, 2);
            break;
        case 198:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay198, 2)) + 0.5*pow(-my198 + y198, 2)/pow(sigmay198, 2);
            break;
        case 199:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay199, 2)) + 0.5*pow(-my199 + y199, 2)/pow(sigmay199, 2);
            break;
        case 200:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay200, 2)) + 0.5*pow(-my200 + y200, 2)/pow(sigmay200, 2);
            break;
        case 201:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay201, 2)) + 0.5*pow(-my201 + y201, 2)/pow(sigmay201, 2);
            break;
        case 202:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay202, 2)) + 0.5*pow(-my202 + y202, 2)/pow(sigmay202, 2);
            break;
        case 203:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay203, 2)) + 0.5*pow(-my203 + y203, 2)/pow(sigmay203, 2);
            break;
        case 204:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay204, 2)) + 0.5*pow(-my204 + y204, 2)/pow(sigmay204, 2);
            break;
        case 205:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay205, 2)) + 0.5*pow(-my205 + y205, 2)/pow(sigmay205, 2);
            break;
        case 206:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay206, 2)) + 0.5*pow(-my206 + y206, 2)/pow(sigmay206, 2);
            break;
        case 207:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay207, 2)) + 0.5*pow(-my207 + y207, 2)/pow(sigmay207, 2);
            break;
        case 208:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay208, 2)) + 0.5*pow(-my208 + y208, 2)/pow(sigmay208, 2);
            break;
        case 209:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay209, 2)) + 0.5*pow(-my209 + y209, 2)/pow(sigmay209, 2);
            break;
        case 210:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay210, 2)) + 0.5*pow(-my210 + y210, 2)/pow(sigmay210, 2);
            break;
        case 211:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay211, 2)) + 0.5*pow(-my211 + y211, 2)/pow(sigmay211, 2);
            break;
        case 212:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay212, 2)) + 0.5*pow(-my212 + y212, 2)/pow(sigmay212, 2);
            break;
        case 213:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay213, 2)) + 0.5*pow(-my213 + y213, 2)/pow(sigmay213, 2);
            break;
        case 214:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay214, 2)) + 0.5*pow(-my214 + y214, 2)/pow(sigmay214, 2);
            break;
        case 215:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay215, 2)) + 0.5*pow(-my215 + y215, 2)/pow(sigmay215, 2);
            break;
        case 216:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay216, 2)) + 0.5*pow(-my216 + y216, 2)/pow(sigmay216, 2);
            break;
        case 217:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay217, 2)) + 0.5*pow(-my217 + y217, 2)/pow(sigmay217, 2);
            break;
        case 218:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay218, 2)) + 0.5*pow(-my218 + y218, 2)/pow(sigmay218, 2);
            break;
        case 219:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay219, 2)) + 0.5*pow(-my219 + y219, 2)/pow(sigmay219, 2);
            break;
        case 220:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay220, 2)) + 0.5*pow(-my220 + y220, 2)/pow(sigmay220, 2);
            break;
        case 221:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay221, 2)) + 0.5*pow(-my221 + y221, 2)/pow(sigmay221, 2);
            break;
        case 222:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay222, 2)) + 0.5*pow(-my222 + y222, 2)/pow(sigmay222, 2);
            break;
        case 223:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay223, 2)) + 0.5*pow(-my223 + y223, 2)/pow(sigmay223, 2);
            break;
        case 224:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay224, 2)) + 0.5*pow(-my224 + y224, 2)/pow(sigmay224, 2);
            break;
        case 225:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay225, 2)) + 0.5*pow(-my225 + y225, 2)/pow(sigmay225, 2);
            break;
        case 226:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay226, 2)) + 0.5*pow(-my226 + y226, 2)/pow(sigmay226, 2);
            break;
        case 227:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay227, 2)) + 0.5*pow(-my227 + y227, 2)/pow(sigmay227, 2);
            break;
        case 228:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay228, 2)) + 0.5*pow(-my228 + y228, 2)/pow(sigmay228, 2);
            break;
        case 229:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay229, 2)) + 0.5*pow(-my229 + y229, 2)/pow(sigmay229, 2);
            break;
        case 230:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay230, 2)) + 0.5*pow(-my230 + y230, 2)/pow(sigmay230, 2);
            break;
        case 231:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay231, 2)) + 0.5*pow(-my231 + y231, 2)/pow(sigmay231, 2);
            break;
        case 232:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay232, 2)) + 0.5*pow(-my232 + y232, 2)/pow(sigmay232, 2);
            break;
        case 233:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay233, 2)) + 0.5*pow(-my233 + y233, 2)/pow(sigmay233, 2);
            break;
        case 234:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay234, 2)) + 0.5*pow(-my234 + y234, 2)/pow(sigmay234, 2);
            break;
        case 235:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay235, 2)) + 0.5*pow(-my235 + y235, 2)/pow(sigmay235, 2);
            break;
        case 236:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay236, 2)) + 0.5*pow(-my236 + y236, 2)/pow(sigmay236, 2);
            break;
        case 237:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay237, 2)) + 0.5*pow(-my237 + y237, 2)/pow(sigmay237, 2);
            break;
        case 238:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay238, 2)) + 0.5*pow(-my238 + y238, 2)/pow(sigmay238, 2);
            break;
        case 239:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay239, 2)) + 0.5*pow(-my239 + y239, 2)/pow(sigmay239, 2);
            break;
        case 240:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay240, 2)) + 0.5*pow(-my240 + y240, 2)/pow(sigmay240, 2);
            break;
        case 241:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay241, 2)) + 0.5*pow(-my241 + y241, 2)/pow(sigmay241, 2);
            break;
        case 242:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay242, 2)) + 0.5*pow(-my242 + y242, 2)/pow(sigmay242, 2);
            break;
        case 243:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay243, 2)) + 0.5*pow(-my243 + y243, 2)/pow(sigmay243, 2);
            break;
        case 244:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay244, 2)) + 0.5*pow(-my244 + y244, 2)/pow(sigmay244, 2);
            break;
        case 245:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay245, 2)) + 0.5*pow(-my245 + y245, 2)/pow(sigmay245, 2);
            break;
        case 246:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay246, 2)) + 0.5*pow(-my246 + y246, 2)/pow(sigmay246, 2);
            break;
        case 247:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay247, 2)) + 0.5*pow(-my247 + y247, 2)/pow(sigmay247, 2);
            break;
        case 248:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay248, 2)) + 0.5*pow(-my248 + y248, 2)/pow(sigmay248, 2);
            break;
        case 249:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay249, 2)) + 0.5*pow(-my249 + y249, 2)/pow(sigmay249, 2);
            break;
        case 250:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay250, 2)) + 0.5*pow(-my250 + y250, 2)/pow(sigmay250, 2);
            break;
        case 251:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay251, 2)) + 0.5*pow(-my251 + y251, 2)/pow(sigmay251, 2);
            break;
        case 252:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay252, 2)) + 0.5*pow(-my252 + y252, 2)/pow(sigmay252, 2);
            break;
        case 253:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay253, 2)) + 0.5*pow(-my253 + y253, 2)/pow(sigmay253, 2);
            break;
        case 254:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay254, 2)) + 0.5*pow(-my254 + y254, 2)/pow(sigmay254, 2);
            break;
        case 255:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay255, 2)) + 0.5*pow(-my255 + y255, 2)/pow(sigmay255, 2);
            break;
        case 256:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay256, 2)) + 0.5*pow(-my256 + y256, 2)/pow(sigmay256, 2);
            break;
        case 257:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay257, 2)) + 0.5*pow(-my257 + y257, 2)/pow(sigmay257, 2);
            break;
        case 258:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay258, 2)) + 0.5*pow(-my258 + y258, 2)/pow(sigmay258, 2);
            break;
        case 259:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay259, 2)) + 0.5*pow(-my259 + y259, 2)/pow(sigmay259, 2);
            break;
        case 260:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay260, 2)) + 0.5*pow(-my260 + y260, 2)/pow(sigmay260, 2);
            break;
        case 261:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay261, 2)) + 0.5*pow(-my261 + y261, 2)/pow(sigmay261, 2);
            break;
        case 262:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay262, 2)) + 0.5*pow(-my262 + y262, 2)/pow(sigmay262, 2);
            break;
        case 263:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay263, 2)) + 0.5*pow(-my263 + y263, 2)/pow(sigmay263, 2);
            break;
        case 264:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay264, 2)) + 0.5*pow(-my264 + y264, 2)/pow(sigmay264, 2);
            break;
        case 265:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay265, 2)) + 0.5*pow(-my265 + y265, 2)/pow(sigmay265, 2);
            break;
        case 266:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay266, 2)) + 0.5*pow(-my266 + y266, 2)/pow(sigmay266, 2);
            break;
        case 267:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay267, 2)) + 0.5*pow(-my267 + y267, 2)/pow(sigmay267, 2);
            break;
        case 268:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay268, 2)) + 0.5*pow(-my268 + y268, 2)/pow(sigmay268, 2);
            break;
        case 269:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay269, 2)) + 0.5*pow(-my269 + y269, 2)/pow(sigmay269, 2);
            break;
        case 270:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay270, 2)) + 0.5*pow(-my270 + y270, 2)/pow(sigmay270, 2);
            break;
        case 271:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay271, 2)) + 0.5*pow(-my271 + y271, 2)/pow(sigmay271, 2);
            break;
        case 272:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay272, 2)) + 0.5*pow(-my272 + y272, 2)/pow(sigmay272, 2);
            break;
        case 273:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay273, 2)) + 0.5*pow(-my273 + y273, 2)/pow(sigmay273, 2);
            break;
        case 274:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay274, 2)) + 0.5*pow(-my274 + y274, 2)/pow(sigmay274, 2);
            break;
        case 275:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay275, 2)) + 0.5*pow(-my275 + y275, 2)/pow(sigmay275, 2);
            break;
        case 276:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay276, 2)) + 0.5*pow(-my276 + y276, 2)/pow(sigmay276, 2);
            break;
        case 277:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay277, 2)) + 0.5*pow(-my277 + y277, 2)/pow(sigmay277, 2);
            break;
        case 278:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay278, 2)) + 0.5*pow(-my278 + y278, 2)/pow(sigmay278, 2);
            break;
        case 279:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay279, 2)) + 0.5*pow(-my279 + y279, 2)/pow(sigmay279, 2);
            break;
        case 280:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay280, 2)) + 0.5*pow(-my280 + y280, 2)/pow(sigmay280, 2);
            break;
        case 281:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay281, 2)) + 0.5*pow(-my281 + y281, 2)/pow(sigmay281, 2);
            break;
        case 282:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay282, 2)) + 0.5*pow(-my282 + y282, 2)/pow(sigmay282, 2);
            break;
        case 283:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay283, 2)) + 0.5*pow(-my283 + y283, 2)/pow(sigmay283, 2);
            break;
        case 284:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay284, 2)) + 0.5*pow(-my284 + y284, 2)/pow(sigmay284, 2);
            break;
        case 285:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay285, 2)) + 0.5*pow(-my285 + y285, 2)/pow(sigmay285, 2);
            break;
        case 286:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay286, 2)) + 0.5*pow(-my286 + y286, 2)/pow(sigmay286, 2);
            break;
        case 287:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay287, 2)) + 0.5*pow(-my287 + y287, 2)/pow(sigmay287, 2);
            break;
        case 288:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay288, 2)) + 0.5*pow(-my288 + y288, 2)/pow(sigmay288, 2);
            break;
        case 289:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay289, 2)) + 0.5*pow(-my289 + y289, 2)/pow(sigmay289, 2);
            break;
        case 290:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay290, 2)) + 0.5*pow(-my290 + y290, 2)/pow(sigmay290, 2);
            break;
        case 291:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay291, 2)) + 0.5*pow(-my291 + y291, 2)/pow(sigmay291, 2);
            break;
        case 292:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay292, 2)) + 0.5*pow(-my292 + y292, 2)/pow(sigmay292, 2);
            break;
        case 293:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay293, 2)) + 0.5*pow(-my293 + y293, 2)/pow(sigmay293, 2);
            break;
        case 294:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay294, 2)) + 0.5*pow(-my294 + y294, 2)/pow(sigmay294, 2);
            break;
        case 295:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay295, 2)) + 0.5*pow(-my295 + y295, 2)/pow(sigmay295, 2);
            break;
        case 296:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay296, 2)) + 0.5*pow(-my296 + y296, 2)/pow(sigmay296, 2);
            break;
        case 297:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay297, 2)) + 0.5*pow(-my297 + y297, 2)/pow(sigmay297, 2);
            break;
        case 298:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay298, 2)) + 0.5*pow(-my298 + y298, 2)/pow(sigmay298, 2);
            break;
        case 299:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay299, 2)) + 0.5*pow(-my299 + y299, 2)/pow(sigmay299, 2);
            break;
        case 300:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay300, 2)) + 0.5*pow(-my300 + y300, 2)/pow(sigmay300, 2);
            break;
        case 301:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay301, 2)) + 0.5*pow(-my301 + y301, 2)/pow(sigmay301, 2);
            break;
        case 302:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay302, 2)) + 0.5*pow(-my302 + y302, 2)/pow(sigmay302, 2);
            break;
        case 303:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay303, 2)) + 0.5*pow(-my303 + y303, 2)/pow(sigmay303, 2);
            break;
        case 304:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay304, 2)) + 0.5*pow(-my304 + y304, 2)/pow(sigmay304, 2);
            break;
        case 305:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay305, 2)) + 0.5*pow(-my305 + y305, 2)/pow(sigmay305, 2);
            break;
        case 306:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay306, 2)) + 0.5*pow(-my306 + y306, 2)/pow(sigmay306, 2);
            break;
        case 307:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay307, 2)) + 0.5*pow(-my307 + y307, 2)/pow(sigmay307, 2);
            break;
        case 308:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay308, 2)) + 0.5*pow(-my308 + y308, 2)/pow(sigmay308, 2);
            break;
        case 309:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay309, 2)) + 0.5*pow(-my309 + y309, 2)/pow(sigmay309, 2);
            break;
        case 310:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay310, 2)) + 0.5*pow(-my310 + y310, 2)/pow(sigmay310, 2);
            break;
        case 311:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay311, 2)) + 0.5*pow(-my311 + y311, 2)/pow(sigmay311, 2);
            break;
        case 312:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay312, 2)) + 0.5*pow(-my312 + y312, 2)/pow(sigmay312, 2);
            break;
        case 313:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay313, 2)) + 0.5*pow(-my313 + y313, 2)/pow(sigmay313, 2);
            break;
        case 314:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay314, 2)) + 0.5*pow(-my314 + y314, 2)/pow(sigmay314, 2);
            break;
        case 315:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay315, 2)) + 0.5*pow(-my315 + y315, 2)/pow(sigmay315, 2);
            break;
        case 316:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay316, 2)) + 0.5*pow(-my316 + y316, 2)/pow(sigmay316, 2);
            break;
        case 317:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay317, 2)) + 0.5*pow(-my317 + y317, 2)/pow(sigmay317, 2);
            break;
        case 318:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay318, 2)) + 0.5*pow(-my318 + y318, 2)/pow(sigmay318, 2);
            break;
        case 319:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay319, 2)) + 0.5*pow(-my319 + y319, 2)/pow(sigmay319, 2);
            break;
        case 320:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay320, 2)) + 0.5*pow(-my320 + y320, 2)/pow(sigmay320, 2);
            break;
        case 321:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay321, 2)) + 0.5*pow(-my321 + y321, 2)/pow(sigmay321, 2);
            break;
        case 322:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay322, 2)) + 0.5*pow(-my322 + y322, 2)/pow(sigmay322, 2);
            break;
        case 323:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay323, 2)) + 0.5*pow(-my323 + y323, 2)/pow(sigmay323, 2);
            break;
        case 324:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay324, 2)) + 0.5*pow(-my324 + y324, 2)/pow(sigmay324, 2);
            break;
        case 325:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay325, 2)) + 0.5*pow(-my325 + y325, 2)/pow(sigmay325, 2);
            break;
        case 326:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay326, 2)) + 0.5*pow(-my326 + y326, 2)/pow(sigmay326, 2);
            break;
        case 327:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay327, 2)) + 0.5*pow(-my327 + y327, 2)/pow(sigmay327, 2);
            break;
        case 328:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay328, 2)) + 0.5*pow(-my328 + y328, 2)/pow(sigmay328, 2);
            break;
        case 329:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay329, 2)) + 0.5*pow(-my329 + y329, 2)/pow(sigmay329, 2);
            break;
        case 330:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay330, 2)) + 0.5*pow(-my330 + y330, 2)/pow(sigmay330, 2);
            break;
        case 331:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay331, 2)) + 0.5*pow(-my331 + y331, 2)/pow(sigmay331, 2);
            break;
        case 332:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay332, 2)) + 0.5*pow(-my332 + y332, 2)/pow(sigmay332, 2);
            break;
        case 333:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay333, 2)) + 0.5*pow(-my333 + y333, 2)/pow(sigmay333, 2);
            break;
        case 334:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay334, 2)) + 0.5*pow(-my334 + y334, 2)/pow(sigmay334, 2);
            break;
        case 335:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay335, 2)) + 0.5*pow(-my335 + y335, 2)/pow(sigmay335, 2);
            break;
        case 336:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay336, 2)) + 0.5*pow(-my336 + y336, 2)/pow(sigmay336, 2);
            break;
        case 337:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay337, 2)) + 0.5*pow(-my337 + y337, 2)/pow(sigmay337, 2);
            break;
        case 338:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay338, 2)) + 0.5*pow(-my338 + y338, 2)/pow(sigmay338, 2);
            break;
        case 339:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay339, 2)) + 0.5*pow(-my339 + y339, 2)/pow(sigmay339, 2);
            break;
        case 340:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay340, 2)) + 0.5*pow(-my340 + y340, 2)/pow(sigmay340, 2);
            break;
        case 341:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay341, 2)) + 0.5*pow(-my341 + y341, 2)/pow(sigmay341, 2);
            break;
        case 342:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay342, 2)) + 0.5*pow(-my342 + y342, 2)/pow(sigmay342, 2);
            break;
        case 343:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay343, 2)) + 0.5*pow(-my343 + y343, 2)/pow(sigmay343, 2);
            break;
        case 344:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay344, 2)) + 0.5*pow(-my344 + y344, 2)/pow(sigmay344, 2);
            break;
        case 345:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay345, 2)) + 0.5*pow(-my345 + y345, 2)/pow(sigmay345, 2);
            break;
        case 346:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay346, 2)) + 0.5*pow(-my346 + y346, 2)/pow(sigmay346, 2);
            break;
        case 347:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay347, 2)) + 0.5*pow(-my347 + y347, 2)/pow(sigmay347, 2);
            break;
        case 348:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay348, 2)) + 0.5*pow(-my348 + y348, 2)/pow(sigmay348, 2);
            break;
        case 349:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay349, 2)) + 0.5*pow(-my349 + y349, 2)/pow(sigmay349, 2);
            break;
        case 350:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay350, 2)) + 0.5*pow(-my350 + y350, 2)/pow(sigmay350, 2);
            break;
        case 351:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay351, 2)) + 0.5*pow(-my351 + y351, 2)/pow(sigmay351, 2);
            break;
        case 352:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay352, 2)) + 0.5*pow(-my352 + y352, 2)/pow(sigmay352, 2);
            break;
        case 353:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay353, 2)) + 0.5*pow(-my353 + y353, 2)/pow(sigmay353, 2);
            break;
        case 354:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay354, 2)) + 0.5*pow(-my354 + y354, 2)/pow(sigmay354, 2);
            break;
        case 355:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay355, 2)) + 0.5*pow(-my355 + y355, 2)/pow(sigmay355, 2);
            break;
        case 356:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay356, 2)) + 0.5*pow(-my356 + y356, 2)/pow(sigmay356, 2);
            break;
        case 357:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay357, 2)) + 0.5*pow(-my357 + y357, 2)/pow(sigmay357, 2);
            break;
        case 358:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay358, 2)) + 0.5*pow(-my358 + y358, 2)/pow(sigmay358, 2);
            break;
        case 359:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay359, 2)) + 0.5*pow(-my359 + y359, 2)/pow(sigmay359, 2);
            break;
        case 360:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay360, 2)) + 0.5*pow(-my360 + y360, 2)/pow(sigmay360, 2);
            break;
        case 361:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay361, 2)) + 0.5*pow(-my361 + y361, 2)/pow(sigmay361, 2);
            break;
        case 362:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay362, 2)) + 0.5*pow(-my362 + y362, 2)/pow(sigmay362, 2);
            break;
        case 363:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay363, 2)) + 0.5*pow(-my363 + y363, 2)/pow(sigmay363, 2);
            break;
        case 364:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay364, 2)) + 0.5*pow(-my364 + y364, 2)/pow(sigmay364, 2);
            break;
        case 365:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay365, 2)) + 0.5*pow(-my365 + y365, 2)/pow(sigmay365, 2);
            break;
        case 366:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay366, 2)) + 0.5*pow(-my366 + y366, 2)/pow(sigmay366, 2);
            break;
        case 367:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay367, 2)) + 0.5*pow(-my367 + y367, 2)/pow(sigmay367, 2);
            break;
        case 368:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay368, 2)) + 0.5*pow(-my368 + y368, 2)/pow(sigmay368, 2);
            break;
        case 369:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay369, 2)) + 0.5*pow(-my369 + y369, 2)/pow(sigmay369, 2);
            break;
        case 370:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay370, 2)) + 0.5*pow(-my370 + y370, 2)/pow(sigmay370, 2);
            break;
        case 371:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay371, 2)) + 0.5*pow(-my371 + y371, 2)/pow(sigmay371, 2);
            break;
        case 372:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay372, 2)) + 0.5*pow(-my372 + y372, 2)/pow(sigmay372, 2);
            break;
        case 373:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay373, 2)) + 0.5*pow(-my373 + y373, 2)/pow(sigmay373, 2);
            break;
        case 374:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay374, 2)) + 0.5*pow(-my374 + y374, 2)/pow(sigmay374, 2);
            break;
        case 375:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay375, 2)) + 0.5*pow(-my375 + y375, 2)/pow(sigmay375, 2);
            break;
        case 376:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay376, 2)) + 0.5*pow(-my376 + y376, 2)/pow(sigmay376, 2);
            break;
        case 377:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay377, 2)) + 0.5*pow(-my377 + y377, 2)/pow(sigmay377, 2);
            break;
        case 378:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay378, 2)) + 0.5*pow(-my378 + y378, 2)/pow(sigmay378, 2);
            break;
        case 379:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay379, 2)) + 0.5*pow(-my379 + y379, 2)/pow(sigmay379, 2);
            break;
        case 380:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay380, 2)) + 0.5*pow(-my380 + y380, 2)/pow(sigmay380, 2);
            break;
        case 381:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay381, 2)) + 0.5*pow(-my381 + y381, 2)/pow(sigmay381, 2);
            break;
        case 382:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay382, 2)) + 0.5*pow(-my382 + y382, 2)/pow(sigmay382, 2);
            break;
        case 383:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay383, 2)) + 0.5*pow(-my383 + y383, 2)/pow(sigmay383, 2);
            break;
        case 384:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay384, 2)) + 0.5*pow(-my384 + y384, 2)/pow(sigmay384, 2);
            break;
        case 385:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay385, 2)) + 0.5*pow(-my385 + y385, 2)/pow(sigmay385, 2);
            break;
        case 386:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay386, 2)) + 0.5*pow(-my386 + y386, 2)/pow(sigmay386, 2);
            break;
        case 387:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay387, 2)) + 0.5*pow(-my387 + y387, 2)/pow(sigmay387, 2);
            break;
        case 388:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay388, 2)) + 0.5*pow(-my388 + y388, 2)/pow(sigmay388, 2);
            break;
        case 389:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay389, 2)) + 0.5*pow(-my389 + y389, 2)/pow(sigmay389, 2);
            break;
        case 390:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay390, 2)) + 0.5*pow(-my390 + y390, 2)/pow(sigmay390, 2);
            break;
        case 391:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay391, 2)) + 0.5*pow(-my391 + y391, 2)/pow(sigmay391, 2);
            break;
        case 392:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay392, 2)) + 0.5*pow(-my392 + y392, 2)/pow(sigmay392, 2);
            break;
        case 393:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay393, 2)) + 0.5*pow(-my393 + y393, 2)/pow(sigmay393, 2);
            break;
        case 394:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay394, 2)) + 0.5*pow(-my394 + y394, 2)/pow(sigmay394, 2);
            break;
        case 395:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay395, 2)) + 0.5*pow(-my395 + y395, 2)/pow(sigmay395, 2);
            break;
        case 396:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay396, 2)) + 0.5*pow(-my396 + y396, 2)/pow(sigmay396, 2);
            break;
        case 397:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay397, 2)) + 0.5*pow(-my397 + y397, 2)/pow(sigmay397, 2);
            break;
        case 398:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay398, 2)) + 0.5*pow(-my398 + y398, 2)/pow(sigmay398, 2);
            break;
        case 399:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay399, 2)) + 0.5*pow(-my399 + y399, 2)/pow(sigmay399, 2);
            break;
        case 400:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay400, 2)) + 0.5*pow(-my400 + y400, 2)/pow(sigmay400, 2);
            break;
        case 401:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay401, 2)) + 0.5*pow(-my401 + y401, 2)/pow(sigmay401, 2);
            break;
        case 402:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay402, 2)) + 0.5*pow(-my402 + y402, 2)/pow(sigmay402, 2);
            break;
        case 403:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay403, 2)) + 0.5*pow(-my403 + y403, 2)/pow(sigmay403, 2);
            break;
        case 404:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay404, 2)) + 0.5*pow(-my404 + y404, 2)/pow(sigmay404, 2);
            break;
        case 405:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay405, 2)) + 0.5*pow(-my405 + y405, 2)/pow(sigmay405, 2);
            break;
        case 406:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay406, 2)) + 0.5*pow(-my406 + y406, 2)/pow(sigmay406, 2);
            break;
        case 407:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay407, 2)) + 0.5*pow(-my407 + y407, 2)/pow(sigmay407, 2);
            break;
        case 408:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay408, 2)) + 0.5*pow(-my408 + y408, 2)/pow(sigmay408, 2);
            break;
        case 409:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay409, 2)) + 0.5*pow(-my409 + y409, 2)/pow(sigmay409, 2);
            break;
        case 410:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay410, 2)) + 0.5*pow(-my410 + y410, 2)/pow(sigmay410, 2);
            break;
        case 411:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay411, 2)) + 0.5*pow(-my411 + y411, 2)/pow(sigmay411, 2);
            break;
        case 412:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay412, 2)) + 0.5*pow(-my412 + y412, 2)/pow(sigmay412, 2);
            break;
        case 413:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay413, 2)) + 0.5*pow(-my413 + y413, 2)/pow(sigmay413, 2);
            break;
        case 414:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay414, 2)) + 0.5*pow(-my414 + y414, 2)/pow(sigmay414, 2);
            break;
        case 415:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay415, 2)) + 0.5*pow(-my415 + y415, 2)/pow(sigmay415, 2);
            break;
        case 416:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay416, 2)) + 0.5*pow(-my416 + y416, 2)/pow(sigmay416, 2);
            break;
        case 417:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay417, 2)) + 0.5*pow(-my417 + y417, 2)/pow(sigmay417, 2);
            break;
        case 418:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay418, 2)) + 0.5*pow(-my418 + y418, 2)/pow(sigmay418, 2);
            break;
        case 419:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay419, 2)) + 0.5*pow(-my419 + y419, 2)/pow(sigmay419, 2);
            break;
        case 420:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay420, 2)) + 0.5*pow(-my420 + y420, 2)/pow(sigmay420, 2);
            break;
        case 421:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay421, 2)) + 0.5*pow(-my421 + y421, 2)/pow(sigmay421, 2);
            break;
        case 422:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay422, 2)) + 0.5*pow(-my422 + y422, 2)/pow(sigmay422, 2);
            break;
        case 423:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay423, 2)) + 0.5*pow(-my423 + y423, 2)/pow(sigmay423, 2);
            break;
        case 424:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay424, 2)) + 0.5*pow(-my424 + y424, 2)/pow(sigmay424, 2);
            break;
        case 425:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay425, 2)) + 0.5*pow(-my425 + y425, 2)/pow(sigmay425, 2);
            break;
        case 426:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay426, 2)) + 0.5*pow(-my426 + y426, 2)/pow(sigmay426, 2);
            break;
        case 427:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay427, 2)) + 0.5*pow(-my427 + y427, 2)/pow(sigmay427, 2);
            break;
        case 428:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay428, 2)) + 0.5*pow(-my428 + y428, 2)/pow(sigmay428, 2);
            break;
        case 429:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay429, 2)) + 0.5*pow(-my429 + y429, 2)/pow(sigmay429, 2);
            break;
        case 430:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay430, 2)) + 0.5*pow(-my430 + y430, 2)/pow(sigmay430, 2);
            break;
        case 431:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay431, 2)) + 0.5*pow(-my431 + y431, 2)/pow(sigmay431, 2);
            break;
        case 432:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay432, 2)) + 0.5*pow(-my432 + y432, 2)/pow(sigmay432, 2);
            break;
        case 433:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay433, 2)) + 0.5*pow(-my433 + y433, 2)/pow(sigmay433, 2);
            break;
        case 434:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay434, 2)) + 0.5*pow(-my434 + y434, 2)/pow(sigmay434, 2);
            break;
        case 435:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay435, 2)) + 0.5*pow(-my435 + y435, 2)/pow(sigmay435, 2);
            break;
        case 436:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay436, 2)) + 0.5*pow(-my436 + y436, 2)/pow(sigmay436, 2);
            break;
        case 437:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay437, 2)) + 0.5*pow(-my437 + y437, 2)/pow(sigmay437, 2);
            break;
        case 438:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay438, 2)) + 0.5*pow(-my438 + y438, 2)/pow(sigmay438, 2);
            break;
        case 439:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay439, 2)) + 0.5*pow(-my439 + y439, 2)/pow(sigmay439, 2);
            break;
        case 440:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay440, 2)) + 0.5*pow(-my440 + y440, 2)/pow(sigmay440, 2);
            break;
        case 441:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay441, 2)) + 0.5*pow(-my441 + y441, 2)/pow(sigmay441, 2);
            break;
        case 442:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay442, 2)) + 0.5*pow(-my442 + y442, 2)/pow(sigmay442, 2);
            break;
        case 443:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay443, 2)) + 0.5*pow(-my443 + y443, 2)/pow(sigmay443, 2);
            break;
        case 444:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay444, 2)) + 0.5*pow(-my444 + y444, 2)/pow(sigmay444, 2);
            break;
        case 445:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay445, 2)) + 0.5*pow(-my445 + y445, 2)/pow(sigmay445, 2);
            break;
        case 446:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay446, 2)) + 0.5*pow(-my446 + y446, 2)/pow(sigmay446, 2);
            break;
        case 447:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay447, 2)) + 0.5*pow(-my447 + y447, 2)/pow(sigmay447, 2);
            break;
        case 448:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay448, 2)) + 0.5*pow(-my448 + y448, 2)/pow(sigmay448, 2);
            break;
        case 449:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay449, 2)) + 0.5*pow(-my449 + y449, 2)/pow(sigmay449, 2);
            break;
        case 450:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay450, 2)) + 0.5*pow(-my450 + y450, 2)/pow(sigmay450, 2);
            break;
        case 451:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay451, 2)) + 0.5*pow(-my451 + y451, 2)/pow(sigmay451, 2);
            break;
        case 452:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay452, 2)) + 0.5*pow(-my452 + y452, 2)/pow(sigmay452, 2);
            break;
        case 453:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay453, 2)) + 0.5*pow(-my453 + y453, 2)/pow(sigmay453, 2);
            break;
        case 454:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay454, 2)) + 0.5*pow(-my454 + y454, 2)/pow(sigmay454, 2);
            break;
        case 455:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay455, 2)) + 0.5*pow(-my455 + y455, 2)/pow(sigmay455, 2);
            break;
        case 456:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay456, 2)) + 0.5*pow(-my456 + y456, 2)/pow(sigmay456, 2);
            break;
        case 457:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay457, 2)) + 0.5*pow(-my457 + y457, 2)/pow(sigmay457, 2);
            break;
        case 458:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay458, 2)) + 0.5*pow(-my458 + y458, 2)/pow(sigmay458, 2);
            break;
        case 459:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay459, 2)) + 0.5*pow(-my459 + y459, 2)/pow(sigmay459, 2);
            break;
        case 460:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay460, 2)) + 0.5*pow(-my460 + y460, 2)/pow(sigmay460, 2);
            break;
        case 461:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay461, 2)) + 0.5*pow(-my461 + y461, 2)/pow(sigmay461, 2);
            break;
        case 462:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay462, 2)) + 0.5*pow(-my462 + y462, 2)/pow(sigmay462, 2);
            break;
        case 463:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay463, 2)) + 0.5*pow(-my463 + y463, 2)/pow(sigmay463, 2);
            break;
        case 464:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay464, 2)) + 0.5*pow(-my464 + y464, 2)/pow(sigmay464, 2);
            break;
        case 465:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay465, 2)) + 0.5*pow(-my465 + y465, 2)/pow(sigmay465, 2);
            break;
        case 466:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay466, 2)) + 0.5*pow(-my466 + y466, 2)/pow(sigmay466, 2);
            break;
        case 467:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay467, 2)) + 0.5*pow(-my467 + y467, 2)/pow(sigmay467, 2);
            break;
        case 468:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay468, 2)) + 0.5*pow(-my468 + y468, 2)/pow(sigmay468, 2);
            break;
        case 469:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay469, 2)) + 0.5*pow(-my469 + y469, 2)/pow(sigmay469, 2);
            break;
        case 470:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay470, 2)) + 0.5*pow(-my470 + y470, 2)/pow(sigmay470, 2);
            break;
        case 471:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay471, 2)) + 0.5*pow(-my471 + y471, 2)/pow(sigmay471, 2);
            break;
        case 472:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay472, 2)) + 0.5*pow(-my472 + y472, 2)/pow(sigmay472, 2);
            break;
        case 473:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay473, 2)) + 0.5*pow(-my473 + y473, 2)/pow(sigmay473, 2);
            break;
        case 474:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay474, 2)) + 0.5*pow(-my474 + y474, 2)/pow(sigmay474, 2);
            break;
        case 475:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay475, 2)) + 0.5*pow(-my475 + y475, 2)/pow(sigmay475, 2);
            break;
        case 476:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay476, 2)) + 0.5*pow(-my476 + y476, 2)/pow(sigmay476, 2);
            break;
        case 477:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay477, 2)) + 0.5*pow(-my477 + y477, 2)/pow(sigmay477, 2);
            break;
        case 478:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay478, 2)) + 0.5*pow(-my478 + y478, 2)/pow(sigmay478, 2);
            break;
        case 479:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay479, 2)) + 0.5*pow(-my479 + y479, 2)/pow(sigmay479, 2);
            break;
        case 480:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay480, 2)) + 0.5*pow(-my480 + y480, 2)/pow(sigmay480, 2);
            break;
        case 481:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay481, 2)) + 0.5*pow(-my481 + y481, 2)/pow(sigmay481, 2);
            break;
        case 482:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay482, 2)) + 0.5*pow(-my482 + y482, 2)/pow(sigmay482, 2);
            break;
        case 483:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay483, 2)) + 0.5*pow(-my483 + y483, 2)/pow(sigmay483, 2);
            break;
        case 484:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay484, 2)) + 0.5*pow(-my484 + y484, 2)/pow(sigmay484, 2);
            break;
        case 485:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay485, 2)) + 0.5*pow(-my485 + y485, 2)/pow(sigmay485, 2);
            break;
        case 486:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay486, 2)) + 0.5*pow(-my486 + y486, 2)/pow(sigmay486, 2);
            break;
        case 487:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay487, 2)) + 0.5*pow(-my487 + y487, 2)/pow(sigmay487, 2);
            break;
        case 488:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay488, 2)) + 0.5*pow(-my488 + y488, 2)/pow(sigmay488, 2);
            break;
        case 489:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay489, 2)) + 0.5*pow(-my489 + y489, 2)/pow(sigmay489, 2);
            break;
        case 490:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay490, 2)) + 0.5*pow(-my490 + y490, 2)/pow(sigmay490, 2);
            break;
        case 491:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay491, 2)) + 0.5*pow(-my491 + y491, 2)/pow(sigmay491, 2);
            break;
        case 492:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay492, 2)) + 0.5*pow(-my492 + y492, 2)/pow(sigmay492, 2);
            break;
        case 493:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay493, 2)) + 0.5*pow(-my493 + y493, 2)/pow(sigmay493, 2);
            break;
        case 494:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay494, 2)) + 0.5*pow(-my494 + y494, 2)/pow(sigmay494, 2);
            break;
        case 495:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay495, 2)) + 0.5*pow(-my495 + y495, 2)/pow(sigmay495, 2);
            break;
        case 496:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay496, 2)) + 0.5*pow(-my496 + y496, 2)/pow(sigmay496, 2);
            break;
        case 497:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay497, 2)) + 0.5*pow(-my497 + y497, 2)/pow(sigmay497, 2);
            break;
        case 498:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay498, 2)) + 0.5*pow(-my498 + y498, 2)/pow(sigmay498, 2);
            break;
        case 499:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay499, 2)) + 0.5*pow(-my499 + y499, 2)/pow(sigmay499, 2);
            break;
        case 500:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay500, 2)) + 0.5*pow(-my500 + y500, 2)/pow(sigmay500, 2);
            break;
        case 501:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay501, 2)) + 0.5*pow(-my501 + y501, 2)/pow(sigmay501, 2);
            break;
        case 502:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay502, 2)) + 0.5*pow(-my502 + y502, 2)/pow(sigmay502, 2);
            break;
        case 503:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay503, 2)) + 0.5*pow(-my503 + y503, 2)/pow(sigmay503, 2);
            break;
        case 504:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay504, 2)) + 0.5*pow(-my504 + y504, 2)/pow(sigmay504, 2);
            break;
        case 505:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay505, 2)) + 0.5*pow(-my505 + y505, 2)/pow(sigmay505, 2);
            break;
        case 506:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay506, 2)) + 0.5*pow(-my506 + y506, 2)/pow(sigmay506, 2);
            break;
        case 507:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay507, 2)) + 0.5*pow(-my507 + y507, 2)/pow(sigmay507, 2);
            break;
        case 508:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay508, 2)) + 0.5*pow(-my508 + y508, 2)/pow(sigmay508, 2);
            break;
        case 509:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay509, 2)) + 0.5*pow(-my509 + y509, 2)/pow(sigmay509, 2);
            break;
        case 510:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay510, 2)) + 0.5*pow(-my510 + y510, 2)/pow(sigmay510, 2);
            break;
        case 511:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay511, 2)) + 0.5*pow(-my511 + y511, 2)/pow(sigmay511, 2);
            break;
        case 512:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay512, 2)) + 0.5*pow(-my512 + y512, 2)/pow(sigmay512, 2);
            break;
        case 513:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay513, 2)) + 0.5*pow(-my513 + y513, 2)/pow(sigmay513, 2);
            break;
        case 514:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay514, 2)) + 0.5*pow(-my514 + y514, 2)/pow(sigmay514, 2);
            break;
        case 515:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay515, 2)) + 0.5*pow(-my515 + y515, 2)/pow(sigmay515, 2);
            break;
        case 516:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay516, 2)) + 0.5*pow(-my516 + y516, 2)/pow(sigmay516, 2);
            break;
        case 517:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay517, 2)) + 0.5*pow(-my517 + y517, 2)/pow(sigmay517, 2);
            break;
        case 518:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay518, 2)) + 0.5*pow(-my518 + y518, 2)/pow(sigmay518, 2);
            break;
        case 519:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay519, 2)) + 0.5*pow(-my519 + y519, 2)/pow(sigmay519, 2);
            break;
        case 520:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay520, 2)) + 0.5*pow(-my520 + y520, 2)/pow(sigmay520, 2);
            break;
        case 521:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay521, 2)) + 0.5*pow(-my521 + y521, 2)/pow(sigmay521, 2);
            break;
        case 522:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay522, 2)) + 0.5*pow(-my522 + y522, 2)/pow(sigmay522, 2);
            break;
        case 523:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay523, 2)) + 0.5*pow(-my523 + y523, 2)/pow(sigmay523, 2);
            break;
        case 524:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay524, 2)) + 0.5*pow(-my524 + y524, 2)/pow(sigmay524, 2);
            break;
        case 525:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay525, 2)) + 0.5*pow(-my525 + y525, 2)/pow(sigmay525, 2);
            break;
        case 526:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay526, 2)) + 0.5*pow(-my526 + y526, 2)/pow(sigmay526, 2);
            break;
        case 527:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay527, 2)) + 0.5*pow(-my527 + y527, 2)/pow(sigmay527, 2);
            break;
        case 528:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay528, 2)) + 0.5*pow(-my528 + y528, 2)/pow(sigmay528, 2);
            break;
        case 529:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay529, 2)) + 0.5*pow(-my529 + y529, 2)/pow(sigmay529, 2);
            break;
        case 530:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay530, 2)) + 0.5*pow(-my530 + y530, 2)/pow(sigmay530, 2);
            break;
        case 531:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay531, 2)) + 0.5*pow(-my531 + y531, 2)/pow(sigmay531, 2);
            break;
        case 532:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay532, 2)) + 0.5*pow(-my532 + y532, 2)/pow(sigmay532, 2);
            break;
        case 533:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay533, 2)) + 0.5*pow(-my533 + y533, 2)/pow(sigmay533, 2);
            break;
        case 534:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay534, 2)) + 0.5*pow(-my534 + y534, 2)/pow(sigmay534, 2);
            break;
        case 535:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay535, 2)) + 0.5*pow(-my535 + y535, 2)/pow(sigmay535, 2);
            break;
        case 536:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay536, 2)) + 0.5*pow(-my536 + y536, 2)/pow(sigmay536, 2);
            break;
        case 537:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay537, 2)) + 0.5*pow(-my537 + y537, 2)/pow(sigmay537, 2);
            break;
        case 538:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay538, 2)) + 0.5*pow(-my538 + y538, 2)/pow(sigmay538, 2);
            break;
        case 539:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay539, 2)) + 0.5*pow(-my539 + y539, 2)/pow(sigmay539, 2);
            break;
        case 540:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay540, 2)) + 0.5*pow(-my540 + y540, 2)/pow(sigmay540, 2);
            break;
        case 541:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay541, 2)) + 0.5*pow(-my541 + y541, 2)/pow(sigmay541, 2);
            break;
        case 542:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay542, 2)) + 0.5*pow(-my542 + y542, 2)/pow(sigmay542, 2);
            break;
        case 543:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay543, 2)) + 0.5*pow(-my543 + y543, 2)/pow(sigmay543, 2);
            break;
        case 544:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay544, 2)) + 0.5*pow(-my544 + y544, 2)/pow(sigmay544, 2);
            break;
        case 545:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay545, 2)) + 0.5*pow(-my545 + y545, 2)/pow(sigmay545, 2);
            break;
        case 546:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay546, 2)) + 0.5*pow(-my546 + y546, 2)/pow(sigmay546, 2);
            break;
        case 547:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay547, 2)) + 0.5*pow(-my547 + y547, 2)/pow(sigmay547, 2);
            break;
        case 548:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay548, 2)) + 0.5*pow(-my548 + y548, 2)/pow(sigmay548, 2);
            break;
        case 549:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay549, 2)) + 0.5*pow(-my549 + y549, 2)/pow(sigmay549, 2);
            break;
        case 550:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay550, 2)) + 0.5*pow(-my550 + y550, 2)/pow(sigmay550, 2);
            break;
        case 551:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay551, 2)) + 0.5*pow(-my551 + y551, 2)/pow(sigmay551, 2);
            break;
        case 552:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay552, 2)) + 0.5*pow(-my552 + y552, 2)/pow(sigmay552, 2);
            break;
        case 553:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay553, 2)) + 0.5*pow(-my553 + y553, 2)/pow(sigmay553, 2);
            break;
        case 554:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay554, 2)) + 0.5*pow(-my554 + y554, 2)/pow(sigmay554, 2);
            break;
        case 555:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay555, 2)) + 0.5*pow(-my555 + y555, 2)/pow(sigmay555, 2);
            break;
        case 556:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay556, 2)) + 0.5*pow(-my556 + y556, 2)/pow(sigmay556, 2);
            break;
        case 557:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay557, 2)) + 0.5*pow(-my557 + y557, 2)/pow(sigmay557, 2);
            break;
        case 558:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay558, 2)) + 0.5*pow(-my558 + y558, 2)/pow(sigmay558, 2);
            break;
        case 559:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay559, 2)) + 0.5*pow(-my559 + y559, 2)/pow(sigmay559, 2);
            break;
        case 560:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay560, 2)) + 0.5*pow(-my560 + y560, 2)/pow(sigmay560, 2);
            break;
        case 561:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay561, 2)) + 0.5*pow(-my561 + y561, 2)/pow(sigmay561, 2);
            break;
        case 562:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay562, 2)) + 0.5*pow(-my562 + y562, 2)/pow(sigmay562, 2);
            break;
        case 563:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay563, 2)) + 0.5*pow(-my563 + y563, 2)/pow(sigmay563, 2);
            break;
        case 564:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay564, 2)) + 0.5*pow(-my564 + y564, 2)/pow(sigmay564, 2);
            break;
        case 565:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay565, 2)) + 0.5*pow(-my565 + y565, 2)/pow(sigmay565, 2);
            break;
        case 566:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay566, 2)) + 0.5*pow(-my566 + y566, 2)/pow(sigmay566, 2);
            break;
        case 567:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay567, 2)) + 0.5*pow(-my567 + y567, 2)/pow(sigmay567, 2);
            break;
        case 568:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay568, 2)) + 0.5*pow(-my568 + y568, 2)/pow(sigmay568, 2);
            break;
        case 569:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay569, 2)) + 0.5*pow(-my569 + y569, 2)/pow(sigmay569, 2);
            break;
        case 570:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay570, 2)) + 0.5*pow(-my570 + y570, 2)/pow(sigmay570, 2);
            break;
        case 571:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay571, 2)) + 0.5*pow(-my571 + y571, 2)/pow(sigmay571, 2);
            break;
        case 572:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay572, 2)) + 0.5*pow(-my572 + y572, 2)/pow(sigmay572, 2);
            break;
        case 573:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay573, 2)) + 0.5*pow(-my573 + y573, 2)/pow(sigmay573, 2);
            break;
        case 574:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay574, 2)) + 0.5*pow(-my574 + y574, 2)/pow(sigmay574, 2);
            break;
        case 575:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay575, 2)) + 0.5*pow(-my575 + y575, 2)/pow(sigmay575, 2);
            break;
        case 576:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay576, 2)) + 0.5*pow(-my576 + y576, 2)/pow(sigmay576, 2);
            break;
        case 577:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay577, 2)) + 0.5*pow(-my577 + y577, 2)/pow(sigmay577, 2);
            break;
        case 578:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay578, 2)) + 0.5*pow(-my578 + y578, 2)/pow(sigmay578, 2);
            break;
        case 579:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay579, 2)) + 0.5*pow(-my579 + y579, 2)/pow(sigmay579, 2);
            break;
        case 580:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay580, 2)) + 0.5*pow(-my580 + y580, 2)/pow(sigmay580, 2);
            break;
        case 581:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay581, 2)) + 0.5*pow(-my581 + y581, 2)/pow(sigmay581, 2);
            break;
        case 582:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay582, 2)) + 0.5*pow(-my582 + y582, 2)/pow(sigmay582, 2);
            break;
        case 583:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay583, 2)) + 0.5*pow(-my583 + y583, 2)/pow(sigmay583, 2);
            break;
        case 584:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay584, 2)) + 0.5*pow(-my584 + y584, 2)/pow(sigmay584, 2);
            break;
        case 585:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay585, 2)) + 0.5*pow(-my585 + y585, 2)/pow(sigmay585, 2);
            break;
        case 586:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay586, 2)) + 0.5*pow(-my586 + y586, 2)/pow(sigmay586, 2);
            break;
        case 587:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay587, 2)) + 0.5*pow(-my587 + y587, 2)/pow(sigmay587, 2);
            break;
        case 588:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay588, 2)) + 0.5*pow(-my588 + y588, 2)/pow(sigmay588, 2);
            break;
        case 589:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay589, 2)) + 0.5*pow(-my589 + y589, 2)/pow(sigmay589, 2);
            break;
        case 590:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay590, 2)) + 0.5*pow(-my590 + y590, 2)/pow(sigmay590, 2);
            break;
        case 591:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay591, 2)) + 0.5*pow(-my591 + y591, 2)/pow(sigmay591, 2);
            break;
        case 592:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay592, 2)) + 0.5*pow(-my592 + y592, 2)/pow(sigmay592, 2);
            break;
        case 593:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay593, 2)) + 0.5*pow(-my593 + y593, 2)/pow(sigmay593, 2);
            break;
        case 594:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay594, 2)) + 0.5*pow(-my594 + y594, 2)/pow(sigmay594, 2);
            break;
        case 595:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay595, 2)) + 0.5*pow(-my595 + y595, 2)/pow(sigmay595, 2);
            break;
        case 596:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay596, 2)) + 0.5*pow(-my596 + y596, 2)/pow(sigmay596, 2);
            break;
        case 597:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay597, 2)) + 0.5*pow(-my597 + y597, 2)/pow(sigmay597, 2);
            break;
        case 598:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay598, 2)) + 0.5*pow(-my598 + y598, 2)/pow(sigmay598, 2);
            break;
        case 599:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay599, 2)) + 0.5*pow(-my599 + y599, 2)/pow(sigmay599, 2);
            break;
        case 600:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay600, 2)) + 0.5*pow(-my600 + y600, 2)/pow(sigmay600, 2);
            break;
        case 601:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay601, 2)) + 0.5*pow(-my601 + y601, 2)/pow(sigmay601, 2);
            break;
        case 602:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay602, 2)) + 0.5*pow(-my602 + y602, 2)/pow(sigmay602, 2);
            break;
        case 603:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay603, 2)) + 0.5*pow(-my603 + y603, 2)/pow(sigmay603, 2);
            break;
        case 604:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay604, 2)) + 0.5*pow(-my604 + y604, 2)/pow(sigmay604, 2);
            break;
        case 605:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay605, 2)) + 0.5*pow(-my605 + y605, 2)/pow(sigmay605, 2);
            break;
        case 606:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay606, 2)) + 0.5*pow(-my606 + y606, 2)/pow(sigmay606, 2);
            break;
        case 607:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay607, 2)) + 0.5*pow(-my607 + y607, 2)/pow(sigmay607, 2);
            break;
        case 608:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay608, 2)) + 0.5*pow(-my608 + y608, 2)/pow(sigmay608, 2);
            break;
        case 609:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay609, 2)) + 0.5*pow(-my609 + y609, 2)/pow(sigmay609, 2);
            break;
        case 610:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay610, 2)) + 0.5*pow(-my610 + y610, 2)/pow(sigmay610, 2);
            break;
        case 611:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay611, 2)) + 0.5*pow(-my611 + y611, 2)/pow(sigmay611, 2);
            break;
        case 612:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay612, 2)) + 0.5*pow(-my612 + y612, 2)/pow(sigmay612, 2);
            break;
        case 613:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay613, 2)) + 0.5*pow(-my613 + y613, 2)/pow(sigmay613, 2);
            break;
        case 614:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay614, 2)) + 0.5*pow(-my614 + y614, 2)/pow(sigmay614, 2);
            break;
        case 615:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay615, 2)) + 0.5*pow(-my615 + y615, 2)/pow(sigmay615, 2);
            break;
        case 616:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay616, 2)) + 0.5*pow(-my616 + y616, 2)/pow(sigmay616, 2);
            break;
        case 617:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay617, 2)) + 0.5*pow(-my617 + y617, 2)/pow(sigmay617, 2);
            break;
        case 618:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay618, 2)) + 0.5*pow(-my618 + y618, 2)/pow(sigmay618, 2);
            break;
        case 619:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay619, 2)) + 0.5*pow(-my619 + y619, 2)/pow(sigmay619, 2);
            break;
        case 620:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay620, 2)) + 0.5*pow(-my620 + y620, 2)/pow(sigmay620, 2);
            break;
        case 621:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay621, 2)) + 0.5*pow(-my621 + y621, 2)/pow(sigmay621, 2);
            break;
        case 622:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay622, 2)) + 0.5*pow(-my622 + y622, 2)/pow(sigmay622, 2);
            break;
        case 623:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay623, 2)) + 0.5*pow(-my623 + y623, 2)/pow(sigmay623, 2);
            break;
        case 624:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay624, 2)) + 0.5*pow(-my624 + y624, 2)/pow(sigmay624, 2);
            break;
        case 625:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay625, 2)) + 0.5*pow(-my625 + y625, 2)/pow(sigmay625, 2);
            break;
        case 626:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay626, 2)) + 0.5*pow(-my626 + y626, 2)/pow(sigmay626, 2);
            break;
        case 627:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay627, 2)) + 0.5*pow(-my627 + y627, 2)/pow(sigmay627, 2);
            break;
        case 628:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay628, 2)) + 0.5*pow(-my628 + y628, 2)/pow(sigmay628, 2);
            break;
        case 629:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay629, 2)) + 0.5*pow(-my629 + y629, 2)/pow(sigmay629, 2);
            break;
        case 630:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay630, 2)) + 0.5*pow(-my630 + y630, 2)/pow(sigmay630, 2);
            break;
        case 631:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay631, 2)) + 0.5*pow(-my631 + y631, 2)/pow(sigmay631, 2);
            break;
        case 632:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay632, 2)) + 0.5*pow(-my632 + y632, 2)/pow(sigmay632, 2);
            break;
        case 633:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay633, 2)) + 0.5*pow(-my633 + y633, 2)/pow(sigmay633, 2);
            break;
        case 634:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay634, 2)) + 0.5*pow(-my634 + y634, 2)/pow(sigmay634, 2);
            break;
        case 635:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay635, 2)) + 0.5*pow(-my635 + y635, 2)/pow(sigmay635, 2);
            break;
        case 636:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay636, 2)) + 0.5*pow(-my636 + y636, 2)/pow(sigmay636, 2);
            break;
        case 637:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay637, 2)) + 0.5*pow(-my637 + y637, 2)/pow(sigmay637, 2);
            break;
        case 638:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay638, 2)) + 0.5*pow(-my638 + y638, 2)/pow(sigmay638, 2);
            break;
        case 639:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay639, 2)) + 0.5*pow(-my639 + y639, 2)/pow(sigmay639, 2);
            break;
        case 640:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay640, 2)) + 0.5*pow(-my640 + y640, 2)/pow(sigmay640, 2);
            break;
        case 641:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay641, 2)) + 0.5*pow(-my641 + y641, 2)/pow(sigmay641, 2);
            break;
        case 642:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay642, 2)) + 0.5*pow(-my642 + y642, 2)/pow(sigmay642, 2);
            break;
        case 643:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay643, 2)) + 0.5*pow(-my643 + y643, 2)/pow(sigmay643, 2);
            break;
        case 644:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay644, 2)) + 0.5*pow(-my644 + y644, 2)/pow(sigmay644, 2);
            break;
        case 645:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay645, 2)) + 0.5*pow(-my645 + y645, 2)/pow(sigmay645, 2);
            break;
        case 646:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay646, 2)) + 0.5*pow(-my646 + y646, 2)/pow(sigmay646, 2);
            break;
        case 647:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay647, 2)) + 0.5*pow(-my647 + y647, 2)/pow(sigmay647, 2);
            break;
        case 648:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay648, 2)) + 0.5*pow(-my648 + y648, 2)/pow(sigmay648, 2);
            break;
        case 649:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay649, 2)) + 0.5*pow(-my649 + y649, 2)/pow(sigmay649, 2);
            break;
        case 650:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay650, 2)) + 0.5*pow(-my650 + y650, 2)/pow(sigmay650, 2);
            break;
        case 651:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay651, 2)) + 0.5*pow(-my651 + y651, 2)/pow(sigmay651, 2);
            break;
        case 652:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay652, 2)) + 0.5*pow(-my652 + y652, 2)/pow(sigmay652, 2);
            break;
        case 653:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay653, 2)) + 0.5*pow(-my653 + y653, 2)/pow(sigmay653, 2);
            break;
        case 654:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay654, 2)) + 0.5*pow(-my654 + y654, 2)/pow(sigmay654, 2);
            break;
        case 655:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay655, 2)) + 0.5*pow(-my655 + y655, 2)/pow(sigmay655, 2);
            break;
        case 656:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay656, 2)) + 0.5*pow(-my656 + y656, 2)/pow(sigmay656, 2);
            break;
        case 657:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay657, 2)) + 0.5*pow(-my657 + y657, 2)/pow(sigmay657, 2);
            break;
        case 658:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay658, 2)) + 0.5*pow(-my658 + y658, 2)/pow(sigmay658, 2);
            break;
        case 659:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay659, 2)) + 0.5*pow(-my659 + y659, 2)/pow(sigmay659, 2);
            break;
        case 660:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay660, 2)) + 0.5*pow(-my660 + y660, 2)/pow(sigmay660, 2);
            break;
        case 661:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay661, 2)) + 0.5*pow(-my661 + y661, 2)/pow(sigmay661, 2);
            break;
        case 662:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay662, 2)) + 0.5*pow(-my662 + y662, 2)/pow(sigmay662, 2);
            break;
        case 663:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay663, 2)) + 0.5*pow(-my663 + y663, 2)/pow(sigmay663, 2);
            break;
        case 664:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay664, 2)) + 0.5*pow(-my664 + y664, 2)/pow(sigmay664, 2);
            break;
        case 665:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay665, 2)) + 0.5*pow(-my665 + y665, 2)/pow(sigmay665, 2);
            break;
        case 666:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay666, 2)) + 0.5*pow(-my666 + y666, 2)/pow(sigmay666, 2);
            break;
        case 667:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay667, 2)) + 0.5*pow(-my667 + y667, 2)/pow(sigmay667, 2);
            break;
        case 668:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay668, 2)) + 0.5*pow(-my668 + y668, 2)/pow(sigmay668, 2);
            break;
        case 669:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay669, 2)) + 0.5*pow(-my669 + y669, 2)/pow(sigmay669, 2);
            break;
        case 670:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay670, 2)) + 0.5*pow(-my670 + y670, 2)/pow(sigmay670, 2);
            break;
        case 671:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay671, 2)) + 0.5*pow(-my671 + y671, 2)/pow(sigmay671, 2);
            break;
        case 672:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay672, 2)) + 0.5*pow(-my672 + y672, 2)/pow(sigmay672, 2);
            break;
        case 673:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay673, 2)) + 0.5*pow(-my673 + y673, 2)/pow(sigmay673, 2);
            break;
        case 674:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay674, 2)) + 0.5*pow(-my674 + y674, 2)/pow(sigmay674, 2);
            break;
        case 675:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay675, 2)) + 0.5*pow(-my675 + y675, 2)/pow(sigmay675, 2);
            break;
        case 676:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay676, 2)) + 0.5*pow(-my676 + y676, 2)/pow(sigmay676, 2);
            break;
        case 677:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay677, 2)) + 0.5*pow(-my677 + y677, 2)/pow(sigmay677, 2);
            break;
        case 678:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay678, 2)) + 0.5*pow(-my678 + y678, 2)/pow(sigmay678, 2);
            break;
        case 679:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay679, 2)) + 0.5*pow(-my679 + y679, 2)/pow(sigmay679, 2);
            break;
        case 680:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay680, 2)) + 0.5*pow(-my680 + y680, 2)/pow(sigmay680, 2);
            break;
        case 681:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay681, 2)) + 0.5*pow(-my681 + y681, 2)/pow(sigmay681, 2);
            break;
        case 682:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay682, 2)) + 0.5*pow(-my682 + y682, 2)/pow(sigmay682, 2);
            break;
        case 683:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay683, 2)) + 0.5*pow(-my683 + y683, 2)/pow(sigmay683, 2);
            break;
        case 684:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay684, 2)) + 0.5*pow(-my684 + y684, 2)/pow(sigmay684, 2);
            break;
        case 685:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay685, 2)) + 0.5*pow(-my685 + y685, 2)/pow(sigmay685, 2);
            break;
        case 686:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay686, 2)) + 0.5*pow(-my686 + y686, 2)/pow(sigmay686, 2);
            break;
        case 687:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay687, 2)) + 0.5*pow(-my687 + y687, 2)/pow(sigmay687, 2);
            break;
        case 688:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay688, 2)) + 0.5*pow(-my688 + y688, 2)/pow(sigmay688, 2);
            break;
        case 689:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay689, 2)) + 0.5*pow(-my689 + y689, 2)/pow(sigmay689, 2);
            break;
        case 690:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay690, 2)) + 0.5*pow(-my690 + y690, 2)/pow(sigmay690, 2);
            break;
        case 691:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay691, 2)) + 0.5*pow(-my691 + y691, 2)/pow(sigmay691, 2);
            break;
        case 692:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay692, 2)) + 0.5*pow(-my692 + y692, 2)/pow(sigmay692, 2);
            break;
        case 693:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay693, 2)) + 0.5*pow(-my693 + y693, 2)/pow(sigmay693, 2);
            break;
        case 694:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay694, 2)) + 0.5*pow(-my694 + y694, 2)/pow(sigmay694, 2);
            break;
        case 695:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay695, 2)) + 0.5*pow(-my695 + y695, 2)/pow(sigmay695, 2);
            break;
        case 696:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay696, 2)) + 0.5*pow(-my696 + y696, 2)/pow(sigmay696, 2);
            break;
        case 697:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay697, 2)) + 0.5*pow(-my697 + y697, 2)/pow(sigmay697, 2);
            break;
        case 698:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay698, 2)) + 0.5*pow(-my698 + y698, 2)/pow(sigmay698, 2);
            break;
        case 699:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay699, 2)) + 0.5*pow(-my699 + y699, 2)/pow(sigmay699, 2);
            break;
        case 700:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay700, 2)) + 0.5*pow(-my700 + y700, 2)/pow(sigmay700, 2);
            break;
        case 701:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay701, 2)) + 0.5*pow(-my701 + y701, 2)/pow(sigmay701, 2);
            break;
        case 702:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay702, 2)) + 0.5*pow(-my702 + y702, 2)/pow(sigmay702, 2);
            break;
        case 703:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay703, 2)) + 0.5*pow(-my703 + y703, 2)/pow(sigmay703, 2);
            break;
        case 704:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay704, 2)) + 0.5*pow(-my704 + y704, 2)/pow(sigmay704, 2);
            break;
        case 705:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay705, 2)) + 0.5*pow(-my705 + y705, 2)/pow(sigmay705, 2);
            break;
        case 706:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay706, 2)) + 0.5*pow(-my706 + y706, 2)/pow(sigmay706, 2);
            break;
        case 707:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay707, 2)) + 0.5*pow(-my707 + y707, 2)/pow(sigmay707, 2);
            break;
        case 708:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay708, 2)) + 0.5*pow(-my708 + y708, 2)/pow(sigmay708, 2);
            break;
        case 709:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay709, 2)) + 0.5*pow(-my709 + y709, 2)/pow(sigmay709, 2);
            break;
        case 710:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay710, 2)) + 0.5*pow(-my710 + y710, 2)/pow(sigmay710, 2);
            break;
        case 711:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay711, 2)) + 0.5*pow(-my711 + y711, 2)/pow(sigmay711, 2);
            break;
        case 712:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay712, 2)) + 0.5*pow(-my712 + y712, 2)/pow(sigmay712, 2);
            break;
        case 713:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay713, 2)) + 0.5*pow(-my713 + y713, 2)/pow(sigmay713, 2);
            break;
        case 714:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay714, 2)) + 0.5*pow(-my714 + y714, 2)/pow(sigmay714, 2);
            break;
        case 715:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay715, 2)) + 0.5*pow(-my715 + y715, 2)/pow(sigmay715, 2);
            break;
        case 716:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay716, 2)) + 0.5*pow(-my716 + y716, 2)/pow(sigmay716, 2);
            break;
        case 717:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay717, 2)) + 0.5*pow(-my717 + y717, 2)/pow(sigmay717, 2);
            break;
        case 718:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay718, 2)) + 0.5*pow(-my718 + y718, 2)/pow(sigmay718, 2);
            break;
        case 719:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay719, 2)) + 0.5*pow(-my719 + y719, 2)/pow(sigmay719, 2);
            break;
        case 720:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay720, 2)) + 0.5*pow(-my720 + y720, 2)/pow(sigmay720, 2);
            break;
        case 721:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay721, 2)) + 0.5*pow(-my721 + y721, 2)/pow(sigmay721, 2);
            break;
        case 722:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay722, 2)) + 0.5*pow(-my722 + y722, 2)/pow(sigmay722, 2);
            break;
        case 723:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay723, 2)) + 0.5*pow(-my723 + y723, 2)/pow(sigmay723, 2);
            break;
        case 724:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay724, 2)) + 0.5*pow(-my724 + y724, 2)/pow(sigmay724, 2);
            break;
        case 725:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay725, 2)) + 0.5*pow(-my725 + y725, 2)/pow(sigmay725, 2);
            break;
        case 726:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay726, 2)) + 0.5*pow(-my726 + y726, 2)/pow(sigmay726, 2);
            break;
        case 727:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay727, 2)) + 0.5*pow(-my727 + y727, 2)/pow(sigmay727, 2);
            break;
        case 728:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay728, 2)) + 0.5*pow(-my728 + y728, 2)/pow(sigmay728, 2);
            break;
        case 729:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay729, 2)) + 0.5*pow(-my729 + y729, 2)/pow(sigmay729, 2);
            break;
        case 730:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay730, 2)) + 0.5*pow(-my730 + y730, 2)/pow(sigmay730, 2);
            break;
        case 731:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay731, 2)) + 0.5*pow(-my731 + y731, 2)/pow(sigmay731, 2);
            break;
        case 732:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay732, 2)) + 0.5*pow(-my732 + y732, 2)/pow(sigmay732, 2);
            break;
        case 733:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay733, 2)) + 0.5*pow(-my733 + y733, 2)/pow(sigmay733, 2);
            break;
        case 734:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay734, 2)) + 0.5*pow(-my734 + y734, 2)/pow(sigmay734, 2);
            break;
        case 735:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay735, 2)) + 0.5*pow(-my735 + y735, 2)/pow(sigmay735, 2);
            break;
        case 736:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay736, 2)) + 0.5*pow(-my736 + y736, 2)/pow(sigmay736, 2);
            break;
        case 737:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay737, 2)) + 0.5*pow(-my737 + y737, 2)/pow(sigmay737, 2);
            break;
        case 738:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay738, 2)) + 0.5*pow(-my738 + y738, 2)/pow(sigmay738, 2);
            break;
        case 739:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay739, 2)) + 0.5*pow(-my739 + y739, 2)/pow(sigmay739, 2);
            break;
        case 740:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay740, 2)) + 0.5*pow(-my740 + y740, 2)/pow(sigmay740, 2);
            break;
        case 741:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay741, 2)) + 0.5*pow(-my741 + y741, 2)/pow(sigmay741, 2);
            break;
        case 742:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay742, 2)) + 0.5*pow(-my742 + y742, 2)/pow(sigmay742, 2);
            break;
        case 743:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay743, 2)) + 0.5*pow(-my743 + y743, 2)/pow(sigmay743, 2);
            break;
        case 744:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay744, 2)) + 0.5*pow(-my744 + y744, 2)/pow(sigmay744, 2);
            break;
        case 745:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay745, 2)) + 0.5*pow(-my745 + y745, 2)/pow(sigmay745, 2);
            break;
        case 746:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay746, 2)) + 0.5*pow(-my746 + y746, 2)/pow(sigmay746, 2);
            break;
        case 747:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay747, 2)) + 0.5*pow(-my747 + y747, 2)/pow(sigmay747, 2);
            break;
        case 748:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay748, 2)) + 0.5*pow(-my748 + y748, 2)/pow(sigmay748, 2);
            break;
        case 749:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay749, 2)) + 0.5*pow(-my749 + y749, 2)/pow(sigmay749, 2);
            break;
        case 750:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay750, 2)) + 0.5*pow(-my750 + y750, 2)/pow(sigmay750, 2);
            break;
        case 751:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay751, 2)) + 0.5*pow(-my751 + y751, 2)/pow(sigmay751, 2);
            break;
        case 752:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay752, 2)) + 0.5*pow(-my752 + y752, 2)/pow(sigmay752, 2);
            break;
        case 753:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay753, 2)) + 0.5*pow(-my753 + y753, 2)/pow(sigmay753, 2);
            break;
        case 754:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay754, 2)) + 0.5*pow(-my754 + y754, 2)/pow(sigmay754, 2);
            break;
        case 755:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay755, 2)) + 0.5*pow(-my755 + y755, 2)/pow(sigmay755, 2);
            break;
        case 756:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay756, 2)) + 0.5*pow(-my756 + y756, 2)/pow(sigmay756, 2);
            break;
        case 757:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay757, 2)) + 0.5*pow(-my757 + y757, 2)/pow(sigmay757, 2);
            break;
        case 758:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay758, 2)) + 0.5*pow(-my758 + y758, 2)/pow(sigmay758, 2);
            break;
        case 759:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay759, 2)) + 0.5*pow(-my759 + y759, 2)/pow(sigmay759, 2);
            break;
        case 760:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay760, 2)) + 0.5*pow(-my760 + y760, 2)/pow(sigmay760, 2);
            break;
        case 761:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay761, 2)) + 0.5*pow(-my761 + y761, 2)/pow(sigmay761, 2);
            break;
        case 762:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay762, 2)) + 0.5*pow(-my762 + y762, 2)/pow(sigmay762, 2);
            break;
        case 763:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay763, 2)) + 0.5*pow(-my763 + y763, 2)/pow(sigmay763, 2);
            break;
        case 764:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay764, 2)) + 0.5*pow(-my764 + y764, 2)/pow(sigmay764, 2);
            break;
        case 765:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay765, 2)) + 0.5*pow(-my765 + y765, 2)/pow(sigmay765, 2);
            break;
        case 766:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay766, 2)) + 0.5*pow(-my766 + y766, 2)/pow(sigmay766, 2);
            break;
        case 767:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay767, 2)) + 0.5*pow(-my767 + y767, 2)/pow(sigmay767, 2);
            break;
        case 768:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay768, 2)) + 0.5*pow(-my768 + y768, 2)/pow(sigmay768, 2);
            break;
        case 769:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay769, 2)) + 0.5*pow(-my769 + y769, 2)/pow(sigmay769, 2);
            break;
        case 770:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay770, 2)) + 0.5*pow(-my770 + y770, 2)/pow(sigmay770, 2);
            break;
        case 771:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay771, 2)) + 0.5*pow(-my771 + y771, 2)/pow(sigmay771, 2);
            break;
        case 772:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay772, 2)) + 0.5*pow(-my772 + y772, 2)/pow(sigmay772, 2);
            break;
        case 773:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay773, 2)) + 0.5*pow(-my773 + y773, 2)/pow(sigmay773, 2);
            break;
        case 774:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay774, 2)) + 0.5*pow(-my774 + y774, 2)/pow(sigmay774, 2);
            break;
        case 775:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay775, 2)) + 0.5*pow(-my775 + y775, 2)/pow(sigmay775, 2);
            break;
        case 776:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay776, 2)) + 0.5*pow(-my776 + y776, 2)/pow(sigmay776, 2);
            break;
        case 777:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay777, 2)) + 0.5*pow(-my777 + y777, 2)/pow(sigmay777, 2);
            break;
        case 778:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay778, 2)) + 0.5*pow(-my778 + y778, 2)/pow(sigmay778, 2);
            break;
        case 779:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay779, 2)) + 0.5*pow(-my779 + y779, 2)/pow(sigmay779, 2);
            break;
        case 780:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay780, 2)) + 0.5*pow(-my780 + y780, 2)/pow(sigmay780, 2);
            break;
        case 781:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay781, 2)) + 0.5*pow(-my781 + y781, 2)/pow(sigmay781, 2);
            break;
        case 782:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay782, 2)) + 0.5*pow(-my782 + y782, 2)/pow(sigmay782, 2);
            break;
        case 783:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay783, 2)) + 0.5*pow(-my783 + y783, 2)/pow(sigmay783, 2);
            break;
        case 784:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay784, 2)) + 0.5*pow(-my784 + y784, 2)/pow(sigmay784, 2);
            break;
        case 785:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay785, 2)) + 0.5*pow(-my785 + y785, 2)/pow(sigmay785, 2);
            break;
        case 786:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay786, 2)) + 0.5*pow(-my786 + y786, 2)/pow(sigmay786, 2);
            break;
        case 787:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay787, 2)) + 0.5*pow(-my787 + y787, 2)/pow(sigmay787, 2);
            break;
        case 788:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay788, 2)) + 0.5*pow(-my788 + y788, 2)/pow(sigmay788, 2);
            break;
        case 789:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay789, 2)) + 0.5*pow(-my789 + y789, 2)/pow(sigmay789, 2);
            break;
        case 790:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay790, 2)) + 0.5*pow(-my790 + y790, 2)/pow(sigmay790, 2);
            break;
        case 791:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay791, 2)) + 0.5*pow(-my791 + y791, 2)/pow(sigmay791, 2);
            break;
        case 792:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay792, 2)) + 0.5*pow(-my792 + y792, 2)/pow(sigmay792, 2);
            break;
        case 793:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay793, 2)) + 0.5*pow(-my793 + y793, 2)/pow(sigmay793, 2);
            break;
        case 794:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay794, 2)) + 0.5*pow(-my794 + y794, 2)/pow(sigmay794, 2);
            break;
        case 795:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay795, 2)) + 0.5*pow(-my795 + y795, 2)/pow(sigmay795, 2);
            break;
        case 796:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay796, 2)) + 0.5*pow(-my796 + y796, 2)/pow(sigmay796, 2);
            break;
        case 797:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay797, 2)) + 0.5*pow(-my797 + y797, 2)/pow(sigmay797, 2);
            break;
        case 798:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay798, 2)) + 0.5*pow(-my798 + y798, 2)/pow(sigmay798, 2);
            break;
        case 799:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay799, 2)) + 0.5*pow(-my799 + y799, 2)/pow(sigmay799, 2);
            break;
        case 800:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay800, 2)) + 0.5*pow(-my800 + y800, 2)/pow(sigmay800, 2);
            break;
        case 801:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay801, 2)) + 0.5*pow(-my801 + y801, 2)/pow(sigmay801, 2);
            break;
        case 802:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay802, 2)) + 0.5*pow(-my802 + y802, 2)/pow(sigmay802, 2);
            break;
        case 803:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay803, 2)) + 0.5*pow(-my803 + y803, 2)/pow(sigmay803, 2);
            break;
        case 804:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay804, 2)) + 0.5*pow(-my804 + y804, 2)/pow(sigmay804, 2);
            break;
        case 805:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay805, 2)) + 0.5*pow(-my805 + y805, 2)/pow(sigmay805, 2);
            break;
        case 806:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay806, 2)) + 0.5*pow(-my806 + y806, 2)/pow(sigmay806, 2);
            break;
        case 807:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay807, 2)) + 0.5*pow(-my807 + y807, 2)/pow(sigmay807, 2);
            break;
        case 808:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay808, 2)) + 0.5*pow(-my808 + y808, 2)/pow(sigmay808, 2);
            break;
        case 809:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay809, 2)) + 0.5*pow(-my809 + y809, 2)/pow(sigmay809, 2);
            break;
        case 810:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay810, 2)) + 0.5*pow(-my810 + y810, 2)/pow(sigmay810, 2);
            break;
        case 811:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay811, 2)) + 0.5*pow(-my811 + y811, 2)/pow(sigmay811, 2);
            break;
        case 812:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay812, 2)) + 0.5*pow(-my812 + y812, 2)/pow(sigmay812, 2);
            break;
        case 813:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay813, 2)) + 0.5*pow(-my813 + y813, 2)/pow(sigmay813, 2);
            break;
        case 814:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay814, 2)) + 0.5*pow(-my814 + y814, 2)/pow(sigmay814, 2);
            break;
        case 815:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay815, 2)) + 0.5*pow(-my815 + y815, 2)/pow(sigmay815, 2);
            break;
        case 816:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay816, 2)) + 0.5*pow(-my816 + y816, 2)/pow(sigmay816, 2);
            break;
        case 817:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay817, 2)) + 0.5*pow(-my817 + y817, 2)/pow(sigmay817, 2);
            break;
        case 818:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay818, 2)) + 0.5*pow(-my818 + y818, 2)/pow(sigmay818, 2);
            break;
        case 819:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay819, 2)) + 0.5*pow(-my819 + y819, 2)/pow(sigmay819, 2);
            break;
        case 820:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay820, 2)) + 0.5*pow(-my820 + y820, 2)/pow(sigmay820, 2);
            break;
        case 821:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay821, 2)) + 0.5*pow(-my821 + y821, 2)/pow(sigmay821, 2);
            break;
        case 822:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay822, 2)) + 0.5*pow(-my822 + y822, 2)/pow(sigmay822, 2);
            break;
        case 823:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay823, 2)) + 0.5*pow(-my823 + y823, 2)/pow(sigmay823, 2);
            break;
        case 824:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay824, 2)) + 0.5*pow(-my824 + y824, 2)/pow(sigmay824, 2);
            break;
        case 825:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay825, 2)) + 0.5*pow(-my825 + y825, 2)/pow(sigmay825, 2);
            break;
        case 826:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay826, 2)) + 0.5*pow(-my826 + y826, 2)/pow(sigmay826, 2);
            break;
        case 827:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay827, 2)) + 0.5*pow(-my827 + y827, 2)/pow(sigmay827, 2);
            break;
        case 828:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay828, 2)) + 0.5*pow(-my828 + y828, 2)/pow(sigmay828, 2);
            break;
        case 829:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay829, 2)) + 0.5*pow(-my829 + y829, 2)/pow(sigmay829, 2);
            break;
        case 830:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay830, 2)) + 0.5*pow(-my830 + y830, 2)/pow(sigmay830, 2);
            break;
        case 831:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay831, 2)) + 0.5*pow(-my831 + y831, 2)/pow(sigmay831, 2);
            break;
        case 832:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay832, 2)) + 0.5*pow(-my832 + y832, 2)/pow(sigmay832, 2);
            break;
        case 833:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay833, 2)) + 0.5*pow(-my833 + y833, 2)/pow(sigmay833, 2);
            break;
        case 834:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay834, 2)) + 0.5*pow(-my834 + y834, 2)/pow(sigmay834, 2);
            break;
        case 835:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay835, 2)) + 0.5*pow(-my835 + y835, 2)/pow(sigmay835, 2);
            break;
        case 836:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay836, 2)) + 0.5*pow(-my836 + y836, 2)/pow(sigmay836, 2);
            break;
        case 837:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay837, 2)) + 0.5*pow(-my837 + y837, 2)/pow(sigmay837, 2);
            break;
        case 838:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay838, 2)) + 0.5*pow(-my838 + y838, 2)/pow(sigmay838, 2);
            break;
        case 839:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay839, 2)) + 0.5*pow(-my839 + y839, 2)/pow(sigmay839, 2);
            break;
        case 840:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay840, 2)) + 0.5*pow(-my840 + y840, 2)/pow(sigmay840, 2);
            break;
        case 841:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay841, 2)) + 0.5*pow(-my841 + y841, 2)/pow(sigmay841, 2);
            break;
        case 842:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay842, 2)) + 0.5*pow(-my842 + y842, 2)/pow(sigmay842, 2);
            break;
        case 843:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay843, 2)) + 0.5*pow(-my843 + y843, 2)/pow(sigmay843, 2);
            break;
        case 844:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay844, 2)) + 0.5*pow(-my844 + y844, 2)/pow(sigmay844, 2);
            break;
        case 845:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay845, 2)) + 0.5*pow(-my845 + y845, 2)/pow(sigmay845, 2);
            break;
        case 846:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay846, 2)) + 0.5*pow(-my846 + y846, 2)/pow(sigmay846, 2);
            break;
        case 847:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay847, 2)) + 0.5*pow(-my847 + y847, 2)/pow(sigmay847, 2);
            break;
        case 848:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay848, 2)) + 0.5*pow(-my848 + y848, 2)/pow(sigmay848, 2);
            break;
        case 849:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay849, 2)) + 0.5*pow(-my849 + y849, 2)/pow(sigmay849, 2);
            break;
        case 850:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay850, 2)) + 0.5*pow(-my850 + y850, 2)/pow(sigmay850, 2);
            break;
        case 851:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay851, 2)) + 0.5*pow(-my851 + y851, 2)/pow(sigmay851, 2);
            break;
        case 852:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay852, 2)) + 0.5*pow(-my852 + y852, 2)/pow(sigmay852, 2);
            break;
        case 853:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay853, 2)) + 0.5*pow(-my853 + y853, 2)/pow(sigmay853, 2);
            break;
        case 854:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay854, 2)) + 0.5*pow(-my854 + y854, 2)/pow(sigmay854, 2);
            break;
        case 855:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay855, 2)) + 0.5*pow(-my855 + y855, 2)/pow(sigmay855, 2);
            break;
        case 856:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay856, 2)) + 0.5*pow(-my856 + y856, 2)/pow(sigmay856, 2);
            break;
        case 857:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay857, 2)) + 0.5*pow(-my857 + y857, 2)/pow(sigmay857, 2);
            break;
        case 858:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay858, 2)) + 0.5*pow(-my858 + y858, 2)/pow(sigmay858, 2);
            break;
        case 859:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay859, 2)) + 0.5*pow(-my859 + y859, 2)/pow(sigmay859, 2);
            break;
        case 860:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay860, 2)) + 0.5*pow(-my860 + y860, 2)/pow(sigmay860, 2);
            break;
        case 861:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay861, 2)) + 0.5*pow(-my861 + y861, 2)/pow(sigmay861, 2);
            break;
        case 862:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay862, 2)) + 0.5*pow(-my862 + y862, 2)/pow(sigmay862, 2);
            break;
        case 863:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay863, 2)) + 0.5*pow(-my863 + y863, 2)/pow(sigmay863, 2);
            break;
        case 864:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay864, 2)) + 0.5*pow(-my864 + y864, 2)/pow(sigmay864, 2);
            break;
        case 865:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay865, 2)) + 0.5*pow(-my865 + y865, 2)/pow(sigmay865, 2);
            break;
        case 866:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay866, 2)) + 0.5*pow(-my866 + y866, 2)/pow(sigmay866, 2);
            break;
        case 867:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay867, 2)) + 0.5*pow(-my867 + y867, 2)/pow(sigmay867, 2);
            break;
        case 868:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay868, 2)) + 0.5*pow(-my868 + y868, 2)/pow(sigmay868, 2);
            break;
        case 869:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay869, 2)) + 0.5*pow(-my869 + y869, 2)/pow(sigmay869, 2);
            break;
        case 870:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay870, 2)) + 0.5*pow(-my870 + y870, 2)/pow(sigmay870, 2);
            break;
        case 871:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay871, 2)) + 0.5*pow(-my871 + y871, 2)/pow(sigmay871, 2);
            break;
        case 872:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay872, 2)) + 0.5*pow(-my872 + y872, 2)/pow(sigmay872, 2);
            break;
        case 873:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay873, 2)) + 0.5*pow(-my873 + y873, 2)/pow(sigmay873, 2);
            break;
        case 874:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay874, 2)) + 0.5*pow(-my874 + y874, 2)/pow(sigmay874, 2);
            break;
        case 875:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay875, 2)) + 0.5*pow(-my875 + y875, 2)/pow(sigmay875, 2);
            break;
        case 876:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay876, 2)) + 0.5*pow(-my876 + y876, 2)/pow(sigmay876, 2);
            break;
        case 877:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay877, 2)) + 0.5*pow(-my877 + y877, 2)/pow(sigmay877, 2);
            break;
        case 878:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay878, 2)) + 0.5*pow(-my878 + y878, 2)/pow(sigmay878, 2);
            break;
        case 879:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay879, 2)) + 0.5*pow(-my879 + y879, 2)/pow(sigmay879, 2);
            break;
        case 880:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay880, 2)) + 0.5*pow(-my880 + y880, 2)/pow(sigmay880, 2);
            break;
        case 881:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay881, 2)) + 0.5*pow(-my881 + y881, 2)/pow(sigmay881, 2);
            break;
        case 882:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay882, 2)) + 0.5*pow(-my882 + y882, 2)/pow(sigmay882, 2);
            break;
        case 883:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay883, 2)) + 0.5*pow(-my883 + y883, 2)/pow(sigmay883, 2);
            break;
        case 884:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay884, 2)) + 0.5*pow(-my884 + y884, 2)/pow(sigmay884, 2);
            break;
        case 885:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay885, 2)) + 0.5*pow(-my885 + y885, 2)/pow(sigmay885, 2);
            break;
        case 886:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay886, 2)) + 0.5*pow(-my886 + y886, 2)/pow(sigmay886, 2);
            break;
        case 887:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay887, 2)) + 0.5*pow(-my887 + y887, 2)/pow(sigmay887, 2);
            break;
        case 888:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay888, 2)) + 0.5*pow(-my888 + y888, 2)/pow(sigmay888, 2);
            break;
        case 889:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay889, 2)) + 0.5*pow(-my889 + y889, 2)/pow(sigmay889, 2);
            break;
        case 890:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay890, 2)) + 0.5*pow(-my890 + y890, 2)/pow(sigmay890, 2);
            break;
        case 891:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay891, 2)) + 0.5*pow(-my891 + y891, 2)/pow(sigmay891, 2);
            break;
        case 892:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay892, 2)) + 0.5*pow(-my892 + y892, 2)/pow(sigmay892, 2);
            break;
        case 893:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay893, 2)) + 0.5*pow(-my893 + y893, 2)/pow(sigmay893, 2);
            break;
        case 894:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay894, 2)) + 0.5*pow(-my894 + y894, 2)/pow(sigmay894, 2);
            break;
        case 895:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay895, 2)) + 0.5*pow(-my895 + y895, 2)/pow(sigmay895, 2);
            break;
        case 896:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay896, 2)) + 0.5*pow(-my896 + y896, 2)/pow(sigmay896, 2);
            break;
        case 897:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay897, 2)) + 0.5*pow(-my897 + y897, 2)/pow(sigmay897, 2);
            break;
        case 898:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay898, 2)) + 0.5*pow(-my898 + y898, 2)/pow(sigmay898, 2);
            break;
        case 899:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay899, 2)) + 0.5*pow(-my899 + y899, 2)/pow(sigmay899, 2);
            break;
        case 900:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay900, 2)) + 0.5*pow(-my900 + y900, 2)/pow(sigmay900, 2);
            break;
        case 901:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay901, 2)) + 0.5*pow(-my901 + y901, 2)/pow(sigmay901, 2);
            break;
        case 902:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay902, 2)) + 0.5*pow(-my902 + y902, 2)/pow(sigmay902, 2);
            break;
        case 903:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay903, 2)) + 0.5*pow(-my903 + y903, 2)/pow(sigmay903, 2);
            break;
        case 904:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay904, 2)) + 0.5*pow(-my904 + y904, 2)/pow(sigmay904, 2);
            break;
        case 905:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay905, 2)) + 0.5*pow(-my905 + y905, 2)/pow(sigmay905, 2);
            break;
        case 906:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay906, 2)) + 0.5*pow(-my906 + y906, 2)/pow(sigmay906, 2);
            break;
        case 907:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay907, 2)) + 0.5*pow(-my907 + y907, 2)/pow(sigmay907, 2);
            break;
        case 908:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay908, 2)) + 0.5*pow(-my908 + y908, 2)/pow(sigmay908, 2);
            break;
        case 909:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay909, 2)) + 0.5*pow(-my909 + y909, 2)/pow(sigmay909, 2);
            break;
        case 910:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay910, 2)) + 0.5*pow(-my910 + y910, 2)/pow(sigmay910, 2);
            break;
        case 911:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay911, 2)) + 0.5*pow(-my911 + y911, 2)/pow(sigmay911, 2);
            break;
        case 912:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay912, 2)) + 0.5*pow(-my912 + y912, 2)/pow(sigmay912, 2);
            break;
        case 913:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay913, 2)) + 0.5*pow(-my913 + y913, 2)/pow(sigmay913, 2);
            break;
        case 914:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay914, 2)) + 0.5*pow(-my914 + y914, 2)/pow(sigmay914, 2);
            break;
        case 915:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay915, 2)) + 0.5*pow(-my915 + y915, 2)/pow(sigmay915, 2);
            break;
        case 916:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay916, 2)) + 0.5*pow(-my916 + y916, 2)/pow(sigmay916, 2);
            break;
        case 917:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay917, 2)) + 0.5*pow(-my917 + y917, 2)/pow(sigmay917, 2);
            break;
        case 918:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay918, 2)) + 0.5*pow(-my918 + y918, 2)/pow(sigmay918, 2);
            break;
        case 919:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay919, 2)) + 0.5*pow(-my919 + y919, 2)/pow(sigmay919, 2);
            break;
        case 920:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay920, 2)) + 0.5*pow(-my920 + y920, 2)/pow(sigmay920, 2);
            break;
        case 921:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay921, 2)) + 0.5*pow(-my921 + y921, 2)/pow(sigmay921, 2);
            break;
        case 922:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay922, 2)) + 0.5*pow(-my922 + y922, 2)/pow(sigmay922, 2);
            break;
        case 923:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay923, 2)) + 0.5*pow(-my923 + y923, 2)/pow(sigmay923, 2);
            break;
        case 924:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay924, 2)) + 0.5*pow(-my924 + y924, 2)/pow(sigmay924, 2);
            break;
        case 925:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay925, 2)) + 0.5*pow(-my925 + y925, 2)/pow(sigmay925, 2);
            break;
        case 926:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay926, 2)) + 0.5*pow(-my926 + y926, 2)/pow(sigmay926, 2);
            break;
        case 927:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay927, 2)) + 0.5*pow(-my927 + y927, 2)/pow(sigmay927, 2);
            break;
        case 928:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay928, 2)) + 0.5*pow(-my928 + y928, 2)/pow(sigmay928, 2);
            break;
        case 929:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay929, 2)) + 0.5*pow(-my929 + y929, 2)/pow(sigmay929, 2);
            break;
        case 930:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay930, 2)) + 0.5*pow(-my930 + y930, 2)/pow(sigmay930, 2);
            break;
        case 931:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay931, 2)) + 0.5*pow(-my931 + y931, 2)/pow(sigmay931, 2);
            break;
        case 932:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay932, 2)) + 0.5*pow(-my932 + y932, 2)/pow(sigmay932, 2);
            break;
        case 933:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay933, 2)) + 0.5*pow(-my933 + y933, 2)/pow(sigmay933, 2);
            break;
        case 934:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay934, 2)) + 0.5*pow(-my934 + y934, 2)/pow(sigmay934, 2);
            break;
        case 935:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay935, 2)) + 0.5*pow(-my935 + y935, 2)/pow(sigmay935, 2);
            break;
        case 936:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay936, 2)) + 0.5*pow(-my936 + y936, 2)/pow(sigmay936, 2);
            break;
        case 937:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay937, 2)) + 0.5*pow(-my937 + y937, 2)/pow(sigmay937, 2);
            break;
        case 938:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay938, 2)) + 0.5*pow(-my938 + y938, 2)/pow(sigmay938, 2);
            break;
        case 939:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay939, 2)) + 0.5*pow(-my939 + y939, 2)/pow(sigmay939, 2);
            break;
        case 940:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay940, 2)) + 0.5*pow(-my940 + y940, 2)/pow(sigmay940, 2);
            break;
        case 941:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay941, 2)) + 0.5*pow(-my941 + y941, 2)/pow(sigmay941, 2);
            break;
        case 942:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay942, 2)) + 0.5*pow(-my942 + y942, 2)/pow(sigmay942, 2);
            break;
        case 943:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay943, 2)) + 0.5*pow(-my943 + y943, 2)/pow(sigmay943, 2);
            break;
        case 944:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay944, 2)) + 0.5*pow(-my944 + y944, 2)/pow(sigmay944, 2);
            break;
        case 945:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay945, 2)) + 0.5*pow(-my945 + y945, 2)/pow(sigmay945, 2);
            break;
        case 946:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay946, 2)) + 0.5*pow(-my946 + y946, 2)/pow(sigmay946, 2);
            break;
        case 947:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay947, 2)) + 0.5*pow(-my947 + y947, 2)/pow(sigmay947, 2);
            break;
        case 948:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay948, 2)) + 0.5*pow(-my948 + y948, 2)/pow(sigmay948, 2);
            break;
        case 949:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay949, 2)) + 0.5*pow(-my949 + y949, 2)/pow(sigmay949, 2);
            break;
        case 950:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay950, 2)) + 0.5*pow(-my950 + y950, 2)/pow(sigmay950, 2);
            break;
        case 951:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay951, 2)) + 0.5*pow(-my951 + y951, 2)/pow(sigmay951, 2);
            break;
        case 952:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay952, 2)) + 0.5*pow(-my952 + y952, 2)/pow(sigmay952, 2);
            break;
        case 953:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay953, 2)) + 0.5*pow(-my953 + y953, 2)/pow(sigmay953, 2);
            break;
        case 954:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay954, 2)) + 0.5*pow(-my954 + y954, 2)/pow(sigmay954, 2);
            break;
        case 955:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay955, 2)) + 0.5*pow(-my955 + y955, 2)/pow(sigmay955, 2);
            break;
        case 956:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay956, 2)) + 0.5*pow(-my956 + y956, 2)/pow(sigmay956, 2);
            break;
        case 957:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay957, 2)) + 0.5*pow(-my957 + y957, 2)/pow(sigmay957, 2);
            break;
        case 958:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay958, 2)) + 0.5*pow(-my958 + y958, 2)/pow(sigmay958, 2);
            break;
        case 959:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay959, 2)) + 0.5*pow(-my959 + y959, 2)/pow(sigmay959, 2);
            break;
        case 960:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay960, 2)) + 0.5*pow(-my960 + y960, 2)/pow(sigmay960, 2);
            break;
        case 961:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay961, 2)) + 0.5*pow(-my961 + y961, 2)/pow(sigmay961, 2);
            break;
        case 962:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay962, 2)) + 0.5*pow(-my962 + y962, 2)/pow(sigmay962, 2);
            break;
        case 963:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay963, 2)) + 0.5*pow(-my963 + y963, 2)/pow(sigmay963, 2);
            break;
        case 964:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay964, 2)) + 0.5*pow(-my964 + y964, 2)/pow(sigmay964, 2);
            break;
        case 965:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay965, 2)) + 0.5*pow(-my965 + y965, 2)/pow(sigmay965, 2);
            break;
        case 966:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay966, 2)) + 0.5*pow(-my966 + y966, 2)/pow(sigmay966, 2);
            break;
        case 967:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay967, 2)) + 0.5*pow(-my967 + y967, 2)/pow(sigmay967, 2);
            break;
        case 968:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay968, 2)) + 0.5*pow(-my968 + y968, 2)/pow(sigmay968, 2);
            break;
        case 969:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay969, 2)) + 0.5*pow(-my969 + y969, 2)/pow(sigmay969, 2);
            break;
        case 970:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay970, 2)) + 0.5*pow(-my970 + y970, 2)/pow(sigmay970, 2);
            break;
        case 971:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay971, 2)) + 0.5*pow(-my971 + y971, 2)/pow(sigmay971, 2);
            break;
        case 972:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay972, 2)) + 0.5*pow(-my972 + y972, 2)/pow(sigmay972, 2);
            break;
        case 973:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay973, 2)) + 0.5*pow(-my973 + y973, 2)/pow(sigmay973, 2);
            break;
        case 974:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay974, 2)) + 0.5*pow(-my974 + y974, 2)/pow(sigmay974, 2);
            break;
        case 975:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay975, 2)) + 0.5*pow(-my975 + y975, 2)/pow(sigmay975, 2);
            break;
        case 976:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay976, 2)) + 0.5*pow(-my976 + y976, 2)/pow(sigmay976, 2);
            break;
        case 977:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay977, 2)) + 0.5*pow(-my977 + y977, 2)/pow(sigmay977, 2);
            break;
        case 978:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay978, 2)) + 0.5*pow(-my978 + y978, 2)/pow(sigmay978, 2);
            break;
        case 979:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay979, 2)) + 0.5*pow(-my979 + y979, 2)/pow(sigmay979, 2);
            break;
        case 980:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay980, 2)) + 0.5*pow(-my980 + y980, 2)/pow(sigmay980, 2);
            break;
        case 981:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay981, 2)) + 0.5*pow(-my981 + y981, 2)/pow(sigmay981, 2);
            break;
        case 982:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay982, 2)) + 0.5*pow(-my982 + y982, 2)/pow(sigmay982, 2);
            break;
        case 983:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay983, 2)) + 0.5*pow(-my983 + y983, 2)/pow(sigmay983, 2);
            break;
        case 984:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay984, 2)) + 0.5*pow(-my984 + y984, 2)/pow(sigmay984, 2);
            break;
        case 985:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay985, 2)) + 0.5*pow(-my985 + y985, 2)/pow(sigmay985, 2);
            break;
        case 986:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay986, 2)) + 0.5*pow(-my986 + y986, 2)/pow(sigmay986, 2);
            break;
        case 987:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay987, 2)) + 0.5*pow(-my987 + y987, 2)/pow(sigmay987, 2);
            break;
        case 988:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay988, 2)) + 0.5*pow(-my988 + y988, 2)/pow(sigmay988, 2);
            break;
        case 989:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay989, 2)) + 0.5*pow(-my989 + y989, 2)/pow(sigmay989, 2);
            break;
        case 990:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay990, 2)) + 0.5*pow(-my990 + y990, 2)/pow(sigmay990, 2);
            break;
        case 991:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay991, 2)) + 0.5*pow(-my991 + y991, 2)/pow(sigmay991, 2);
            break;
        case 992:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay992, 2)) + 0.5*pow(-my992 + y992, 2)/pow(sigmay992, 2);
            break;
        case 993:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay993, 2)) + 0.5*pow(-my993 + y993, 2)/pow(sigmay993, 2);
            break;
        case 994:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay994, 2)) + 0.5*pow(-my994 + y994, 2)/pow(sigmay994, 2);
            break;
        case 995:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay995, 2)) + 0.5*pow(-my995 + y995, 2)/pow(sigmay995, 2);
            break;
        case 996:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay996, 2)) + 0.5*pow(-my996 + y996, 2)/pow(sigmay996, 2);
            break;
        case 997:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay997, 2)) + 0.5*pow(-my997 + y997, 2)/pow(sigmay997, 2);
            break;
        case 998:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay998, 2)) + 0.5*pow(-my998 + y998, 2)/pow(sigmay998, 2);
            break;
        case 999:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay999, 2)) + 0.5*pow(-my999 + y999, 2)/pow(sigmay999, 2);
            break;
        case 1000:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1000, 2)) + 0.5*pow(-my1000 + y1000, 2)/pow(sigmay1000, 2);
            break;
        case 1001:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1001, 2)) + 0.5*pow(-my1001 + y1001, 2)/pow(sigmay1001, 2);
            break;
        case 1002:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1002, 2)) + 0.5*pow(-my1002 + y1002, 2)/pow(sigmay1002, 2);
            break;
        case 1003:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1003, 2)) + 0.5*pow(-my1003 + y1003, 2)/pow(sigmay1003, 2);
            break;
        case 1004:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1004, 2)) + 0.5*pow(-my1004 + y1004, 2)/pow(sigmay1004, 2);
            break;
        case 1005:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1005, 2)) + 0.5*pow(-my1005 + y1005, 2)/pow(sigmay1005, 2);
            break;
        case 1006:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1006, 2)) + 0.5*pow(-my1006 + y1006, 2)/pow(sigmay1006, 2);
            break;
        case 1007:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1007, 2)) + 0.5*pow(-my1007 + y1007, 2)/pow(sigmay1007, 2);
            break;
        case 1008:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1008, 2)) + 0.5*pow(-my1008 + y1008, 2)/pow(sigmay1008, 2);
            break;
        case 1009:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1009, 2)) + 0.5*pow(-my1009 + y1009, 2)/pow(sigmay1009, 2);
            break;
        case 1010:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1010, 2)) + 0.5*pow(-my1010 + y1010, 2)/pow(sigmay1010, 2);
            break;
        case 1011:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1011, 2)) + 0.5*pow(-my1011 + y1011, 2)/pow(sigmay1011, 2);
            break;
        case 1012:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1012, 2)) + 0.5*pow(-my1012 + y1012, 2)/pow(sigmay1012, 2);
            break;
        case 1013:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1013, 2)) + 0.5*pow(-my1013 + y1013, 2)/pow(sigmay1013, 2);
            break;
        case 1014:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1014, 2)) + 0.5*pow(-my1014 + y1014, 2)/pow(sigmay1014, 2);
            break;
        case 1015:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1015, 2)) + 0.5*pow(-my1015 + y1015, 2)/pow(sigmay1015, 2);
            break;
        case 1016:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1016, 2)) + 0.5*pow(-my1016 + y1016, 2)/pow(sigmay1016, 2);
            break;
        case 1017:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1017, 2)) + 0.5*pow(-my1017 + y1017, 2)/pow(sigmay1017, 2);
            break;
        case 1018:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1018, 2)) + 0.5*pow(-my1018 + y1018, 2)/pow(sigmay1018, 2);
            break;
        case 1019:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1019, 2)) + 0.5*pow(-my1019 + y1019, 2)/pow(sigmay1019, 2);
            break;
        case 1020:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1020, 2)) + 0.5*pow(-my1020 + y1020, 2)/pow(sigmay1020, 2);
            break;
        case 1021:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1021, 2)) + 0.5*pow(-my1021 + y1021, 2)/pow(sigmay1021, 2);
            break;
        case 1022:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1022, 2)) + 0.5*pow(-my1022 + y1022, 2)/pow(sigmay1022, 2);
            break;
        case 1023:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1023, 2)) + 0.5*pow(-my1023 + y1023, 2)/pow(sigmay1023, 2);
            break;
        case 1024:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1024, 2)) + 0.5*pow(-my1024 + y1024, 2)/pow(sigmay1024, 2);
            break;
        case 1025:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1025, 2)) + 0.5*pow(-my1025 + y1025, 2)/pow(sigmay1025, 2);
            break;
        case 1026:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1026, 2)) + 0.5*pow(-my1026 + y1026, 2)/pow(sigmay1026, 2);
            break;
        case 1027:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1027, 2)) + 0.5*pow(-my1027 + y1027, 2)/pow(sigmay1027, 2);
            break;
        case 1028:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1028, 2)) + 0.5*pow(-my1028 + y1028, 2)/pow(sigmay1028, 2);
            break;
        case 1029:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1029, 2)) + 0.5*pow(-my1029 + y1029, 2)/pow(sigmay1029, 2);
            break;
        case 1030:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1030, 2)) + 0.5*pow(-my1030 + y1030, 2)/pow(sigmay1030, 2);
            break;
        case 1031:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1031, 2)) + 0.5*pow(-my1031 + y1031, 2)/pow(sigmay1031, 2);
            break;
        case 1032:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1032, 2)) + 0.5*pow(-my1032 + y1032, 2)/pow(sigmay1032, 2);
            break;
        case 1033:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1033, 2)) + 0.5*pow(-my1033 + y1033, 2)/pow(sigmay1033, 2);
            break;
        case 1034:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1034, 2)) + 0.5*pow(-my1034 + y1034, 2)/pow(sigmay1034, 2);
            break;
        case 1035:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1035, 2)) + 0.5*pow(-my1035 + y1035, 2)/pow(sigmay1035, 2);
            break;
        case 1036:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1036, 2)) + 0.5*pow(-my1036 + y1036, 2)/pow(sigmay1036, 2);
            break;
        case 1037:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1037, 2)) + 0.5*pow(-my1037 + y1037, 2)/pow(sigmay1037, 2);
            break;
        case 1038:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1038, 2)) + 0.5*pow(-my1038 + y1038, 2)/pow(sigmay1038, 2);
            break;
        case 1039:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1039, 2)) + 0.5*pow(-my1039 + y1039, 2)/pow(sigmay1039, 2);
            break;
        case 1040:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1040, 2)) + 0.5*pow(-my1040 + y1040, 2)/pow(sigmay1040, 2);
            break;
        case 1041:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1041, 2)) + 0.5*pow(-my1041 + y1041, 2)/pow(sigmay1041, 2);
            break;
        case 1042:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1042, 2)) + 0.5*pow(-my1042 + y1042, 2)/pow(sigmay1042, 2);
            break;
        case 1043:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1043, 2)) + 0.5*pow(-my1043 + y1043, 2)/pow(sigmay1043, 2);
            break;
        case 1044:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1044, 2)) + 0.5*pow(-my1044 + y1044, 2)/pow(sigmay1044, 2);
            break;
        case 1045:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1045, 2)) + 0.5*pow(-my1045 + y1045, 2)/pow(sigmay1045, 2);
            break;
        case 1046:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1046, 2)) + 0.5*pow(-my1046 + y1046, 2)/pow(sigmay1046, 2);
            break;
        case 1047:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1047, 2)) + 0.5*pow(-my1047 + y1047, 2)/pow(sigmay1047, 2);
            break;
        case 1048:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1048, 2)) + 0.5*pow(-my1048 + y1048, 2)/pow(sigmay1048, 2);
            break;
        case 1049:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1049, 2)) + 0.5*pow(-my1049 + y1049, 2)/pow(sigmay1049, 2);
            break;
        case 1050:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1050, 2)) + 0.5*pow(-my1050 + y1050, 2)/pow(sigmay1050, 2);
            break;
        case 1051:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1051, 2)) + 0.5*pow(-my1051 + y1051, 2)/pow(sigmay1051, 2);
            break;
        case 1052:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1052, 2)) + 0.5*pow(-my1052 + y1052, 2)/pow(sigmay1052, 2);
            break;
        case 1053:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1053, 2)) + 0.5*pow(-my1053 + y1053, 2)/pow(sigmay1053, 2);
            break;
        case 1054:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1054, 2)) + 0.5*pow(-my1054 + y1054, 2)/pow(sigmay1054, 2);
            break;
        case 1055:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1055, 2)) + 0.5*pow(-my1055 + y1055, 2)/pow(sigmay1055, 2);
            break;
        case 1056:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1056, 2)) + 0.5*pow(-my1056 + y1056, 2)/pow(sigmay1056, 2);
            break;
        case 1057:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1057, 2)) + 0.5*pow(-my1057 + y1057, 2)/pow(sigmay1057, 2);
            break;
        case 1058:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1058, 2)) + 0.5*pow(-my1058 + y1058, 2)/pow(sigmay1058, 2);
            break;
        case 1059:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1059, 2)) + 0.5*pow(-my1059 + y1059, 2)/pow(sigmay1059, 2);
            break;
        case 1060:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1060, 2)) + 0.5*pow(-my1060 + y1060, 2)/pow(sigmay1060, 2);
            break;
        case 1061:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1061, 2)) + 0.5*pow(-my1061 + y1061, 2)/pow(sigmay1061, 2);
            break;
        case 1062:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1062, 2)) + 0.5*pow(-my1062 + y1062, 2)/pow(sigmay1062, 2);
            break;
        case 1063:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1063, 2)) + 0.5*pow(-my1063 + y1063, 2)/pow(sigmay1063, 2);
            break;
        case 1064:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1064, 2)) + 0.5*pow(-my1064 + y1064, 2)/pow(sigmay1064, 2);
            break;
        case 1065:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1065, 2)) + 0.5*pow(-my1065 + y1065, 2)/pow(sigmay1065, 2);
            break;
        case 1066:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1066, 2)) + 0.5*pow(-my1066 + y1066, 2)/pow(sigmay1066, 2);
            break;
        case 1067:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1067, 2)) + 0.5*pow(-my1067 + y1067, 2)/pow(sigmay1067, 2);
            break;
        case 1068:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1068, 2)) + 0.5*pow(-my1068 + y1068, 2)/pow(sigmay1068, 2);
            break;
        case 1069:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1069, 2)) + 0.5*pow(-my1069 + y1069, 2)/pow(sigmay1069, 2);
            break;
        case 1070:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1070, 2)) + 0.5*pow(-my1070 + y1070, 2)/pow(sigmay1070, 2);
            break;
        case 1071:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1071, 2)) + 0.5*pow(-my1071 + y1071, 2)/pow(sigmay1071, 2);
            break;
        case 1072:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1072, 2)) + 0.5*pow(-my1072 + y1072, 2)/pow(sigmay1072, 2);
            break;
        case 1073:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1073, 2)) + 0.5*pow(-my1073 + y1073, 2)/pow(sigmay1073, 2);
            break;
        case 1074:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1074, 2)) + 0.5*pow(-my1074 + y1074, 2)/pow(sigmay1074, 2);
            break;
        case 1075:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1075, 2)) + 0.5*pow(-my1075 + y1075, 2)/pow(sigmay1075, 2);
            break;
        case 1076:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1076, 2)) + 0.5*pow(-my1076 + y1076, 2)/pow(sigmay1076, 2);
            break;
        case 1077:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1077, 2)) + 0.5*pow(-my1077 + y1077, 2)/pow(sigmay1077, 2);
            break;
        case 1078:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1078, 2)) + 0.5*pow(-my1078 + y1078, 2)/pow(sigmay1078, 2);
            break;
        case 1079:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1079, 2)) + 0.5*pow(-my1079 + y1079, 2)/pow(sigmay1079, 2);
            break;
        case 1080:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1080, 2)) + 0.5*pow(-my1080 + y1080, 2)/pow(sigmay1080, 2);
            break;
        case 1081:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1081, 2)) + 0.5*pow(-my1081 + y1081, 2)/pow(sigmay1081, 2);
            break;
        case 1082:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1082, 2)) + 0.5*pow(-my1082 + y1082, 2)/pow(sigmay1082, 2);
            break;
        case 1083:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1083, 2)) + 0.5*pow(-my1083 + y1083, 2)/pow(sigmay1083, 2);
            break;
        case 1084:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1084, 2)) + 0.5*pow(-my1084 + y1084, 2)/pow(sigmay1084, 2);
            break;
        case 1085:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1085, 2)) + 0.5*pow(-my1085 + y1085, 2)/pow(sigmay1085, 2);
            break;
        case 1086:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1086, 2)) + 0.5*pow(-my1086 + y1086, 2)/pow(sigmay1086, 2);
            break;
        case 1087:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1087, 2)) + 0.5*pow(-my1087 + y1087, 2)/pow(sigmay1087, 2);
            break;
        case 1088:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1088, 2)) + 0.5*pow(-my1088 + y1088, 2)/pow(sigmay1088, 2);
            break;
        case 1089:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1089, 2)) + 0.5*pow(-my1089 + y1089, 2)/pow(sigmay1089, 2);
            break;
        case 1090:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1090, 2)) + 0.5*pow(-my1090 + y1090, 2)/pow(sigmay1090, 2);
            break;
        case 1091:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1091, 2)) + 0.5*pow(-my1091 + y1091, 2)/pow(sigmay1091, 2);
            break;
        case 1092:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1092, 2)) + 0.5*pow(-my1092 + y1092, 2)/pow(sigmay1092, 2);
            break;
        case 1093:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1093, 2)) + 0.5*pow(-my1093 + y1093, 2)/pow(sigmay1093, 2);
            break;
        case 1094:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1094, 2)) + 0.5*pow(-my1094 + y1094, 2)/pow(sigmay1094, 2);
            break;
        case 1095:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1095, 2)) + 0.5*pow(-my1095 + y1095, 2)/pow(sigmay1095, 2);
            break;
        case 1096:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1096, 2)) + 0.5*pow(-my1096 + y1096, 2)/pow(sigmay1096, 2);
            break;
        case 1097:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1097, 2)) + 0.5*pow(-my1097 + y1097, 2)/pow(sigmay1097, 2);
            break;
        case 1098:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1098, 2)) + 0.5*pow(-my1098 + y1098, 2)/pow(sigmay1098, 2);
            break;
        case 1099:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1099, 2)) + 0.5*pow(-my1099 + y1099, 2)/pow(sigmay1099, 2);
            break;
        case 1100:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1100, 2)) + 0.5*pow(-my1100 + y1100, 2)/pow(sigmay1100, 2);
            break;
        case 1101:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1101, 2)) + 0.5*pow(-my1101 + y1101, 2)/pow(sigmay1101, 2);
            break;
        case 1102:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1102, 2)) + 0.5*pow(-my1102 + y1102, 2)/pow(sigmay1102, 2);
            break;
        case 1103:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1103, 2)) + 0.5*pow(-my1103 + y1103, 2)/pow(sigmay1103, 2);
            break;
        case 1104:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1104, 2)) + 0.5*pow(-my1104 + y1104, 2)/pow(sigmay1104, 2);
            break;
        case 1105:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1105, 2)) + 0.5*pow(-my1105 + y1105, 2)/pow(sigmay1105, 2);
            break;
        case 1106:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1106, 2)) + 0.5*pow(-my1106 + y1106, 2)/pow(sigmay1106, 2);
            break;
        case 1107:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1107, 2)) + 0.5*pow(-my1107 + y1107, 2)/pow(sigmay1107, 2);
            break;
        case 1108:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1108, 2)) + 0.5*pow(-my1108 + y1108, 2)/pow(sigmay1108, 2);
            break;
        case 1109:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1109, 2)) + 0.5*pow(-my1109 + y1109, 2)/pow(sigmay1109, 2);
            break;
        case 1110:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1110, 2)) + 0.5*pow(-my1110 + y1110, 2)/pow(sigmay1110, 2);
            break;
        case 1111:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1111, 2)) + 0.5*pow(-my1111 + y1111, 2)/pow(sigmay1111, 2);
            break;
        case 1112:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1112, 2)) + 0.5*pow(-my1112 + y1112, 2)/pow(sigmay1112, 2);
            break;
        case 1113:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1113, 2)) + 0.5*pow(-my1113 + y1113, 2)/pow(sigmay1113, 2);
            break;
        case 1114:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1114, 2)) + 0.5*pow(-my1114 + y1114, 2)/pow(sigmay1114, 2);
            break;
        case 1115:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1115, 2)) + 0.5*pow(-my1115 + y1115, 2)/pow(sigmay1115, 2);
            break;
        case 1116:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1116, 2)) + 0.5*pow(-my1116 + y1116, 2)/pow(sigmay1116, 2);
            break;
        case 1117:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1117, 2)) + 0.5*pow(-my1117 + y1117, 2)/pow(sigmay1117, 2);
            break;
        case 1118:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1118, 2)) + 0.5*pow(-my1118 + y1118, 2)/pow(sigmay1118, 2);
            break;
        case 1119:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1119, 2)) + 0.5*pow(-my1119 + y1119, 2)/pow(sigmay1119, 2);
            break;
        case 1120:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1120, 2)) + 0.5*pow(-my1120 + y1120, 2)/pow(sigmay1120, 2);
            break;
        case 1121:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1121, 2)) + 0.5*pow(-my1121 + y1121, 2)/pow(sigmay1121, 2);
            break;
        case 1122:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1122, 2)) + 0.5*pow(-my1122 + y1122, 2)/pow(sigmay1122, 2);
            break;
        case 1123:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1123, 2)) + 0.5*pow(-my1123 + y1123, 2)/pow(sigmay1123, 2);
            break;
        case 1124:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1124, 2)) + 0.5*pow(-my1124 + y1124, 2)/pow(sigmay1124, 2);
            break;
        case 1125:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1125, 2)) + 0.5*pow(-my1125 + y1125, 2)/pow(sigmay1125, 2);
            break;
        case 1126:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1126, 2)) + 0.5*pow(-my1126 + y1126, 2)/pow(sigmay1126, 2);
            break;
        case 1127:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1127, 2)) + 0.5*pow(-my1127 + y1127, 2)/pow(sigmay1127, 2);
            break;
        case 1128:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1128, 2)) + 0.5*pow(-my1128 + y1128, 2)/pow(sigmay1128, 2);
            break;
        case 1129:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1129, 2)) + 0.5*pow(-my1129 + y1129, 2)/pow(sigmay1129, 2);
            break;
        case 1130:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1130, 2)) + 0.5*pow(-my1130 + y1130, 2)/pow(sigmay1130, 2);
            break;
        case 1131:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1131, 2)) + 0.5*pow(-my1131 + y1131, 2)/pow(sigmay1131, 2);
            break;
        case 1132:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1132, 2)) + 0.5*pow(-my1132 + y1132, 2)/pow(sigmay1132, 2);
            break;
        case 1133:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1133, 2)) + 0.5*pow(-my1133 + y1133, 2)/pow(sigmay1133, 2);
            break;
        case 1134:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1134, 2)) + 0.5*pow(-my1134 + y1134, 2)/pow(sigmay1134, 2);
            break;
        case 1135:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1135, 2)) + 0.5*pow(-my1135 + y1135, 2)/pow(sigmay1135, 2);
            break;
        case 1136:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1136, 2)) + 0.5*pow(-my1136 + y1136, 2)/pow(sigmay1136, 2);
            break;
        case 1137:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1137, 2)) + 0.5*pow(-my1137 + y1137, 2)/pow(sigmay1137, 2);
            break;
        case 1138:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1138, 2)) + 0.5*pow(-my1138 + y1138, 2)/pow(sigmay1138, 2);
            break;
        case 1139:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1139, 2)) + 0.5*pow(-my1139 + y1139, 2)/pow(sigmay1139, 2);
            break;
        case 1140:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1140, 2)) + 0.5*pow(-my1140 + y1140, 2)/pow(sigmay1140, 2);
            break;
        case 1141:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1141, 2)) + 0.5*pow(-my1141 + y1141, 2)/pow(sigmay1141, 2);
            break;
        case 1142:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1142, 2)) + 0.5*pow(-my1142 + y1142, 2)/pow(sigmay1142, 2);
            break;
        case 1143:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1143, 2)) + 0.5*pow(-my1143 + y1143, 2)/pow(sigmay1143, 2);
            break;
        case 1144:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1144, 2)) + 0.5*pow(-my1144 + y1144, 2)/pow(sigmay1144, 2);
            break;
        case 1145:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1145, 2)) + 0.5*pow(-my1145 + y1145, 2)/pow(sigmay1145, 2);
            break;
        case 1146:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1146, 2)) + 0.5*pow(-my1146 + y1146, 2)/pow(sigmay1146, 2);
            break;
        case 1147:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1147, 2)) + 0.5*pow(-my1147 + y1147, 2)/pow(sigmay1147, 2);
            break;
        case 1148:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1148, 2)) + 0.5*pow(-my1148 + y1148, 2)/pow(sigmay1148, 2);
            break;
        case 1149:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1149, 2)) + 0.5*pow(-my1149 + y1149, 2)/pow(sigmay1149, 2);
            break;
        case 1150:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1150, 2)) + 0.5*pow(-my1150 + y1150, 2)/pow(sigmay1150, 2);
            break;
        case 1151:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1151, 2)) + 0.5*pow(-my1151 + y1151, 2)/pow(sigmay1151, 2);
            break;
        case 1152:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1152, 2)) + 0.5*pow(-my1152 + y1152, 2)/pow(sigmay1152, 2);
            break;
        case 1153:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1153, 2)) + 0.5*pow(-my1153 + y1153, 2)/pow(sigmay1153, 2);
            break;
        case 1154:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1154, 2)) + 0.5*pow(-my1154 + y1154, 2)/pow(sigmay1154, 2);
            break;
        case 1155:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1155, 2)) + 0.5*pow(-my1155 + y1155, 2)/pow(sigmay1155, 2);
            break;
        case 1156:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1156, 2)) + 0.5*pow(-my1156 + y1156, 2)/pow(sigmay1156, 2);
            break;
        case 1157:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1157, 2)) + 0.5*pow(-my1157 + y1157, 2)/pow(sigmay1157, 2);
            break;
        case 1158:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1158, 2)) + 0.5*pow(-my1158 + y1158, 2)/pow(sigmay1158, 2);
            break;
        case 1159:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1159, 2)) + 0.5*pow(-my1159 + y1159, 2)/pow(sigmay1159, 2);
            break;
        case 1160:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1160, 2)) + 0.5*pow(-my1160 + y1160, 2)/pow(sigmay1160, 2);
            break;
        case 1161:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1161, 2)) + 0.5*pow(-my1161 + y1161, 2)/pow(sigmay1161, 2);
            break;
        case 1162:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1162, 2)) + 0.5*pow(-my1162 + y1162, 2)/pow(sigmay1162, 2);
            break;
        case 1163:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1163, 2)) + 0.5*pow(-my1163 + y1163, 2)/pow(sigmay1163, 2);
            break;
        case 1164:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1164, 2)) + 0.5*pow(-my1164 + y1164, 2)/pow(sigmay1164, 2);
            break;
        case 1165:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1165, 2)) + 0.5*pow(-my1165 + y1165, 2)/pow(sigmay1165, 2);
            break;
        case 1166:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1166, 2)) + 0.5*pow(-my1166 + y1166, 2)/pow(sigmay1166, 2);
            break;
        case 1167:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1167, 2)) + 0.5*pow(-my1167 + y1167, 2)/pow(sigmay1167, 2);
            break;
        case 1168:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1168, 2)) + 0.5*pow(-my1168 + y1168, 2)/pow(sigmay1168, 2);
            break;
        case 1169:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1169, 2)) + 0.5*pow(-my1169 + y1169, 2)/pow(sigmay1169, 2);
            break;
        case 1170:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1170, 2)) + 0.5*pow(-my1170 + y1170, 2)/pow(sigmay1170, 2);
            break;
        case 1171:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1171, 2)) + 0.5*pow(-my1171 + y1171, 2)/pow(sigmay1171, 2);
            break;
        case 1172:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1172, 2)) + 0.5*pow(-my1172 + y1172, 2)/pow(sigmay1172, 2);
            break;
        case 1173:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1173, 2)) + 0.5*pow(-my1173 + y1173, 2)/pow(sigmay1173, 2);
            break;
        case 1174:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1174, 2)) + 0.5*pow(-my1174 + y1174, 2)/pow(sigmay1174, 2);
            break;
        case 1175:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1175, 2)) + 0.5*pow(-my1175 + y1175, 2)/pow(sigmay1175, 2);
            break;
        case 1176:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1176, 2)) + 0.5*pow(-my1176 + y1176, 2)/pow(sigmay1176, 2);
            break;
        case 1177:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1177, 2)) + 0.5*pow(-my1177 + y1177, 2)/pow(sigmay1177, 2);
            break;
        case 1178:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1178, 2)) + 0.5*pow(-my1178 + y1178, 2)/pow(sigmay1178, 2);
            break;
        case 1179:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1179, 2)) + 0.5*pow(-my1179 + y1179, 2)/pow(sigmay1179, 2);
            break;
        case 1180:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1180, 2)) + 0.5*pow(-my1180 + y1180, 2)/pow(sigmay1180, 2);
            break;
        case 1181:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1181, 2)) + 0.5*pow(-my1181 + y1181, 2)/pow(sigmay1181, 2);
            break;
        case 1182:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1182, 2)) + 0.5*pow(-my1182 + y1182, 2)/pow(sigmay1182, 2);
            break;
        case 1183:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1183, 2)) + 0.5*pow(-my1183 + y1183, 2)/pow(sigmay1183, 2);
            break;
        case 1184:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1184, 2)) + 0.5*pow(-my1184 + y1184, 2)/pow(sigmay1184, 2);
            break;
        case 1185:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1185, 2)) + 0.5*pow(-my1185 + y1185, 2)/pow(sigmay1185, 2);
            break;
        case 1186:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1186, 2)) + 0.5*pow(-my1186 + y1186, 2)/pow(sigmay1186, 2);
            break;
        case 1187:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1187, 2)) + 0.5*pow(-my1187 + y1187, 2)/pow(sigmay1187, 2);
            break;
        case 1188:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1188, 2)) + 0.5*pow(-my1188 + y1188, 2)/pow(sigmay1188, 2);
            break;
        case 1189:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1189, 2)) + 0.5*pow(-my1189 + y1189, 2)/pow(sigmay1189, 2);
            break;
        case 1190:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1190, 2)) + 0.5*pow(-my1190 + y1190, 2)/pow(sigmay1190, 2);
            break;
        case 1191:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1191, 2)) + 0.5*pow(-my1191 + y1191, 2)/pow(sigmay1191, 2);
            break;
        case 1192:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1192, 2)) + 0.5*pow(-my1192 + y1192, 2)/pow(sigmay1192, 2);
            break;
        case 1193:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1193, 2)) + 0.5*pow(-my1193 + y1193, 2)/pow(sigmay1193, 2);
            break;
        case 1194:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1194, 2)) + 0.5*pow(-my1194 + y1194, 2)/pow(sigmay1194, 2);
            break;
        case 1195:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1195, 2)) + 0.5*pow(-my1195 + y1195, 2)/pow(sigmay1195, 2);
            break;
        case 1196:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1196, 2)) + 0.5*pow(-my1196 + y1196, 2)/pow(sigmay1196, 2);
            break;
        case 1197:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1197, 2)) + 0.5*pow(-my1197 + y1197, 2)/pow(sigmay1197, 2);
            break;
        case 1198:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1198, 2)) + 0.5*pow(-my1198 + y1198, 2)/pow(sigmay1198, 2);
            break;
        case 1199:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1199, 2)) + 0.5*pow(-my1199 + y1199, 2)/pow(sigmay1199, 2);
            break;
        case 1200:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1200, 2)) + 0.5*pow(-my1200 + y1200, 2)/pow(sigmay1200, 2);
            break;
        case 1201:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1201, 2)) + 0.5*pow(-my1201 + y1201, 2)/pow(sigmay1201, 2);
            break;
        case 1202:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1202, 2)) + 0.5*pow(-my1202 + y1202, 2)/pow(sigmay1202, 2);
            break;
        case 1203:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1203, 2)) + 0.5*pow(-my1203 + y1203, 2)/pow(sigmay1203, 2);
            break;
        case 1204:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1204, 2)) + 0.5*pow(-my1204 + y1204, 2)/pow(sigmay1204, 2);
            break;
        case 1205:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1205, 2)) + 0.5*pow(-my1205 + y1205, 2)/pow(sigmay1205, 2);
            break;
        case 1206:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1206, 2)) + 0.5*pow(-my1206 + y1206, 2)/pow(sigmay1206, 2);
            break;
        case 1207:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1207, 2)) + 0.5*pow(-my1207 + y1207, 2)/pow(sigmay1207, 2);
            break;
        case 1208:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1208, 2)) + 0.5*pow(-my1208 + y1208, 2)/pow(sigmay1208, 2);
            break;
        case 1209:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1209, 2)) + 0.5*pow(-my1209 + y1209, 2)/pow(sigmay1209, 2);
            break;
        case 1210:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1210, 2)) + 0.5*pow(-my1210 + y1210, 2)/pow(sigmay1210, 2);
            break;
        case 1211:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1211, 2)) + 0.5*pow(-my1211 + y1211, 2)/pow(sigmay1211, 2);
            break;
        case 1212:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1212, 2)) + 0.5*pow(-my1212 + y1212, 2)/pow(sigmay1212, 2);
            break;
        case 1213:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1213, 2)) + 0.5*pow(-my1213 + y1213, 2)/pow(sigmay1213, 2);
            break;
        case 1214:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1214, 2)) + 0.5*pow(-my1214 + y1214, 2)/pow(sigmay1214, 2);
            break;
        case 1215:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1215, 2)) + 0.5*pow(-my1215 + y1215, 2)/pow(sigmay1215, 2);
            break;
        case 1216:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1216, 2)) + 0.5*pow(-my1216 + y1216, 2)/pow(sigmay1216, 2);
            break;
        case 1217:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1217, 2)) + 0.5*pow(-my1217 + y1217, 2)/pow(sigmay1217, 2);
            break;
        case 1218:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1218, 2)) + 0.5*pow(-my1218 + y1218, 2)/pow(sigmay1218, 2);
            break;
        case 1219:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1219, 2)) + 0.5*pow(-my1219 + y1219, 2)/pow(sigmay1219, 2);
            break;
        case 1220:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1220, 2)) + 0.5*pow(-my1220 + y1220, 2)/pow(sigmay1220, 2);
            break;
        case 1221:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1221, 2)) + 0.5*pow(-my1221 + y1221, 2)/pow(sigmay1221, 2);
            break;
        case 1222:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1222, 2)) + 0.5*pow(-my1222 + y1222, 2)/pow(sigmay1222, 2);
            break;
        case 1223:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1223, 2)) + 0.5*pow(-my1223 + y1223, 2)/pow(sigmay1223, 2);
            break;
        case 1224:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1224, 2)) + 0.5*pow(-my1224 + y1224, 2)/pow(sigmay1224, 2);
            break;
        case 1225:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1225, 2)) + 0.5*pow(-my1225 + y1225, 2)/pow(sigmay1225, 2);
            break;
        case 1226:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1226, 2)) + 0.5*pow(-my1226 + y1226, 2)/pow(sigmay1226, 2);
            break;
        case 1227:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1227, 2)) + 0.5*pow(-my1227 + y1227, 2)/pow(sigmay1227, 2);
            break;
        case 1228:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1228, 2)) + 0.5*pow(-my1228 + y1228, 2)/pow(sigmay1228, 2);
            break;
        case 1229:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1229, 2)) + 0.5*pow(-my1229 + y1229, 2)/pow(sigmay1229, 2);
            break;
        case 1230:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1230, 2)) + 0.5*pow(-my1230 + y1230, 2)/pow(sigmay1230, 2);
            break;
        case 1231:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1231, 2)) + 0.5*pow(-my1231 + y1231, 2)/pow(sigmay1231, 2);
            break;
        case 1232:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1232, 2)) + 0.5*pow(-my1232 + y1232, 2)/pow(sigmay1232, 2);
            break;
        case 1233:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1233, 2)) + 0.5*pow(-my1233 + y1233, 2)/pow(sigmay1233, 2);
            break;
        case 1234:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1234, 2)) + 0.5*pow(-my1234 + y1234, 2)/pow(sigmay1234, 2);
            break;
        case 1235:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1235, 2)) + 0.5*pow(-my1235 + y1235, 2)/pow(sigmay1235, 2);
            break;
        case 1236:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1236, 2)) + 0.5*pow(-my1236 + y1236, 2)/pow(sigmay1236, 2);
            break;
        case 1237:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1237, 2)) + 0.5*pow(-my1237 + y1237, 2)/pow(sigmay1237, 2);
            break;
        case 1238:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1238, 2)) + 0.5*pow(-my1238 + y1238, 2)/pow(sigmay1238, 2);
            break;
        case 1239:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1239, 2)) + 0.5*pow(-my1239 + y1239, 2)/pow(sigmay1239, 2);
            break;
        case 1240:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1240, 2)) + 0.5*pow(-my1240 + y1240, 2)/pow(sigmay1240, 2);
            break;
        case 1241:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1241, 2)) + 0.5*pow(-my1241 + y1241, 2)/pow(sigmay1241, 2);
            break;
        case 1242:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1242, 2)) + 0.5*pow(-my1242 + y1242, 2)/pow(sigmay1242, 2);
            break;
        case 1243:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1243, 2)) + 0.5*pow(-my1243 + y1243, 2)/pow(sigmay1243, 2);
            break;
        case 1244:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1244, 2)) + 0.5*pow(-my1244 + y1244, 2)/pow(sigmay1244, 2);
            break;
        case 1245:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1245, 2)) + 0.5*pow(-my1245 + y1245, 2)/pow(sigmay1245, 2);
            break;
        case 1246:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1246, 2)) + 0.5*pow(-my1246 + y1246, 2)/pow(sigmay1246, 2);
            break;
        case 1247:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1247, 2)) + 0.5*pow(-my1247 + y1247, 2)/pow(sigmay1247, 2);
            break;
        case 1248:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1248, 2)) + 0.5*pow(-my1248 + y1248, 2)/pow(sigmay1248, 2);
            break;
        case 1249:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1249, 2)) + 0.5*pow(-my1249 + y1249, 2)/pow(sigmay1249, 2);
            break;
        case 1250:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1250, 2)) + 0.5*pow(-my1250 + y1250, 2)/pow(sigmay1250, 2);
            break;
        case 1251:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1251, 2)) + 0.5*pow(-my1251 + y1251, 2)/pow(sigmay1251, 2);
            break;
        case 1252:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1252, 2)) + 0.5*pow(-my1252 + y1252, 2)/pow(sigmay1252, 2);
            break;
        case 1253:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1253, 2)) + 0.5*pow(-my1253 + y1253, 2)/pow(sigmay1253, 2);
            break;
        case 1254:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1254, 2)) + 0.5*pow(-my1254 + y1254, 2)/pow(sigmay1254, 2);
            break;
        case 1255:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1255, 2)) + 0.5*pow(-my1255 + y1255, 2)/pow(sigmay1255, 2);
            break;
        case 1256:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1256, 2)) + 0.5*pow(-my1256 + y1256, 2)/pow(sigmay1256, 2);
            break;
        case 1257:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1257, 2)) + 0.5*pow(-my1257 + y1257, 2)/pow(sigmay1257, 2);
            break;
        case 1258:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1258, 2)) + 0.5*pow(-my1258 + y1258, 2)/pow(sigmay1258, 2);
            break;
        case 1259:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1259, 2)) + 0.5*pow(-my1259 + y1259, 2)/pow(sigmay1259, 2);
            break;
        case 1260:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1260, 2)) + 0.5*pow(-my1260 + y1260, 2)/pow(sigmay1260, 2);
            break;
        case 1261:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1261, 2)) + 0.5*pow(-my1261 + y1261, 2)/pow(sigmay1261, 2);
            break;
        case 1262:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1262, 2)) + 0.5*pow(-my1262 + y1262, 2)/pow(sigmay1262, 2);
            break;
        case 1263:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1263, 2)) + 0.5*pow(-my1263 + y1263, 2)/pow(sigmay1263, 2);
            break;
        case 1264:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1264, 2)) + 0.5*pow(-my1264 + y1264, 2)/pow(sigmay1264, 2);
            break;
        case 1265:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1265, 2)) + 0.5*pow(-my1265 + y1265, 2)/pow(sigmay1265, 2);
            break;
        case 1266:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1266, 2)) + 0.5*pow(-my1266 + y1266, 2)/pow(sigmay1266, 2);
            break;
        case 1267:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1267, 2)) + 0.5*pow(-my1267 + y1267, 2)/pow(sigmay1267, 2);
            break;
        case 1268:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1268, 2)) + 0.5*pow(-my1268 + y1268, 2)/pow(sigmay1268, 2);
            break;
        case 1269:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1269, 2)) + 0.5*pow(-my1269 + y1269, 2)/pow(sigmay1269, 2);
            break;
        case 1270:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1270, 2)) + 0.5*pow(-my1270 + y1270, 2)/pow(sigmay1270, 2);
            break;
        case 1271:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1271, 2)) + 0.5*pow(-my1271 + y1271, 2)/pow(sigmay1271, 2);
            break;
        case 1272:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1272, 2)) + 0.5*pow(-my1272 + y1272, 2)/pow(sigmay1272, 2);
            break;
        case 1273:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1273, 2)) + 0.5*pow(-my1273 + y1273, 2)/pow(sigmay1273, 2);
            break;
        case 1274:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1274, 2)) + 0.5*pow(-my1274 + y1274, 2)/pow(sigmay1274, 2);
            break;
        case 1275:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1275, 2)) + 0.5*pow(-my1275 + y1275, 2)/pow(sigmay1275, 2);
            break;
        case 1276:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1276, 2)) + 0.5*pow(-my1276 + y1276, 2)/pow(sigmay1276, 2);
            break;
        case 1277:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1277, 2)) + 0.5*pow(-my1277 + y1277, 2)/pow(sigmay1277, 2);
            break;
        case 1278:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1278, 2)) + 0.5*pow(-my1278 + y1278, 2)/pow(sigmay1278, 2);
            break;
        case 1279:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1279, 2)) + 0.5*pow(-my1279 + y1279, 2)/pow(sigmay1279, 2);
            break;
        case 1280:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1280, 2)) + 0.5*pow(-my1280 + y1280, 2)/pow(sigmay1280, 2);
            break;
        case 1281:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1281, 2)) + 0.5*pow(-my1281 + y1281, 2)/pow(sigmay1281, 2);
            break;
        case 1282:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1282, 2)) + 0.5*pow(-my1282 + y1282, 2)/pow(sigmay1282, 2);
            break;
        case 1283:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1283, 2)) + 0.5*pow(-my1283 + y1283, 2)/pow(sigmay1283, 2);
            break;
        case 1284:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1284, 2)) + 0.5*pow(-my1284 + y1284, 2)/pow(sigmay1284, 2);
            break;
        case 1285:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1285, 2)) + 0.5*pow(-my1285 + y1285, 2)/pow(sigmay1285, 2);
            break;
        case 1286:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1286, 2)) + 0.5*pow(-my1286 + y1286, 2)/pow(sigmay1286, 2);
            break;
        case 1287:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1287, 2)) + 0.5*pow(-my1287 + y1287, 2)/pow(sigmay1287, 2);
            break;
        case 1288:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1288, 2)) + 0.5*pow(-my1288 + y1288, 2)/pow(sigmay1288, 2);
            break;
        case 1289:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1289, 2)) + 0.5*pow(-my1289 + y1289, 2)/pow(sigmay1289, 2);
            break;
        case 1290:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1290, 2)) + 0.5*pow(-my1290 + y1290, 2)/pow(sigmay1290, 2);
            break;
        case 1291:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1291, 2)) + 0.5*pow(-my1291 + y1291, 2)/pow(sigmay1291, 2);
            break;
        case 1292:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1292, 2)) + 0.5*pow(-my1292 + y1292, 2)/pow(sigmay1292, 2);
            break;
        case 1293:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1293, 2)) + 0.5*pow(-my1293 + y1293, 2)/pow(sigmay1293, 2);
            break;
        case 1294:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1294, 2)) + 0.5*pow(-my1294 + y1294, 2)/pow(sigmay1294, 2);
            break;
        case 1295:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1295, 2)) + 0.5*pow(-my1295 + y1295, 2)/pow(sigmay1295, 2);
            break;
        case 1296:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1296, 2)) + 0.5*pow(-my1296 + y1296, 2)/pow(sigmay1296, 2);
            break;
        case 1297:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1297, 2)) + 0.5*pow(-my1297 + y1297, 2)/pow(sigmay1297, 2);
            break;
        case 1298:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1298, 2)) + 0.5*pow(-my1298 + y1298, 2)/pow(sigmay1298, 2);
            break;
        case 1299:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1299, 2)) + 0.5*pow(-my1299 + y1299, 2)/pow(sigmay1299, 2);
            break;
        case 1300:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1300, 2)) + 0.5*pow(-my1300 + y1300, 2)/pow(sigmay1300, 2);
            break;
        case 1301:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1301, 2)) + 0.5*pow(-my1301 + y1301, 2)/pow(sigmay1301, 2);
            break;
        case 1302:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1302, 2)) + 0.5*pow(-my1302 + y1302, 2)/pow(sigmay1302, 2);
            break;
        case 1303:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1303, 2)) + 0.5*pow(-my1303 + y1303, 2)/pow(sigmay1303, 2);
            break;
        case 1304:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1304, 2)) + 0.5*pow(-my1304 + y1304, 2)/pow(sigmay1304, 2);
            break;
        case 1305:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1305, 2)) + 0.5*pow(-my1305 + y1305, 2)/pow(sigmay1305, 2);
            break;
        case 1306:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1306, 2)) + 0.5*pow(-my1306 + y1306, 2)/pow(sigmay1306, 2);
            break;
        case 1307:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1307, 2)) + 0.5*pow(-my1307 + y1307, 2)/pow(sigmay1307, 2);
            break;
        case 1308:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1308, 2)) + 0.5*pow(-my1308 + y1308, 2)/pow(sigmay1308, 2);
            break;
        case 1309:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1309, 2)) + 0.5*pow(-my1309 + y1309, 2)/pow(sigmay1309, 2);
            break;
        case 1310:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1310, 2)) + 0.5*pow(-my1310 + y1310, 2)/pow(sigmay1310, 2);
            break;
        case 1311:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1311, 2)) + 0.5*pow(-my1311 + y1311, 2)/pow(sigmay1311, 2);
            break;
        case 1312:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1312, 2)) + 0.5*pow(-my1312 + y1312, 2)/pow(sigmay1312, 2);
            break;
        case 1313:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1313, 2)) + 0.5*pow(-my1313 + y1313, 2)/pow(sigmay1313, 2);
            break;
        case 1314:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1314, 2)) + 0.5*pow(-my1314 + y1314, 2)/pow(sigmay1314, 2);
            break;
        case 1315:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1315, 2)) + 0.5*pow(-my1315 + y1315, 2)/pow(sigmay1315, 2);
            break;
        case 1316:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1316, 2)) + 0.5*pow(-my1316 + y1316, 2)/pow(sigmay1316, 2);
            break;
        case 1317:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1317, 2)) + 0.5*pow(-my1317 + y1317, 2)/pow(sigmay1317, 2);
            break;
        case 1318:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1318, 2)) + 0.5*pow(-my1318 + y1318, 2)/pow(sigmay1318, 2);
            break;
        case 1319:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1319, 2)) + 0.5*pow(-my1319 + y1319, 2)/pow(sigmay1319, 2);
            break;
        case 1320:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1320, 2)) + 0.5*pow(-my1320 + y1320, 2)/pow(sigmay1320, 2);
            break;
        case 1321:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1321, 2)) + 0.5*pow(-my1321 + y1321, 2)/pow(sigmay1321, 2);
            break;
        case 1322:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1322, 2)) + 0.5*pow(-my1322 + y1322, 2)/pow(sigmay1322, 2);
            break;
        case 1323:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1323, 2)) + 0.5*pow(-my1323 + y1323, 2)/pow(sigmay1323, 2);
            break;
        case 1324:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1324, 2)) + 0.5*pow(-my1324 + y1324, 2)/pow(sigmay1324, 2);
            break;
        case 1325:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1325, 2)) + 0.5*pow(-my1325 + y1325, 2)/pow(sigmay1325, 2);
            break;
        case 1326:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1326, 2)) + 0.5*pow(-my1326 + y1326, 2)/pow(sigmay1326, 2);
            break;
        case 1327:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1327, 2)) + 0.5*pow(-my1327 + y1327, 2)/pow(sigmay1327, 2);
            break;
        case 1328:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1328, 2)) + 0.5*pow(-my1328 + y1328, 2)/pow(sigmay1328, 2);
            break;
        case 1329:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1329, 2)) + 0.5*pow(-my1329 + y1329, 2)/pow(sigmay1329, 2);
            break;
        case 1330:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1330, 2)) + 0.5*pow(-my1330 + y1330, 2)/pow(sigmay1330, 2);
            break;
        case 1331:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1331, 2)) + 0.5*pow(-my1331 + y1331, 2)/pow(sigmay1331, 2);
            break;
        case 1332:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1332, 2)) + 0.5*pow(-my1332 + y1332, 2)/pow(sigmay1332, 2);
            break;
        case 1333:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1333, 2)) + 0.5*pow(-my1333 + y1333, 2)/pow(sigmay1333, 2);
            break;
        case 1334:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1334, 2)) + 0.5*pow(-my1334 + y1334, 2)/pow(sigmay1334, 2);
            break;
        case 1335:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1335, 2)) + 0.5*pow(-my1335 + y1335, 2)/pow(sigmay1335, 2);
            break;
        case 1336:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1336, 2)) + 0.5*pow(-my1336 + y1336, 2)/pow(sigmay1336, 2);
            break;
        case 1337:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1337, 2)) + 0.5*pow(-my1337 + y1337, 2)/pow(sigmay1337, 2);
            break;
        case 1338:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1338, 2)) + 0.5*pow(-my1338 + y1338, 2)/pow(sigmay1338, 2);
            break;
        case 1339:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1339, 2)) + 0.5*pow(-my1339 + y1339, 2)/pow(sigmay1339, 2);
            break;
        case 1340:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1340, 2)) + 0.5*pow(-my1340 + y1340, 2)/pow(sigmay1340, 2);
            break;
        case 1341:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1341, 2)) + 0.5*pow(-my1341 + y1341, 2)/pow(sigmay1341, 2);
            break;
        case 1342:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1342, 2)) + 0.5*pow(-my1342 + y1342, 2)/pow(sigmay1342, 2);
            break;
        case 1343:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1343, 2)) + 0.5*pow(-my1343 + y1343, 2)/pow(sigmay1343, 2);
            break;
        case 1344:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1344, 2)) + 0.5*pow(-my1344 + y1344, 2)/pow(sigmay1344, 2);
            break;
        case 1345:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1345, 2)) + 0.5*pow(-my1345 + y1345, 2)/pow(sigmay1345, 2);
            break;
        case 1346:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1346, 2)) + 0.5*pow(-my1346 + y1346, 2)/pow(sigmay1346, 2);
            break;
        case 1347:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1347, 2)) + 0.5*pow(-my1347 + y1347, 2)/pow(sigmay1347, 2);
            break;
        case 1348:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1348, 2)) + 0.5*pow(-my1348 + y1348, 2)/pow(sigmay1348, 2);
            break;
        case 1349:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1349, 2)) + 0.5*pow(-my1349 + y1349, 2)/pow(sigmay1349, 2);
            break;
        case 1350:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1350, 2)) + 0.5*pow(-my1350 + y1350, 2)/pow(sigmay1350, 2);
            break;
        case 1351:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1351, 2)) + 0.5*pow(-my1351 + y1351, 2)/pow(sigmay1351, 2);
            break;
        case 1352:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1352, 2)) + 0.5*pow(-my1352 + y1352, 2)/pow(sigmay1352, 2);
            break;
        case 1353:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1353, 2)) + 0.5*pow(-my1353 + y1353, 2)/pow(sigmay1353, 2);
            break;
        case 1354:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1354, 2)) + 0.5*pow(-my1354 + y1354, 2)/pow(sigmay1354, 2);
            break;
        case 1355:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1355, 2)) + 0.5*pow(-my1355 + y1355, 2)/pow(sigmay1355, 2);
            break;
        case 1356:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1356, 2)) + 0.5*pow(-my1356 + y1356, 2)/pow(sigmay1356, 2);
            break;
        case 1357:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1357, 2)) + 0.5*pow(-my1357 + y1357, 2)/pow(sigmay1357, 2);
            break;
        case 1358:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1358, 2)) + 0.5*pow(-my1358 + y1358, 2)/pow(sigmay1358, 2);
            break;
        case 1359:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1359, 2)) + 0.5*pow(-my1359 + y1359, 2)/pow(sigmay1359, 2);
            break;
        case 1360:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1360, 2)) + 0.5*pow(-my1360 + y1360, 2)/pow(sigmay1360, 2);
            break;
        case 1361:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1361, 2)) + 0.5*pow(-my1361 + y1361, 2)/pow(sigmay1361, 2);
            break;
        case 1362:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1362, 2)) + 0.5*pow(-my1362 + y1362, 2)/pow(sigmay1362, 2);
            break;
        case 1363:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1363, 2)) + 0.5*pow(-my1363 + y1363, 2)/pow(sigmay1363, 2);
            break;
        case 1364:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1364, 2)) + 0.5*pow(-my1364 + y1364, 2)/pow(sigmay1364, 2);
            break;
        case 1365:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1365, 2)) + 0.5*pow(-my1365 + y1365, 2)/pow(sigmay1365, 2);
            break;
        case 1366:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1366, 2)) + 0.5*pow(-my1366 + y1366, 2)/pow(sigmay1366, 2);
            break;
        case 1367:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1367, 2)) + 0.5*pow(-my1367 + y1367, 2)/pow(sigmay1367, 2);
            break;
        case 1368:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1368, 2)) + 0.5*pow(-my1368 + y1368, 2)/pow(sigmay1368, 2);
            break;
        case 1369:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1369, 2)) + 0.5*pow(-my1369 + y1369, 2)/pow(sigmay1369, 2);
            break;
        case 1370:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1370, 2)) + 0.5*pow(-my1370 + y1370, 2)/pow(sigmay1370, 2);
            break;
        case 1371:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1371, 2)) + 0.5*pow(-my1371 + y1371, 2)/pow(sigmay1371, 2);
            break;
        case 1372:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1372, 2)) + 0.5*pow(-my1372 + y1372, 2)/pow(sigmay1372, 2);
            break;
        case 1373:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1373, 2)) + 0.5*pow(-my1373 + y1373, 2)/pow(sigmay1373, 2);
            break;
        case 1374:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1374, 2)) + 0.5*pow(-my1374 + y1374, 2)/pow(sigmay1374, 2);
            break;
        case 1375:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1375, 2)) + 0.5*pow(-my1375 + y1375, 2)/pow(sigmay1375, 2);
            break;
        case 1376:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1376, 2)) + 0.5*pow(-my1376 + y1376, 2)/pow(sigmay1376, 2);
            break;
        case 1377:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1377, 2)) + 0.5*pow(-my1377 + y1377, 2)/pow(sigmay1377, 2);
            break;
        case 1378:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1378, 2)) + 0.5*pow(-my1378 + y1378, 2)/pow(sigmay1378, 2);
            break;
        case 1379:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1379, 2)) + 0.5*pow(-my1379 + y1379, 2)/pow(sigmay1379, 2);
            break;
        case 1380:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1380, 2)) + 0.5*pow(-my1380 + y1380, 2)/pow(sigmay1380, 2);
            break;
        case 1381:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1381, 2)) + 0.5*pow(-my1381 + y1381, 2)/pow(sigmay1381, 2);
            break;
        case 1382:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1382, 2)) + 0.5*pow(-my1382 + y1382, 2)/pow(sigmay1382, 2);
            break;
        case 1383:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1383, 2)) + 0.5*pow(-my1383 + y1383, 2)/pow(sigmay1383, 2);
            break;
        case 1384:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1384, 2)) + 0.5*pow(-my1384 + y1384, 2)/pow(sigmay1384, 2);
            break;
        case 1385:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1385, 2)) + 0.5*pow(-my1385 + y1385, 2)/pow(sigmay1385, 2);
            break;
        case 1386:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1386, 2)) + 0.5*pow(-my1386 + y1386, 2)/pow(sigmay1386, 2);
            break;
        case 1387:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1387, 2)) + 0.5*pow(-my1387 + y1387, 2)/pow(sigmay1387, 2);
            break;
        case 1388:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1388, 2)) + 0.5*pow(-my1388 + y1388, 2)/pow(sigmay1388, 2);
            break;
        case 1389:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1389, 2)) + 0.5*pow(-my1389 + y1389, 2)/pow(sigmay1389, 2);
            break;
        case 1390:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1390, 2)) + 0.5*pow(-my1390 + y1390, 2)/pow(sigmay1390, 2);
            break;
        case 1391:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1391, 2)) + 0.5*pow(-my1391 + y1391, 2)/pow(sigmay1391, 2);
            break;
        case 1392:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1392, 2)) + 0.5*pow(-my1392 + y1392, 2)/pow(sigmay1392, 2);
            break;
        case 1393:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1393, 2)) + 0.5*pow(-my1393 + y1393, 2)/pow(sigmay1393, 2);
            break;
        case 1394:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1394, 2)) + 0.5*pow(-my1394 + y1394, 2)/pow(sigmay1394, 2);
            break;
        case 1395:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1395, 2)) + 0.5*pow(-my1395 + y1395, 2)/pow(sigmay1395, 2);
            break;
    }
}