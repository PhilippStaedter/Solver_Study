#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydy_Proctor2010a(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = 0.5*(-2*my0 + 2*y0)/pow(sigmay0, 2);
            break;
        case 1:
            dJydy[0] = 0.5*(-2*my1 + 2*y1)/pow(sigmay1, 2);
            break;
        case 2:
            dJydy[0] = 0.5*(-2*my2 + 2*y2)/pow(sigmay2, 2);
            break;
        case 3:
            dJydy[0] = 0.5*(-2*my3 + 2*y3)/pow(sigmay3, 2);
            break;
        case 4:
            dJydy[0] = 0.5*(-2*my4 + 2*y4)/pow(sigmay4, 2);
            break;
        case 5:
            dJydy[0] = 0.5*(-2*my5 + 2*y5)/pow(sigmay5, 2);
            break;
        case 6:
            dJydy[0] = 0.5*(-2*my6 + 2*y6)/pow(sigmay6, 2);
            break;
        case 7:
            dJydy[0] = 0.5*(-2*my7 + 2*y7)/pow(sigmay7, 2);
            break;
        case 8:
            dJydy[0] = 0.5*(-2*my8 + 2*y8)/pow(sigmay8, 2);
            break;
        case 9:
            dJydy[0] = 0.5*(-2*my9 + 2*y9)/pow(sigmay9, 2);
            break;
        case 10:
            dJydy[0] = 0.5*(-2*my10 + 2*y10)/pow(sigmay10, 2);
            break;
        case 11:
            dJydy[0] = 0.5*(-2*my11 + 2*y11)/pow(sigmay11, 2);
            break;
        case 12:
            dJydy[0] = 0.5*(-2*my12 + 2*y12)/pow(sigmay12, 2);
            break;
        case 13:
            dJydy[0] = 0.5*(-2*my13 + 2*y13)/pow(sigmay13, 2);
            break;
        case 14:
            dJydy[0] = 0.5*(-2*my14 + 2*y14)/pow(sigmay14, 2);
            break;
        case 15:
            dJydy[0] = 0.5*(-2*my15 + 2*y15)/pow(sigmay15, 2);
            break;
        case 16:
            dJydy[0] = 0.5*(-2*my16 + 2*y16)/pow(sigmay16, 2);
            break;
        case 17:
            dJydy[0] = 0.5*(-2*my17 + 2*y17)/pow(sigmay17, 2);
            break;
        case 18:
            dJydy[0] = 0.5*(-2*my18 + 2*y18)/pow(sigmay18, 2);
            break;
        case 19:
            dJydy[0] = 0.5*(-2*my19 + 2*y19)/pow(sigmay19, 2);
            break;
        case 20:
            dJydy[0] = 0.5*(-2*my20 + 2*y20)/pow(sigmay20, 2);
            break;
        case 21:
            dJydy[0] = 0.5*(-2*my21 + 2*y21)/pow(sigmay21, 2);
            break;
        case 22:
            dJydy[0] = 0.5*(-2*my22 + 2*y22)/pow(sigmay22, 2);
            break;
        case 23:
            dJydy[0] = 0.5*(-2*my23 + 2*y23)/pow(sigmay23, 2);
            break;
        case 24:
            dJydy[0] = 0.5*(-2*my24 + 2*y24)/pow(sigmay24, 2);
            break;
        case 25:
            dJydy[0] = 0.5*(-2*my25 + 2*y25)/pow(sigmay25, 2);
            break;
        case 26:
            dJydy[0] = 0.5*(-2*my26 + 2*y26)/pow(sigmay26, 2);
            break;
        case 27:
            dJydy[0] = 0.5*(-2*my27 + 2*y27)/pow(sigmay27, 2);
            break;
        case 28:
            dJydy[0] = 0.5*(-2*my28 + 2*y28)/pow(sigmay28, 2);
            break;
        case 29:
            dJydy[0] = 0.5*(-2*my29 + 2*y29)/pow(sigmay29, 2);
            break;
        case 30:
            dJydy[0] = 0.5*(-2*my30 + 2*y30)/pow(sigmay30, 2);
            break;
        case 31:
            dJydy[0] = 0.5*(-2*my31 + 2*y31)/pow(sigmay31, 2);
            break;
        case 32:
            dJydy[0] = 0.5*(-2*my32 + 2*y32)/pow(sigmay32, 2);
            break;
        case 33:
            dJydy[0] = 0.5*(-2*my33 + 2*y33)/pow(sigmay33, 2);
            break;
        case 34:
            dJydy[0] = 0.5*(-2*my34 + 2*y34)/pow(sigmay34, 2);
            break;
        case 35:
            dJydy[0] = 0.5*(-2*my35 + 2*y35)/pow(sigmay35, 2);
            break;
        case 36:
            dJydy[0] = 0.5*(-2*my36 + 2*y36)/pow(sigmay36, 2);
            break;
        case 37:
            dJydy[0] = 0.5*(-2*my37 + 2*y37)/pow(sigmay37, 2);
            break;
        case 38:
            dJydy[0] = 0.5*(-2*my38 + 2*y38)/pow(sigmay38, 2);
            break;
        case 39:
            dJydy[0] = 0.5*(-2*my39 + 2*y39)/pow(sigmay39, 2);
            break;
        case 40:
            dJydy[0] = 0.5*(-2*my40 + 2*y40)/pow(sigmay40, 2);
            break;
        case 41:
            dJydy[0] = 0.5*(-2*my41 + 2*y41)/pow(sigmay41, 2);
            break;
        case 42:
            dJydy[0] = 0.5*(-2*my42 + 2*y42)/pow(sigmay42, 2);
            break;
        case 43:
            dJydy[0] = 0.5*(-2*my43 + 2*y43)/pow(sigmay43, 2);
            break;
        case 44:
            dJydy[0] = 0.5*(-2*my44 + 2*y44)/pow(sigmay44, 2);
            break;
        case 45:
            dJydy[0] = 0.5*(-2*my45 + 2*y45)/pow(sigmay45, 2);
            break;
        case 46:
            dJydy[0] = 0.5*(-2*my46 + 2*y46)/pow(sigmay46, 2);
            break;
        case 47:
            dJydy[0] = 0.5*(-2*my47 + 2*y47)/pow(sigmay47, 2);
            break;
        case 48:
            dJydy[0] = 0.5*(-2*my48 + 2*y48)/pow(sigmay48, 2);
            break;
        case 49:
            dJydy[0] = 0.5*(-2*my49 + 2*y49)/pow(sigmay49, 2);
            break;
        case 50:
            dJydy[0] = 0.5*(-2*my50 + 2*y50)/pow(sigmay50, 2);
            break;
        case 51:
            dJydy[0] = 0.5*(-2*my51 + 2*y51)/pow(sigmay51, 2);
            break;
        case 52:
            dJydy[0] = 0.5*(-2*my52 + 2*y52)/pow(sigmay52, 2);
            break;
        case 53:
            dJydy[0] = 0.5*(-2*my53 + 2*y53)/pow(sigmay53, 2);
            break;
        case 54:
            dJydy[0] = 0.5*(-2*my54 + 2*y54)/pow(sigmay54, 2);
            break;
        case 55:
            dJydy[0] = 0.5*(-2*my55 + 2*y55)/pow(sigmay55, 2);
            break;
        case 56:
            dJydy[0] = 0.5*(-2*my56 + 2*y56)/pow(sigmay56, 2);
            break;
        case 57:
            dJydy[0] = 0.5*(-2*my57 + 2*y57)/pow(sigmay57, 2);
            break;
        case 58:
            dJydy[0] = 0.5*(-2*my58 + 2*y58)/pow(sigmay58, 2);
            break;
        case 59:
            dJydy[0] = 0.5*(-2*my59 + 2*y59)/pow(sigmay59, 2);
            break;
        case 60:
            dJydy[0] = 0.5*(-2*my60 + 2*y60)/pow(sigmay60, 2);
            break;
        case 61:
            dJydy[0] = 0.5*(-2*my61 + 2*y61)/pow(sigmay61, 2);
            break;
        case 62:
            dJydy[0] = 0.5*(-2*my62 + 2*y62)/pow(sigmay62, 2);
            break;
        case 63:
            dJydy[0] = 0.5*(-2*my63 + 2*y63)/pow(sigmay63, 2);
            break;
        case 64:
            dJydy[0] = 0.5*(-2*my64 + 2*y64)/pow(sigmay64, 2);
            break;
        case 65:
            dJydy[0] = 0.5*(-2*my65 + 2*y65)/pow(sigmay65, 2);
            break;
        case 66:
            dJydy[0] = 0.5*(-2*my66 + 2*y66)/pow(sigmay66, 2);
            break;
        case 67:
            dJydy[0] = 0.5*(-2*my67 + 2*y67)/pow(sigmay67, 2);
            break;
        case 68:
            dJydy[0] = 0.5*(-2*my68 + 2*y68)/pow(sigmay68, 2);
            break;
        case 69:
            dJydy[0] = 0.5*(-2*my69 + 2*y69)/pow(sigmay69, 2);
            break;
        case 70:
            dJydy[0] = 0.5*(-2*my70 + 2*y70)/pow(sigmay70, 2);
            break;
        case 71:
            dJydy[0] = 0.5*(-2*my71 + 2*y71)/pow(sigmay71, 2);
            break;
        case 72:
            dJydy[0] = 0.5*(-2*my72 + 2*y72)/pow(sigmay72, 2);
            break;
        case 73:
            dJydy[0] = 0.5*(-2*my73 + 2*y73)/pow(sigmay73, 2);
            break;
        case 74:
            dJydy[0] = 0.5*(-2*my74 + 2*y74)/pow(sigmay74, 2);
            break;
        case 75:
            dJydy[0] = 0.5*(-2*my75 + 2*y75)/pow(sigmay75, 2);
            break;
        case 76:
            dJydy[0] = 0.5*(-2*my76 + 2*y76)/pow(sigmay76, 2);
            break;
        case 77:
            dJydy[0] = 0.5*(-2*my77 + 2*y77)/pow(sigmay77, 2);
            break;
        case 78:
            dJydy[0] = 0.5*(-2*my78 + 2*y78)/pow(sigmay78, 2);
            break;
        case 79:
            dJydy[0] = 0.5*(-2*my79 + 2*y79)/pow(sigmay79, 2);
            break;
        case 80:
            dJydy[0] = 0.5*(-2*my80 + 2*y80)/pow(sigmay80, 2);
            break;
        case 81:
            dJydy[0] = 0.5*(-2*my81 + 2*y81)/pow(sigmay81, 2);
            break;
        case 82:
            dJydy[0] = 0.5*(-2*my82 + 2*y82)/pow(sigmay82, 2);
            break;
        case 83:
            dJydy[0] = 0.5*(-2*my83 + 2*y83)/pow(sigmay83, 2);
            break;
        case 84:
            dJydy[0] = 0.5*(-2*my84 + 2*y84)/pow(sigmay84, 2);
            break;
        case 85:
            dJydy[0] = 0.5*(-2*my85 + 2*y85)/pow(sigmay85, 2);
            break;
        case 86:
            dJydy[0] = 0.5*(-2*my86 + 2*y86)/pow(sigmay86, 2);
            break;
        case 87:
            dJydy[0] = 0.5*(-2*my87 + 2*y87)/pow(sigmay87, 2);
            break;
        case 88:
            dJydy[0] = 0.5*(-2*my88 + 2*y88)/pow(sigmay88, 2);
            break;
        case 89:
            dJydy[0] = 0.5*(-2*my89 + 2*y89)/pow(sigmay89, 2);
            break;
        case 90:
            dJydy[0] = 0.5*(-2*my90 + 2*y90)/pow(sigmay90, 2);
            break;
        case 91:
            dJydy[0] = 0.5*(-2*my91 + 2*y91)/pow(sigmay91, 2);
            break;
        case 92:
            dJydy[0] = 0.5*(-2*my92 + 2*y92)/pow(sigmay92, 2);
            break;
        case 93:
            dJydy[0] = 0.5*(-2*my93 + 2*y93)/pow(sigmay93, 2);
            break;
        case 94:
            dJydy[0] = 0.5*(-2*my94 + 2*y94)/pow(sigmay94, 2);
            break;
        case 95:
            dJydy[0] = 0.5*(-2*my95 + 2*y95)/pow(sigmay95, 2);
            break;
        case 96:
            dJydy[0] = 0.5*(-2*my96 + 2*y96)/pow(sigmay96, 2);
            break;
        case 97:
            dJydy[0] = 0.5*(-2*my97 + 2*y97)/pow(sigmay97, 2);
            break;
        case 98:
            dJydy[0] = 0.5*(-2*my98 + 2*y98)/pow(sigmay98, 2);
            break;
        case 99:
            dJydy[0] = 0.5*(-2*my99 + 2*y99)/pow(sigmay99, 2);
            break;
        case 100:
            dJydy[0] = 0.5*(-2*my100 + 2*y100)/pow(sigmay100, 2);
            break;
        case 101:
            dJydy[0] = 0.5*(-2*my101 + 2*y101)/pow(sigmay101, 2);
            break;
        case 102:
            dJydy[0] = 0.5*(-2*my102 + 2*y102)/pow(sigmay102, 2);
            break;
        case 103:
            dJydy[0] = 0.5*(-2*my103 + 2*y103)/pow(sigmay103, 2);
            break;
        case 104:
            dJydy[0] = 0.5*(-2*my104 + 2*y104)/pow(sigmay104, 2);
            break;
        case 105:
            dJydy[0] = 0.5*(-2*my105 + 2*y105)/pow(sigmay105, 2);
            break;
        case 106:
            dJydy[0] = 0.5*(-2*my106 + 2*y106)/pow(sigmay106, 2);
            break;
        case 107:
            dJydy[0] = 0.5*(-2*my107 + 2*y107)/pow(sigmay107, 2);
            break;
        case 108:
            dJydy[0] = 0.5*(-2*my108 + 2*y108)/pow(sigmay108, 2);
            break;
        case 109:
            dJydy[0] = 0.5*(-2*my109 + 2*y109)/pow(sigmay109, 2);
            break;
        case 110:
            dJydy[0] = 0.5*(-2*my110 + 2*y110)/pow(sigmay110, 2);
            break;
        case 111:
            dJydy[0] = 0.5*(-2*my111 + 2*y111)/pow(sigmay111, 2);
            break;
        case 112:
            dJydy[0] = 0.5*(-2*my112 + 2*y112)/pow(sigmay112, 2);
            break;
        case 113:
            dJydy[0] = 0.5*(-2*my113 + 2*y113)/pow(sigmay113, 2);
            break;
        case 114:
            dJydy[0] = 0.5*(-2*my114 + 2*y114)/pow(sigmay114, 2);
            break;
        case 115:
            dJydy[0] = 0.5*(-2*my115 + 2*y115)/pow(sigmay115, 2);
            break;
        case 116:
            dJydy[0] = 0.5*(-2*my116 + 2*y116)/pow(sigmay116, 2);
            break;
        case 117:
            dJydy[0] = 0.5*(-2*my117 + 2*y117)/pow(sigmay117, 2);
            break;
        case 118:
            dJydy[0] = 0.5*(-2*my118 + 2*y118)/pow(sigmay118, 2);
            break;
        case 119:
            dJydy[0] = 0.5*(-2*my119 + 2*y119)/pow(sigmay119, 2);
            break;
        case 120:
            dJydy[0] = 0.5*(-2*my120 + 2*y120)/pow(sigmay120, 2);
            break;
        case 121:
            dJydy[0] = 0.5*(-2*my121 + 2*y121)/pow(sigmay121, 2);
            break;
        case 122:
            dJydy[0] = 0.5*(-2*my122 + 2*y122)/pow(sigmay122, 2);
            break;
        case 123:
            dJydy[0] = 0.5*(-2*my123 + 2*y123)/pow(sigmay123, 2);
            break;
        case 124:
            dJydy[0] = 0.5*(-2*my124 + 2*y124)/pow(sigmay124, 2);
            break;
        case 125:
            dJydy[0] = 0.5*(-2*my125 + 2*y125)/pow(sigmay125, 2);
            break;
        case 126:
            dJydy[0] = 0.5*(-2*my126 + 2*y126)/pow(sigmay126, 2);
            break;
        case 127:
            dJydy[0] = 0.5*(-2*my127 + 2*y127)/pow(sigmay127, 2);
            break;
        case 128:
            dJydy[0] = 0.5*(-2*my128 + 2*y128)/pow(sigmay128, 2);
            break;
        case 129:
            dJydy[0] = 0.5*(-2*my129 + 2*y129)/pow(sigmay129, 2);
            break;
        case 130:
            dJydy[0] = 0.5*(-2*my130 + 2*y130)/pow(sigmay130, 2);
            break;
        case 131:
            dJydy[0] = 0.5*(-2*my131 + 2*y131)/pow(sigmay131, 2);
            break;
        case 132:
            dJydy[0] = 0.5*(-2*my132 + 2*y132)/pow(sigmay132, 2);
            break;
        case 133:
            dJydy[0] = 0.5*(-2*my133 + 2*y133)/pow(sigmay133, 2);
            break;
        case 134:
            dJydy[0] = 0.5*(-2*my134 + 2*y134)/pow(sigmay134, 2);
            break;
        case 135:
            dJydy[0] = 0.5*(-2*my135 + 2*y135)/pow(sigmay135, 2);
            break;
    }
}