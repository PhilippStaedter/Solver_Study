#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydy_Froehlich2018(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
        case 136:
            dJydy[0] = 0.5*(-2*my136 + 2*y136)/pow(sigmay136, 2);
            break;
        case 137:
            dJydy[0] = 0.5*(-2*my137 + 2*y137)/pow(sigmay137, 2);
            break;
        case 138:
            dJydy[0] = 0.5*(-2*my138 + 2*y138)/pow(sigmay138, 2);
            break;
        case 139:
            dJydy[0] = 0.5*(-2*my139 + 2*y139)/pow(sigmay139, 2);
            break;
        case 140:
            dJydy[0] = 0.5*(-2*my140 + 2*y140)/pow(sigmay140, 2);
            break;
        case 141:
            dJydy[0] = 0.5*(-2*my141 + 2*y141)/pow(sigmay141, 2);
            break;
        case 142:
            dJydy[0] = 0.5*(-2*my142 + 2*y142)/pow(sigmay142, 2);
            break;
        case 143:
            dJydy[0] = 0.5*(-2*my143 + 2*y143)/pow(sigmay143, 2);
            break;
        case 144:
            dJydy[0] = 0.5*(-2*my144 + 2*y144)/pow(sigmay144, 2);
            break;
        case 145:
            dJydy[0] = 0.5*(-2*my145 + 2*y145)/pow(sigmay145, 2);
            break;
        case 146:
            dJydy[0] = 0.5*(-2*my146 + 2*y146)/pow(sigmay146, 2);
            break;
        case 147:
            dJydy[0] = 0.5*(-2*my147 + 2*y147)/pow(sigmay147, 2);
            break;
        case 148:
            dJydy[0] = 0.5*(-2*my148 + 2*y148)/pow(sigmay148, 2);
            break;
        case 149:
            dJydy[0] = 0.5*(-2*my149 + 2*y149)/pow(sigmay149, 2);
            break;
        case 150:
            dJydy[0] = 0.5*(-2*my150 + 2*y150)/pow(sigmay150, 2);
            break;
        case 151:
            dJydy[0] = 0.5*(-2*my151 + 2*y151)/pow(sigmay151, 2);
            break;
        case 152:
            dJydy[0] = 0.5*(-2*my152 + 2*y152)/pow(sigmay152, 2);
            break;
        case 153:
            dJydy[0] = 0.5*(-2*my153 + 2*y153)/pow(sigmay153, 2);
            break;
        case 154:
            dJydy[0] = 0.5*(-2*my154 + 2*y154)/pow(sigmay154, 2);
            break;
        case 155:
            dJydy[0] = 0.5*(-2*my155 + 2*y155)/pow(sigmay155, 2);
            break;
        case 156:
            dJydy[0] = 0.5*(-2*my156 + 2*y156)/pow(sigmay156, 2);
            break;
        case 157:
            dJydy[0] = 0.5*(-2*my157 + 2*y157)/pow(sigmay157, 2);
            break;
        case 158:
            dJydy[0] = 0.5*(-2*my158 + 2*y158)/pow(sigmay158, 2);
            break;
        case 159:
            dJydy[0] = 0.5*(-2*my159 + 2*y159)/pow(sigmay159, 2);
            break;
        case 160:
            dJydy[0] = 0.5*(-2*my160 + 2*y160)/pow(sigmay160, 2);
            break;
        case 161:
            dJydy[0] = 0.5*(-2*my161 + 2*y161)/pow(sigmay161, 2);
            break;
        case 162:
            dJydy[0] = 0.5*(-2*my162 + 2*y162)/pow(sigmay162, 2);
            break;
        case 163:
            dJydy[0] = 0.5*(-2*my163 + 2*y163)/pow(sigmay163, 2);
            break;
        case 164:
            dJydy[0] = 0.5*(-2*my164 + 2*y164)/pow(sigmay164, 2);
            break;
        case 165:
            dJydy[0] = 0.5*(-2*my165 + 2*y165)/pow(sigmay165, 2);
            break;
        case 166:
            dJydy[0] = 0.5*(-2*my166 + 2*y166)/pow(sigmay166, 2);
            break;
        case 167:
            dJydy[0] = 0.5*(-2*my167 + 2*y167)/pow(sigmay167, 2);
            break;
        case 168:
            dJydy[0] = 0.5*(-2*my168 + 2*y168)/pow(sigmay168, 2);
            break;
        case 169:
            dJydy[0] = 0.5*(-2*my169 + 2*y169)/pow(sigmay169, 2);
            break;
        case 170:
            dJydy[0] = 0.5*(-2*my170 + 2*y170)/pow(sigmay170, 2);
            break;
        case 171:
            dJydy[0] = 0.5*(-2*my171 + 2*y171)/pow(sigmay171, 2);
            break;
        case 172:
            dJydy[0] = 0.5*(-2*my172 + 2*y172)/pow(sigmay172, 2);
            break;
        case 173:
            dJydy[0] = 0.5*(-2*my173 + 2*y173)/pow(sigmay173, 2);
            break;
        case 174:
            dJydy[0] = 0.5*(-2*my174 + 2*y174)/pow(sigmay174, 2);
            break;
        case 175:
            dJydy[0] = 0.5*(-2*my175 + 2*y175)/pow(sigmay175, 2);
            break;
        case 176:
            dJydy[0] = 0.5*(-2*my176 + 2*y176)/pow(sigmay176, 2);
            break;
        case 177:
            dJydy[0] = 0.5*(-2*my177 + 2*y177)/pow(sigmay177, 2);
            break;
        case 178:
            dJydy[0] = 0.5*(-2*my178 + 2*y178)/pow(sigmay178, 2);
            break;
        case 179:
            dJydy[0] = 0.5*(-2*my179 + 2*y179)/pow(sigmay179, 2);
            break;
        case 180:
            dJydy[0] = 0.5*(-2*my180 + 2*y180)/pow(sigmay180, 2);
            break;
        case 181:
            dJydy[0] = 0.5*(-2*my181 + 2*y181)/pow(sigmay181, 2);
            break;
        case 182:
            dJydy[0] = 0.5*(-2*my182 + 2*y182)/pow(sigmay182, 2);
            break;
        case 183:
            dJydy[0] = 0.5*(-2*my183 + 2*y183)/pow(sigmay183, 2);
            break;
        case 184:
            dJydy[0] = 0.5*(-2*my184 + 2*y184)/pow(sigmay184, 2);
            break;
        case 185:
            dJydy[0] = 0.5*(-2*my185 + 2*y185)/pow(sigmay185, 2);
            break;
        case 186:
            dJydy[0] = 0.5*(-2*my186 + 2*y186)/pow(sigmay186, 2);
            break;
        case 187:
            dJydy[0] = 0.5*(-2*my187 + 2*y187)/pow(sigmay187, 2);
            break;
        case 188:
            dJydy[0] = 0.5*(-2*my188 + 2*y188)/pow(sigmay188, 2);
            break;
        case 189:
            dJydy[0] = 0.5*(-2*my189 + 2*y189)/pow(sigmay189, 2);
            break;
        case 190:
            dJydy[0] = 0.5*(-2*my190 + 2*y190)/pow(sigmay190, 2);
            break;
        case 191:
            dJydy[0] = 0.5*(-2*my191 + 2*y191)/pow(sigmay191, 2);
            break;
        case 192:
            dJydy[0] = 0.5*(-2*my192 + 2*y192)/pow(sigmay192, 2);
            break;
        case 193:
            dJydy[0] = 0.5*(-2*my193 + 2*y193)/pow(sigmay193, 2);
            break;
        case 194:
            dJydy[0] = 0.5*(-2*my194 + 2*y194)/pow(sigmay194, 2);
            break;
        case 195:
            dJydy[0] = 0.5*(-2*my195 + 2*y195)/pow(sigmay195, 2);
            break;
        case 196:
            dJydy[0] = 0.5*(-2*my196 + 2*y196)/pow(sigmay196, 2);
            break;
        case 197:
            dJydy[0] = 0.5*(-2*my197 + 2*y197)/pow(sigmay197, 2);
            break;
        case 198:
            dJydy[0] = 0.5*(-2*my198 + 2*y198)/pow(sigmay198, 2);
            break;
        case 199:
            dJydy[0] = 0.5*(-2*my199 + 2*y199)/pow(sigmay199, 2);
            break;
        case 200:
            dJydy[0] = 0.5*(-2*my200 + 2*y200)/pow(sigmay200, 2);
            break;
        case 201:
            dJydy[0] = 0.5*(-2*my201 + 2*y201)/pow(sigmay201, 2);
            break;
        case 202:
            dJydy[0] = 0.5*(-2*my202 + 2*y202)/pow(sigmay202, 2);
            break;
        case 203:
            dJydy[0] = 0.5*(-2*my203 + 2*y203)/pow(sigmay203, 2);
            break;
        case 204:
            dJydy[0] = 0.5*(-2*my204 + 2*y204)/pow(sigmay204, 2);
            break;
        case 205:
            dJydy[0] = 0.5*(-2*my205 + 2*y205)/pow(sigmay205, 2);
            break;
        case 206:
            dJydy[0] = 0.5*(-2*my206 + 2*y206)/pow(sigmay206, 2);
            break;
        case 207:
            dJydy[0] = 0.5*(-2*my207 + 2*y207)/pow(sigmay207, 2);
            break;
        case 208:
            dJydy[0] = 0.5*(-2*my208 + 2*y208)/pow(sigmay208, 2);
            break;
        case 209:
            dJydy[0] = 0.5*(-2*my209 + 2*y209)/pow(sigmay209, 2);
            break;
        case 210:
            dJydy[0] = 0.5*(-2*my210 + 2*y210)/pow(sigmay210, 2);
            break;
        case 211:
            dJydy[0] = 0.5*(-2*my211 + 2*y211)/pow(sigmay211, 2);
            break;
        case 212:
            dJydy[0] = 0.5*(-2*my212 + 2*y212)/pow(sigmay212, 2);
            break;
        case 213:
            dJydy[0] = 0.5*(-2*my213 + 2*y213)/pow(sigmay213, 2);
            break;
        case 214:
            dJydy[0] = 0.5*(-2*my214 + 2*y214)/pow(sigmay214, 2);
            break;
        case 215:
            dJydy[0] = 0.5*(-2*my215 + 2*y215)/pow(sigmay215, 2);
            break;
        case 216:
            dJydy[0] = 0.5*(-2*my216 + 2*y216)/pow(sigmay216, 2);
            break;
        case 217:
            dJydy[0] = 0.5*(-2*my217 + 2*y217)/pow(sigmay217, 2);
            break;
        case 218:
            dJydy[0] = 0.5*(-2*my218 + 2*y218)/pow(sigmay218, 2);
            break;
        case 219:
            dJydy[0] = 0.5*(-2*my219 + 2*y219)/pow(sigmay219, 2);
            break;
        case 220:
            dJydy[0] = 0.5*(-2*my220 + 2*y220)/pow(sigmay220, 2);
            break;
        case 221:
            dJydy[0] = 0.5*(-2*my221 + 2*y221)/pow(sigmay221, 2);
            break;
        case 222:
            dJydy[0] = 0.5*(-2*my222 + 2*y222)/pow(sigmay222, 2);
            break;
        case 223:
            dJydy[0] = 0.5*(-2*my223 + 2*y223)/pow(sigmay223, 2);
            break;
        case 224:
            dJydy[0] = 0.5*(-2*my224 + 2*y224)/pow(sigmay224, 2);
            break;
        case 225:
            dJydy[0] = 0.5*(-2*my225 + 2*y225)/pow(sigmay225, 2);
            break;
        case 226:
            dJydy[0] = 0.5*(-2*my226 + 2*y226)/pow(sigmay226, 2);
            break;
        case 227:
            dJydy[0] = 0.5*(-2*my227 + 2*y227)/pow(sigmay227, 2);
            break;
        case 228:
            dJydy[0] = 0.5*(-2*my228 + 2*y228)/pow(sigmay228, 2);
            break;
        case 229:
            dJydy[0] = 0.5*(-2*my229 + 2*y229)/pow(sigmay229, 2);
            break;
        case 230:
            dJydy[0] = 0.5*(-2*my230 + 2*y230)/pow(sigmay230, 2);
            break;
        case 231:
            dJydy[0] = 0.5*(-2*my231 + 2*y231)/pow(sigmay231, 2);
            break;
        case 232:
            dJydy[0] = 0.5*(-2*my232 + 2*y232)/pow(sigmay232, 2);
            break;
        case 233:
            dJydy[0] = 0.5*(-2*my233 + 2*y233)/pow(sigmay233, 2);
            break;
        case 234:
            dJydy[0] = 0.5*(-2*my234 + 2*y234)/pow(sigmay234, 2);
            break;
        case 235:
            dJydy[0] = 0.5*(-2*my235 + 2*y235)/pow(sigmay235, 2);
            break;
        case 236:
            dJydy[0] = 0.5*(-2*my236 + 2*y236)/pow(sigmay236, 2);
            break;
        case 237:
            dJydy[0] = 0.5*(-2*my237 + 2*y237)/pow(sigmay237, 2);
            break;
        case 238:
            dJydy[0] = 0.5*(-2*my238 + 2*y238)/pow(sigmay238, 2);
            break;
        case 239:
            dJydy[0] = 0.5*(-2*my239 + 2*y239)/pow(sigmay239, 2);
            break;
        case 240:
            dJydy[0] = 0.5*(-2*my240 + 2*y240)/pow(sigmay240, 2);
            break;
        case 241:
            dJydy[0] = 0.5*(-2*my241 + 2*y241)/pow(sigmay241, 2);
            break;
        case 242:
            dJydy[0] = 0.5*(-2*my242 + 2*y242)/pow(sigmay242, 2);
            break;
        case 243:
            dJydy[0] = 0.5*(-2*my243 + 2*y243)/pow(sigmay243, 2);
            break;
        case 244:
            dJydy[0] = 0.5*(-2*my244 + 2*y244)/pow(sigmay244, 2);
            break;
        case 245:
            dJydy[0] = 0.5*(-2*my245 + 2*y245)/pow(sigmay245, 2);
            break;
        case 246:
            dJydy[0] = 0.5*(-2*my246 + 2*y246)/pow(sigmay246, 2);
            break;
        case 247:
            dJydy[0] = 0.5*(-2*my247 + 2*y247)/pow(sigmay247, 2);
            break;
        case 248:
            dJydy[0] = 0.5*(-2*my248 + 2*y248)/pow(sigmay248, 2);
            break;
        case 249:
            dJydy[0] = 0.5*(-2*my249 + 2*y249)/pow(sigmay249, 2);
            break;
        case 250:
            dJydy[0] = 0.5*(-2*my250 + 2*y250)/pow(sigmay250, 2);
            break;
        case 251:
            dJydy[0] = 0.5*(-2*my251 + 2*y251)/pow(sigmay251, 2);
            break;
        case 252:
            dJydy[0] = 0.5*(-2*my252 + 2*y252)/pow(sigmay252, 2);
            break;
        case 253:
            dJydy[0] = 0.5*(-2*my253 + 2*y253)/pow(sigmay253, 2);
            break;
        case 254:
            dJydy[0] = 0.5*(-2*my254 + 2*y254)/pow(sigmay254, 2);
            break;
        case 255:
            dJydy[0] = 0.5*(-2*my255 + 2*y255)/pow(sigmay255, 2);
            break;
        case 256:
            dJydy[0] = 0.5*(-2*my256 + 2*y256)/pow(sigmay256, 2);
            break;
        case 257:
            dJydy[0] = 0.5*(-2*my257 + 2*y257)/pow(sigmay257, 2);
            break;
        case 258:
            dJydy[0] = 0.5*(-2*my258 + 2*y258)/pow(sigmay258, 2);
            break;
        case 259:
            dJydy[0] = 0.5*(-2*my259 + 2*y259)/pow(sigmay259, 2);
            break;
        case 260:
            dJydy[0] = 0.5*(-2*my260 + 2*y260)/pow(sigmay260, 2);
            break;
        case 261:
            dJydy[0] = 0.5*(-2*my261 + 2*y261)/pow(sigmay261, 2);
            break;
        case 262:
            dJydy[0] = 0.5*(-2*my262 + 2*y262)/pow(sigmay262, 2);
            break;
        case 263:
            dJydy[0] = 0.5*(-2*my263 + 2*y263)/pow(sigmay263, 2);
            break;
        case 264:
            dJydy[0] = 0.5*(-2*my264 + 2*y264)/pow(sigmay264, 2);
            break;
        case 265:
            dJydy[0] = 0.5*(-2*my265 + 2*y265)/pow(sigmay265, 2);
            break;
        case 266:
            dJydy[0] = 0.5*(-2*my266 + 2*y266)/pow(sigmay266, 2);
            break;
        case 267:
            dJydy[0] = 0.5*(-2*my267 + 2*y267)/pow(sigmay267, 2);
            break;
        case 268:
            dJydy[0] = 0.5*(-2*my268 + 2*y268)/pow(sigmay268, 2);
            break;
        case 269:
            dJydy[0] = 0.5*(-2*my269 + 2*y269)/pow(sigmay269, 2);
            break;
        case 270:
            dJydy[0] = 0.5*(-2*my270 + 2*y270)/pow(sigmay270, 2);
            break;
        case 271:
            dJydy[0] = 0.5*(-2*my271 + 2*y271)/pow(sigmay271, 2);
            break;
        case 272:
            dJydy[0] = 0.5*(-2*my272 + 2*y272)/pow(sigmay272, 2);
            break;
        case 273:
            dJydy[0] = 0.5*(-2*my273 + 2*y273)/pow(sigmay273, 2);
            break;
        case 274:
            dJydy[0] = 0.5*(-2*my274 + 2*y274)/pow(sigmay274, 2);
            break;
        case 275:
            dJydy[0] = 0.5*(-2*my275 + 2*y275)/pow(sigmay275, 2);
            break;
        case 276:
            dJydy[0] = 0.5*(-2*my276 + 2*y276)/pow(sigmay276, 2);
            break;
        case 277:
            dJydy[0] = 0.5*(-2*my277 + 2*y277)/pow(sigmay277, 2);
            break;
        case 278:
            dJydy[0] = 0.5*(-2*my278 + 2*y278)/pow(sigmay278, 2);
            break;
        case 279:
            dJydy[0] = 0.5*(-2*my279 + 2*y279)/pow(sigmay279, 2);
            break;
        case 280:
            dJydy[0] = 0.5*(-2*my280 + 2*y280)/pow(sigmay280, 2);
            break;
        case 281:
            dJydy[0] = 0.5*(-2*my281 + 2*y281)/pow(sigmay281, 2);
            break;
        case 282:
            dJydy[0] = 0.5*(-2*my282 + 2*y282)/pow(sigmay282, 2);
            break;
        case 283:
            dJydy[0] = 0.5*(-2*my283 + 2*y283)/pow(sigmay283, 2);
            break;
        case 284:
            dJydy[0] = 0.5*(-2*my284 + 2*y284)/pow(sigmay284, 2);
            break;
        case 285:
            dJydy[0] = 0.5*(-2*my285 + 2*y285)/pow(sigmay285, 2);
            break;
        case 286:
            dJydy[0] = 0.5*(-2*my286 + 2*y286)/pow(sigmay286, 2);
            break;
        case 287:
            dJydy[0] = 0.5*(-2*my287 + 2*y287)/pow(sigmay287, 2);
            break;
        case 288:
            dJydy[0] = 0.5*(-2*my288 + 2*y288)/pow(sigmay288, 2);
            break;
        case 289:
            dJydy[0] = 0.5*(-2*my289 + 2*y289)/pow(sigmay289, 2);
            break;
        case 290:
            dJydy[0] = 0.5*(-2*my290 + 2*y290)/pow(sigmay290, 2);
            break;
        case 291:
            dJydy[0] = 0.5*(-2*my291 + 2*y291)/pow(sigmay291, 2);
            break;
        case 292:
            dJydy[0] = 0.5*(-2*my292 + 2*y292)/pow(sigmay292, 2);
            break;
        case 293:
            dJydy[0] = 0.5*(-2*my293 + 2*y293)/pow(sigmay293, 2);
            break;
        case 294:
            dJydy[0] = 0.5*(-2*my294 + 2*y294)/pow(sigmay294, 2);
            break;
        case 295:
            dJydy[0] = 0.5*(-2*my295 + 2*y295)/pow(sigmay295, 2);
            break;
        case 296:
            dJydy[0] = 0.5*(-2*my296 + 2*y296)/pow(sigmay296, 2);
            break;
        case 297:
            dJydy[0] = 0.5*(-2*my297 + 2*y297)/pow(sigmay297, 2);
            break;
        case 298:
            dJydy[0] = 0.5*(-2*my298 + 2*y298)/pow(sigmay298, 2);
            break;
        case 299:
            dJydy[0] = 0.5*(-2*my299 + 2*y299)/pow(sigmay299, 2);
            break;
        case 300:
            dJydy[0] = 0.5*(-2*my300 + 2*y300)/pow(sigmay300, 2);
            break;
        case 301:
            dJydy[0] = 0.5*(-2*my301 + 2*y301)/pow(sigmay301, 2);
            break;
        case 302:
            dJydy[0] = 0.5*(-2*my302 + 2*y302)/pow(sigmay302, 2);
            break;
        case 303:
            dJydy[0] = 0.5*(-2*my303 + 2*y303)/pow(sigmay303, 2);
            break;
        case 304:
            dJydy[0] = 0.5*(-2*my304 + 2*y304)/pow(sigmay304, 2);
            break;
        case 305:
            dJydy[0] = 0.5*(-2*my305 + 2*y305)/pow(sigmay305, 2);
            break;
        case 306:
            dJydy[0] = 0.5*(-2*my306 + 2*y306)/pow(sigmay306, 2);
            break;
        case 307:
            dJydy[0] = 0.5*(-2*my307 + 2*y307)/pow(sigmay307, 2);
            break;
        case 308:
            dJydy[0] = 0.5*(-2*my308 + 2*y308)/pow(sigmay308, 2);
            break;
        case 309:
            dJydy[0] = 0.5*(-2*my309 + 2*y309)/pow(sigmay309, 2);
            break;
        case 310:
            dJydy[0] = 0.5*(-2*my310 + 2*y310)/pow(sigmay310, 2);
            break;
        case 311:
            dJydy[0] = 0.5*(-2*my311 + 2*y311)/pow(sigmay311, 2);
            break;
        case 312:
            dJydy[0] = 0.5*(-2*my312 + 2*y312)/pow(sigmay312, 2);
            break;
        case 313:
            dJydy[0] = 0.5*(-2*my313 + 2*y313)/pow(sigmay313, 2);
            break;
        case 314:
            dJydy[0] = 0.5*(-2*my314 + 2*y314)/pow(sigmay314, 2);
            break;
        case 315:
            dJydy[0] = 0.5*(-2*my315 + 2*y315)/pow(sigmay315, 2);
            break;
        case 316:
            dJydy[0] = 0.5*(-2*my316 + 2*y316)/pow(sigmay316, 2);
            break;
        case 317:
            dJydy[0] = 0.5*(-2*my317 + 2*y317)/pow(sigmay317, 2);
            break;
        case 318:
            dJydy[0] = 0.5*(-2*my318 + 2*y318)/pow(sigmay318, 2);
            break;
        case 319:
            dJydy[0] = 0.5*(-2*my319 + 2*y319)/pow(sigmay319, 2);
            break;
        case 320:
            dJydy[0] = 0.5*(-2*my320 + 2*y320)/pow(sigmay320, 2);
            break;
        case 321:
            dJydy[0] = 0.5*(-2*my321 + 2*y321)/pow(sigmay321, 2);
            break;
        case 322:
            dJydy[0] = 0.5*(-2*my322 + 2*y322)/pow(sigmay322, 2);
            break;
        case 323:
            dJydy[0] = 0.5*(-2*my323 + 2*y323)/pow(sigmay323, 2);
            break;
        case 324:
            dJydy[0] = 0.5*(-2*my324 + 2*y324)/pow(sigmay324, 2);
            break;
        case 325:
            dJydy[0] = 0.5*(-2*my325 + 2*y325)/pow(sigmay325, 2);
            break;
        case 326:
            dJydy[0] = 0.5*(-2*my326 + 2*y326)/pow(sigmay326, 2);
            break;
        case 327:
            dJydy[0] = 0.5*(-2*my327 + 2*y327)/pow(sigmay327, 2);
            break;
        case 328:
            dJydy[0] = 0.5*(-2*my328 + 2*y328)/pow(sigmay328, 2);
            break;
        case 329:
            dJydy[0] = 0.5*(-2*my329 + 2*y329)/pow(sigmay329, 2);
            break;
        case 330:
            dJydy[0] = 0.5*(-2*my330 + 2*y330)/pow(sigmay330, 2);
            break;
        case 331:
            dJydy[0] = 0.5*(-2*my331 + 2*y331)/pow(sigmay331, 2);
            break;
        case 332:
            dJydy[0] = 0.5*(-2*my332 + 2*y332)/pow(sigmay332, 2);
            break;
        case 333:
            dJydy[0] = 0.5*(-2*my333 + 2*y333)/pow(sigmay333, 2);
            break;
        case 334:
            dJydy[0] = 0.5*(-2*my334 + 2*y334)/pow(sigmay334, 2);
            break;
        case 335:
            dJydy[0] = 0.5*(-2*my335 + 2*y335)/pow(sigmay335, 2);
            break;
        case 336:
            dJydy[0] = 0.5*(-2*my336 + 2*y336)/pow(sigmay336, 2);
            break;
        case 337:
            dJydy[0] = 0.5*(-2*my337 + 2*y337)/pow(sigmay337, 2);
            break;
        case 338:
            dJydy[0] = 0.5*(-2*my338 + 2*y338)/pow(sigmay338, 2);
            break;
        case 339:
            dJydy[0] = 0.5*(-2*my339 + 2*y339)/pow(sigmay339, 2);
            break;
        case 340:
            dJydy[0] = 0.5*(-2*my340 + 2*y340)/pow(sigmay340, 2);
            break;
        case 341:
            dJydy[0] = 0.5*(-2*my341 + 2*y341)/pow(sigmay341, 2);
            break;
        case 342:
            dJydy[0] = 0.5*(-2*my342 + 2*y342)/pow(sigmay342, 2);
            break;
        case 343:
            dJydy[0] = 0.5*(-2*my343 + 2*y343)/pow(sigmay343, 2);
            break;
        case 344:
            dJydy[0] = 0.5*(-2*my344 + 2*y344)/pow(sigmay344, 2);
            break;
        case 345:
            dJydy[0] = 0.5*(-2*my345 + 2*y345)/pow(sigmay345, 2);
            break;
        case 346:
            dJydy[0] = 0.5*(-2*my346 + 2*y346)/pow(sigmay346, 2);
            break;
        case 347:
            dJydy[0] = 0.5*(-2*my347 + 2*y347)/pow(sigmay347, 2);
            break;
        case 348:
            dJydy[0] = 0.5*(-2*my348 + 2*y348)/pow(sigmay348, 2);
            break;
        case 349:
            dJydy[0] = 0.5*(-2*my349 + 2*y349)/pow(sigmay349, 2);
            break;
        case 350:
            dJydy[0] = 0.5*(-2*my350 + 2*y350)/pow(sigmay350, 2);
            break;
        case 351:
            dJydy[0] = 0.5*(-2*my351 + 2*y351)/pow(sigmay351, 2);
            break;
        case 352:
            dJydy[0] = 0.5*(-2*my352 + 2*y352)/pow(sigmay352, 2);
            break;
        case 353:
            dJydy[0] = 0.5*(-2*my353 + 2*y353)/pow(sigmay353, 2);
            break;
        case 354:
            dJydy[0] = 0.5*(-2*my354 + 2*y354)/pow(sigmay354, 2);
            break;
        case 355:
            dJydy[0] = 0.5*(-2*my355 + 2*y355)/pow(sigmay355, 2);
            break;
        case 356:
            dJydy[0] = 0.5*(-2*my356 + 2*y356)/pow(sigmay356, 2);
            break;
        case 357:
            dJydy[0] = 0.5*(-2*my357 + 2*y357)/pow(sigmay357, 2);
            break;
        case 358:
            dJydy[0] = 0.5*(-2*my358 + 2*y358)/pow(sigmay358, 2);
            break;
        case 359:
            dJydy[0] = 0.5*(-2*my359 + 2*y359)/pow(sigmay359, 2);
            break;
        case 360:
            dJydy[0] = 0.5*(-2*my360 + 2*y360)/pow(sigmay360, 2);
            break;
        case 361:
            dJydy[0] = 0.5*(-2*my361 + 2*y361)/pow(sigmay361, 2);
            break;
        case 362:
            dJydy[0] = 0.5*(-2*my362 + 2*y362)/pow(sigmay362, 2);
            break;
        case 363:
            dJydy[0] = 0.5*(-2*my363 + 2*y363)/pow(sigmay363, 2);
            break;
        case 364:
            dJydy[0] = 0.5*(-2*my364 + 2*y364)/pow(sigmay364, 2);
            break;
        case 365:
            dJydy[0] = 0.5*(-2*my365 + 2*y365)/pow(sigmay365, 2);
            break;
        case 366:
            dJydy[0] = 0.5*(-2*my366 + 2*y366)/pow(sigmay366, 2);
            break;
        case 367:
            dJydy[0] = 0.5*(-2*my367 + 2*y367)/pow(sigmay367, 2);
            break;
        case 368:
            dJydy[0] = 0.5*(-2*my368 + 2*y368)/pow(sigmay368, 2);
            break;
        case 369:
            dJydy[0] = 0.5*(-2*my369 + 2*y369)/pow(sigmay369, 2);
            break;
        case 370:
            dJydy[0] = 0.5*(-2*my370 + 2*y370)/pow(sigmay370, 2);
            break;
        case 371:
            dJydy[0] = 0.5*(-2*my371 + 2*y371)/pow(sigmay371, 2);
            break;
        case 372:
            dJydy[0] = 0.5*(-2*my372 + 2*y372)/pow(sigmay372, 2);
            break;
        case 373:
            dJydy[0] = 0.5*(-2*my373 + 2*y373)/pow(sigmay373, 2);
            break;
        case 374:
            dJydy[0] = 0.5*(-2*my374 + 2*y374)/pow(sigmay374, 2);
            break;
        case 375:
            dJydy[0] = 0.5*(-2*my375 + 2*y375)/pow(sigmay375, 2);
            break;
        case 376:
            dJydy[0] = 0.5*(-2*my376 + 2*y376)/pow(sigmay376, 2);
            break;
        case 377:
            dJydy[0] = 0.5*(-2*my377 + 2*y377)/pow(sigmay377, 2);
            break;
        case 378:
            dJydy[0] = 0.5*(-2*my378 + 2*y378)/pow(sigmay378, 2);
            break;
        case 379:
            dJydy[0] = 0.5*(-2*my379 + 2*y379)/pow(sigmay379, 2);
            break;
        case 380:
            dJydy[0] = 0.5*(-2*my380 + 2*y380)/pow(sigmay380, 2);
            break;
        case 381:
            dJydy[0] = 0.5*(-2*my381 + 2*y381)/pow(sigmay381, 2);
            break;
        case 382:
            dJydy[0] = 0.5*(-2*my382 + 2*y382)/pow(sigmay382, 2);
            break;
        case 383:
            dJydy[0] = 0.5*(-2*my383 + 2*y383)/pow(sigmay383, 2);
            break;
        case 384:
            dJydy[0] = 0.5*(-2*my384 + 2*y384)/pow(sigmay384, 2);
            break;
        case 385:
            dJydy[0] = 0.5*(-2*my385 + 2*y385)/pow(sigmay385, 2);
            break;
        case 386:
            dJydy[0] = 0.5*(-2*my386 + 2*y386)/pow(sigmay386, 2);
            break;
        case 387:
            dJydy[0] = 0.5*(-2*my387 + 2*y387)/pow(sigmay387, 2);
            break;
        case 388:
            dJydy[0] = 0.5*(-2*my388 + 2*y388)/pow(sigmay388, 2);
            break;
        case 389:
            dJydy[0] = 0.5*(-2*my389 + 2*y389)/pow(sigmay389, 2);
            break;
        case 390:
            dJydy[0] = 0.5*(-2*my390 + 2*y390)/pow(sigmay390, 2);
            break;
        case 391:
            dJydy[0] = 0.5*(-2*my391 + 2*y391)/pow(sigmay391, 2);
            break;
        case 392:
            dJydy[0] = 0.5*(-2*my392 + 2*y392)/pow(sigmay392, 2);
            break;
        case 393:
            dJydy[0] = 0.5*(-2*my393 + 2*y393)/pow(sigmay393, 2);
            break;
        case 394:
            dJydy[0] = 0.5*(-2*my394 + 2*y394)/pow(sigmay394, 2);
            break;
        case 395:
            dJydy[0] = 0.5*(-2*my395 + 2*y395)/pow(sigmay395, 2);
            break;
        case 396:
            dJydy[0] = 0.5*(-2*my396 + 2*y396)/pow(sigmay396, 2);
            break;
        case 397:
            dJydy[0] = 0.5*(-2*my397 + 2*y397)/pow(sigmay397, 2);
            break;
        case 398:
            dJydy[0] = 0.5*(-2*my398 + 2*y398)/pow(sigmay398, 2);
            break;
        case 399:
            dJydy[0] = 0.5*(-2*my399 + 2*y399)/pow(sigmay399, 2);
            break;
        case 400:
            dJydy[0] = 0.5*(-2*my400 + 2*y400)/pow(sigmay400, 2);
            break;
        case 401:
            dJydy[0] = 0.5*(-2*my401 + 2*y401)/pow(sigmay401, 2);
            break;
        case 402:
            dJydy[0] = 0.5*(-2*my402 + 2*y402)/pow(sigmay402, 2);
            break;
        case 403:
            dJydy[0] = 0.5*(-2*my403 + 2*y403)/pow(sigmay403, 2);
            break;
        case 404:
            dJydy[0] = 0.5*(-2*my404 + 2*y404)/pow(sigmay404, 2);
            break;
        case 405:
            dJydy[0] = 0.5*(-2*my405 + 2*y405)/pow(sigmay405, 2);
            break;
        case 406:
            dJydy[0] = 0.5*(-2*my406 + 2*y406)/pow(sigmay406, 2);
            break;
        case 407:
            dJydy[0] = 0.5*(-2*my407 + 2*y407)/pow(sigmay407, 2);
            break;
        case 408:
            dJydy[0] = 0.5*(-2*my408 + 2*y408)/pow(sigmay408, 2);
            break;
        case 409:
            dJydy[0] = 0.5*(-2*my409 + 2*y409)/pow(sigmay409, 2);
            break;
        case 410:
            dJydy[0] = 0.5*(-2*my410 + 2*y410)/pow(sigmay410, 2);
            break;
        case 411:
            dJydy[0] = 0.5*(-2*my411 + 2*y411)/pow(sigmay411, 2);
            break;
        case 412:
            dJydy[0] = 0.5*(-2*my412 + 2*y412)/pow(sigmay412, 2);
            break;
        case 413:
            dJydy[0] = 0.5*(-2*my413 + 2*y413)/pow(sigmay413, 2);
            break;
        case 414:
            dJydy[0] = 0.5*(-2*my414 + 2*y414)/pow(sigmay414, 2);
            break;
        case 415:
            dJydy[0] = 0.5*(-2*my415 + 2*y415)/pow(sigmay415, 2);
            break;
        case 416:
            dJydy[0] = 0.5*(-2*my416 + 2*y416)/pow(sigmay416, 2);
            break;
        case 417:
            dJydy[0] = 0.5*(-2*my417 + 2*y417)/pow(sigmay417, 2);
            break;
        case 418:
            dJydy[0] = 0.5*(-2*my418 + 2*y418)/pow(sigmay418, 2);
            break;
        case 419:
            dJydy[0] = 0.5*(-2*my419 + 2*y419)/pow(sigmay419, 2);
            break;
        case 420:
            dJydy[0] = 0.5*(-2*my420 + 2*y420)/pow(sigmay420, 2);
            break;
        case 421:
            dJydy[0] = 0.5*(-2*my421 + 2*y421)/pow(sigmay421, 2);
            break;
        case 422:
            dJydy[0] = 0.5*(-2*my422 + 2*y422)/pow(sigmay422, 2);
            break;
        case 423:
            dJydy[0] = 0.5*(-2*my423 + 2*y423)/pow(sigmay423, 2);
            break;
        case 424:
            dJydy[0] = 0.5*(-2*my424 + 2*y424)/pow(sigmay424, 2);
            break;
        case 425:
            dJydy[0] = 0.5*(-2*my425 + 2*y425)/pow(sigmay425, 2);
            break;
        case 426:
            dJydy[0] = 0.5*(-2*my426 + 2*y426)/pow(sigmay426, 2);
            break;
        case 427:
            dJydy[0] = 0.5*(-2*my427 + 2*y427)/pow(sigmay427, 2);
            break;
        case 428:
            dJydy[0] = 0.5*(-2*my428 + 2*y428)/pow(sigmay428, 2);
            break;
        case 429:
            dJydy[0] = 0.5*(-2*my429 + 2*y429)/pow(sigmay429, 2);
            break;
        case 430:
            dJydy[0] = 0.5*(-2*my430 + 2*y430)/pow(sigmay430, 2);
            break;
        case 431:
            dJydy[0] = 0.5*(-2*my431 + 2*y431)/pow(sigmay431, 2);
            break;
        case 432:
            dJydy[0] = 0.5*(-2*my432 + 2*y432)/pow(sigmay432, 2);
            break;
        case 433:
            dJydy[0] = 0.5*(-2*my433 + 2*y433)/pow(sigmay433, 2);
            break;
        case 434:
            dJydy[0] = 0.5*(-2*my434 + 2*y434)/pow(sigmay434, 2);
            break;
        case 435:
            dJydy[0] = 0.5*(-2*my435 + 2*y435)/pow(sigmay435, 2);
            break;
        case 436:
            dJydy[0] = 0.5*(-2*my436 + 2*y436)/pow(sigmay436, 2);
            break;
        case 437:
            dJydy[0] = 0.5*(-2*my437 + 2*y437)/pow(sigmay437, 2);
            break;
        case 438:
            dJydy[0] = 0.5*(-2*my438 + 2*y438)/pow(sigmay438, 2);
            break;
        case 439:
            dJydy[0] = 0.5*(-2*my439 + 2*y439)/pow(sigmay439, 2);
            break;
        case 440:
            dJydy[0] = 0.5*(-2*my440 + 2*y440)/pow(sigmay440, 2);
            break;
        case 441:
            dJydy[0] = 0.5*(-2*my441 + 2*y441)/pow(sigmay441, 2);
            break;
        case 442:
            dJydy[0] = 0.5*(-2*my442 + 2*y442)/pow(sigmay442, 2);
            break;
        case 443:
            dJydy[0] = 0.5*(-2*my443 + 2*y443)/pow(sigmay443, 2);
            break;
        case 444:
            dJydy[0] = 0.5*(-2*my444 + 2*y444)/pow(sigmay444, 2);
            break;
        case 445:
            dJydy[0] = 0.5*(-2*my445 + 2*y445)/pow(sigmay445, 2);
            break;
        case 446:
            dJydy[0] = 0.5*(-2*my446 + 2*y446)/pow(sigmay446, 2);
            break;
        case 447:
            dJydy[0] = 0.5*(-2*my447 + 2*y447)/pow(sigmay447, 2);
            break;
        case 448:
            dJydy[0] = 0.5*(-2*my448 + 2*y448)/pow(sigmay448, 2);
            break;
        case 449:
            dJydy[0] = 0.5*(-2*my449 + 2*y449)/pow(sigmay449, 2);
            break;
        case 450:
            dJydy[0] = 0.5*(-2*my450 + 2*y450)/pow(sigmay450, 2);
            break;
        case 451:
            dJydy[0] = 0.5*(-2*my451 + 2*y451)/pow(sigmay451, 2);
            break;
        case 452:
            dJydy[0] = 0.5*(-2*my452 + 2*y452)/pow(sigmay452, 2);
            break;
        case 453:
            dJydy[0] = 0.5*(-2*my453 + 2*y453)/pow(sigmay453, 2);
            break;
        case 454:
            dJydy[0] = 0.5*(-2*my454 + 2*y454)/pow(sigmay454, 2);
            break;
        case 455:
            dJydy[0] = 0.5*(-2*my455 + 2*y455)/pow(sigmay455, 2);
            break;
        case 456:
            dJydy[0] = 0.5*(-2*my456 + 2*y456)/pow(sigmay456, 2);
            break;
        case 457:
            dJydy[0] = 0.5*(-2*my457 + 2*y457)/pow(sigmay457, 2);
            break;
        case 458:
            dJydy[0] = 0.5*(-2*my458 + 2*y458)/pow(sigmay458, 2);
            break;
        case 459:
            dJydy[0] = 0.5*(-2*my459 + 2*y459)/pow(sigmay459, 2);
            break;
        case 460:
            dJydy[0] = 0.5*(-2*my460 + 2*y460)/pow(sigmay460, 2);
            break;
        case 461:
            dJydy[0] = 0.5*(-2*my461 + 2*y461)/pow(sigmay461, 2);
            break;
        case 462:
            dJydy[0] = 0.5*(-2*my462 + 2*y462)/pow(sigmay462, 2);
            break;
        case 463:
            dJydy[0] = 0.5*(-2*my463 + 2*y463)/pow(sigmay463, 2);
            break;
        case 464:
            dJydy[0] = 0.5*(-2*my464 + 2*y464)/pow(sigmay464, 2);
            break;
        case 465:
            dJydy[0] = 0.5*(-2*my465 + 2*y465)/pow(sigmay465, 2);
            break;
        case 466:
            dJydy[0] = 0.5*(-2*my466 + 2*y466)/pow(sigmay466, 2);
            break;
        case 467:
            dJydy[0] = 0.5*(-2*my467 + 2*y467)/pow(sigmay467, 2);
            break;
        case 468:
            dJydy[0] = 0.5*(-2*my468 + 2*y468)/pow(sigmay468, 2);
            break;
        case 469:
            dJydy[0] = 0.5*(-2*my469 + 2*y469)/pow(sigmay469, 2);
            break;
        case 470:
            dJydy[0] = 0.5*(-2*my470 + 2*y470)/pow(sigmay470, 2);
            break;
        case 471:
            dJydy[0] = 0.5*(-2*my471 + 2*y471)/pow(sigmay471, 2);
            break;
        case 472:
            dJydy[0] = 0.5*(-2*my472 + 2*y472)/pow(sigmay472, 2);
            break;
        case 473:
            dJydy[0] = 0.5*(-2*my473 + 2*y473)/pow(sigmay473, 2);
            break;
        case 474:
            dJydy[0] = 0.5*(-2*my474 + 2*y474)/pow(sigmay474, 2);
            break;
        case 475:
            dJydy[0] = 0.5*(-2*my475 + 2*y475)/pow(sigmay475, 2);
            break;
        case 476:
            dJydy[0] = 0.5*(-2*my476 + 2*y476)/pow(sigmay476, 2);
            break;
        case 477:
            dJydy[0] = 0.5*(-2*my477 + 2*y477)/pow(sigmay477, 2);
            break;
        case 478:
            dJydy[0] = 0.5*(-2*my478 + 2*y478)/pow(sigmay478, 2);
            break;
        case 479:
            dJydy[0] = 0.5*(-2*my479 + 2*y479)/pow(sigmay479, 2);
            break;
        case 480:
            dJydy[0] = 0.5*(-2*my480 + 2*y480)/pow(sigmay480, 2);
            break;
        case 481:
            dJydy[0] = 0.5*(-2*my481 + 2*y481)/pow(sigmay481, 2);
            break;
        case 482:
            dJydy[0] = 0.5*(-2*my482 + 2*y482)/pow(sigmay482, 2);
            break;
        case 483:
            dJydy[0] = 0.5*(-2*my483 + 2*y483)/pow(sigmay483, 2);
            break;
        case 484:
            dJydy[0] = 0.5*(-2*my484 + 2*y484)/pow(sigmay484, 2);
            break;
        case 485:
            dJydy[0] = 0.5*(-2*my485 + 2*y485)/pow(sigmay485, 2);
            break;
        case 486:
            dJydy[0] = 0.5*(-2*my486 + 2*y486)/pow(sigmay486, 2);
            break;
        case 487:
            dJydy[0] = 0.5*(-2*my487 + 2*y487)/pow(sigmay487, 2);
            break;
        case 488:
            dJydy[0] = 0.5*(-2*my488 + 2*y488)/pow(sigmay488, 2);
            break;
        case 489:
            dJydy[0] = 0.5*(-2*my489 + 2*y489)/pow(sigmay489, 2);
            break;
        case 490:
            dJydy[0] = 0.5*(-2*my490 + 2*y490)/pow(sigmay490, 2);
            break;
        case 491:
            dJydy[0] = 0.5*(-2*my491 + 2*y491)/pow(sigmay491, 2);
            break;
        case 492:
            dJydy[0] = 0.5*(-2*my492 + 2*y492)/pow(sigmay492, 2);
            break;
        case 493:
            dJydy[0] = 0.5*(-2*my493 + 2*y493)/pow(sigmay493, 2);
            break;
        case 494:
            dJydy[0] = 0.5*(-2*my494 + 2*y494)/pow(sigmay494, 2);
            break;
        case 495:
            dJydy[0] = 0.5*(-2*my495 + 2*y495)/pow(sigmay495, 2);
            break;
        case 496:
            dJydy[0] = 0.5*(-2*my496 + 2*y496)/pow(sigmay496, 2);
            break;
        case 497:
            dJydy[0] = 0.5*(-2*my497 + 2*y497)/pow(sigmay497, 2);
            break;
        case 498:
            dJydy[0] = 0.5*(-2*my498 + 2*y498)/pow(sigmay498, 2);
            break;
        case 499:
            dJydy[0] = 0.5*(-2*my499 + 2*y499)/pow(sigmay499, 2);
            break;
        case 500:
            dJydy[0] = 0.5*(-2*my500 + 2*y500)/pow(sigmay500, 2);
            break;
        case 501:
            dJydy[0] = 0.5*(-2*my501 + 2*y501)/pow(sigmay501, 2);
            break;
        case 502:
            dJydy[0] = 0.5*(-2*my502 + 2*y502)/pow(sigmay502, 2);
            break;
        case 503:
            dJydy[0] = 0.5*(-2*my503 + 2*y503)/pow(sigmay503, 2);
            break;
        case 504:
            dJydy[0] = 0.5*(-2*my504 + 2*y504)/pow(sigmay504, 2);
            break;
        case 505:
            dJydy[0] = 0.5*(-2*my505 + 2*y505)/pow(sigmay505, 2);
            break;
        case 506:
            dJydy[0] = 0.5*(-2*my506 + 2*y506)/pow(sigmay506, 2);
            break;
        case 507:
            dJydy[0] = 0.5*(-2*my507 + 2*y507)/pow(sigmay507, 2);
            break;
        case 508:
            dJydy[0] = 0.5*(-2*my508 + 2*y508)/pow(sigmay508, 2);
            break;
        case 509:
            dJydy[0] = 0.5*(-2*my509 + 2*y509)/pow(sigmay509, 2);
            break;
        case 510:
            dJydy[0] = 0.5*(-2*my510 + 2*y510)/pow(sigmay510, 2);
            break;
        case 511:
            dJydy[0] = 0.5*(-2*my511 + 2*y511)/pow(sigmay511, 2);
            break;
        case 512:
            dJydy[0] = 0.5*(-2*my512 + 2*y512)/pow(sigmay512, 2);
            break;
        case 513:
            dJydy[0] = 0.5*(-2*my513 + 2*y513)/pow(sigmay513, 2);
            break;
        case 514:
            dJydy[0] = 0.5*(-2*my514 + 2*y514)/pow(sigmay514, 2);
            break;
        case 515:
            dJydy[0] = 0.5*(-2*my515 + 2*y515)/pow(sigmay515, 2);
            break;
        case 516:
            dJydy[0] = 0.5*(-2*my516 + 2*y516)/pow(sigmay516, 2);
            break;
        case 517:
            dJydy[0] = 0.5*(-2*my517 + 2*y517)/pow(sigmay517, 2);
            break;
        case 518:
            dJydy[0] = 0.5*(-2*my518 + 2*y518)/pow(sigmay518, 2);
            break;
        case 519:
            dJydy[0] = 0.5*(-2*my519 + 2*y519)/pow(sigmay519, 2);
            break;
        case 520:
            dJydy[0] = 0.5*(-2*my520 + 2*y520)/pow(sigmay520, 2);
            break;
        case 521:
            dJydy[0] = 0.5*(-2*my521 + 2*y521)/pow(sigmay521, 2);
            break;
        case 522:
            dJydy[0] = 0.5*(-2*my522 + 2*y522)/pow(sigmay522, 2);
            break;
        case 523:
            dJydy[0] = 0.5*(-2*my523 + 2*y523)/pow(sigmay523, 2);
            break;
        case 524:
            dJydy[0] = 0.5*(-2*my524 + 2*y524)/pow(sigmay524, 2);
            break;
        case 525:
            dJydy[0] = 0.5*(-2*my525 + 2*y525)/pow(sigmay525, 2);
            break;
        case 526:
            dJydy[0] = 0.5*(-2*my526 + 2*y526)/pow(sigmay526, 2);
            break;
        case 527:
            dJydy[0] = 0.5*(-2*my527 + 2*y527)/pow(sigmay527, 2);
            break;
        case 528:
            dJydy[0] = 0.5*(-2*my528 + 2*y528)/pow(sigmay528, 2);
            break;
        case 529:
            dJydy[0] = 0.5*(-2*my529 + 2*y529)/pow(sigmay529, 2);
            break;
        case 530:
            dJydy[0] = 0.5*(-2*my530 + 2*y530)/pow(sigmay530, 2);
            break;
        case 531:
            dJydy[0] = 0.5*(-2*my531 + 2*y531)/pow(sigmay531, 2);
            break;
        case 532:
            dJydy[0] = 0.5*(-2*my532 + 2*y532)/pow(sigmay532, 2);
            break;
        case 533:
            dJydy[0] = 0.5*(-2*my533 + 2*y533)/pow(sigmay533, 2);
            break;
        case 534:
            dJydy[0] = 0.5*(-2*my534 + 2*y534)/pow(sigmay534, 2);
            break;
        case 535:
            dJydy[0] = 0.5*(-2*my535 + 2*y535)/pow(sigmay535, 2);
            break;
        case 536:
            dJydy[0] = 0.5*(-2*my536 + 2*y536)/pow(sigmay536, 2);
            break;
        case 537:
            dJydy[0] = 0.5*(-2*my537 + 2*y537)/pow(sigmay537, 2);
            break;
        case 538:
            dJydy[0] = 0.5*(-2*my538 + 2*y538)/pow(sigmay538, 2);
            break;
        case 539:
            dJydy[0] = 0.5*(-2*my539 + 2*y539)/pow(sigmay539, 2);
            break;
        case 540:
            dJydy[0] = 0.5*(-2*my540 + 2*y540)/pow(sigmay540, 2);
            break;
        case 541:
            dJydy[0] = 0.5*(-2*my541 + 2*y541)/pow(sigmay541, 2);
            break;
        case 542:
            dJydy[0] = 0.5*(-2*my542 + 2*y542)/pow(sigmay542, 2);
            break;
        case 543:
            dJydy[0] = 0.5*(-2*my543 + 2*y543)/pow(sigmay543, 2);
            break;
        case 544:
            dJydy[0] = 0.5*(-2*my544 + 2*y544)/pow(sigmay544, 2);
            break;
        case 545:
            dJydy[0] = 0.5*(-2*my545 + 2*y545)/pow(sigmay545, 2);
            break;
        case 546:
            dJydy[0] = 0.5*(-2*my546 + 2*y546)/pow(sigmay546, 2);
            break;
        case 547:
            dJydy[0] = 0.5*(-2*my547 + 2*y547)/pow(sigmay547, 2);
            break;
        case 548:
            dJydy[0] = 0.5*(-2*my548 + 2*y548)/pow(sigmay548, 2);
            break;
        case 549:
            dJydy[0] = 0.5*(-2*my549 + 2*y549)/pow(sigmay549, 2);
            break;
        case 550:
            dJydy[0] = 0.5*(-2*my550 + 2*y550)/pow(sigmay550, 2);
            break;
        case 551:
            dJydy[0] = 0.5*(-2*my551 + 2*y551)/pow(sigmay551, 2);
            break;
        case 552:
            dJydy[0] = 0.5*(-2*my552 + 2*y552)/pow(sigmay552, 2);
            break;
        case 553:
            dJydy[0] = 0.5*(-2*my553 + 2*y553)/pow(sigmay553, 2);
            break;
        case 554:
            dJydy[0] = 0.5*(-2*my554 + 2*y554)/pow(sigmay554, 2);
            break;
        case 555:
            dJydy[0] = 0.5*(-2*my555 + 2*y555)/pow(sigmay555, 2);
            break;
        case 556:
            dJydy[0] = 0.5*(-2*my556 + 2*y556)/pow(sigmay556, 2);
            break;
        case 557:
            dJydy[0] = 0.5*(-2*my557 + 2*y557)/pow(sigmay557, 2);
            break;
        case 558:
            dJydy[0] = 0.5*(-2*my558 + 2*y558)/pow(sigmay558, 2);
            break;
        case 559:
            dJydy[0] = 0.5*(-2*my559 + 2*y559)/pow(sigmay559, 2);
            break;
        case 560:
            dJydy[0] = 0.5*(-2*my560 + 2*y560)/pow(sigmay560, 2);
            break;
        case 561:
            dJydy[0] = 0.5*(-2*my561 + 2*y561)/pow(sigmay561, 2);
            break;
        case 562:
            dJydy[0] = 0.5*(-2*my562 + 2*y562)/pow(sigmay562, 2);
            break;
        case 563:
            dJydy[0] = 0.5*(-2*my563 + 2*y563)/pow(sigmay563, 2);
            break;
        case 564:
            dJydy[0] = 0.5*(-2*my564 + 2*y564)/pow(sigmay564, 2);
            break;
        case 565:
            dJydy[0] = 0.5*(-2*my565 + 2*y565)/pow(sigmay565, 2);
            break;
        case 566:
            dJydy[0] = 0.5*(-2*my566 + 2*y566)/pow(sigmay566, 2);
            break;
        case 567:
            dJydy[0] = 0.5*(-2*my567 + 2*y567)/pow(sigmay567, 2);
            break;
        case 568:
            dJydy[0] = 0.5*(-2*my568 + 2*y568)/pow(sigmay568, 2);
            break;
        case 569:
            dJydy[0] = 0.5*(-2*my569 + 2*y569)/pow(sigmay569, 2);
            break;
        case 570:
            dJydy[0] = 0.5*(-2*my570 + 2*y570)/pow(sigmay570, 2);
            break;
        case 571:
            dJydy[0] = 0.5*(-2*my571 + 2*y571)/pow(sigmay571, 2);
            break;
        case 572:
            dJydy[0] = 0.5*(-2*my572 + 2*y572)/pow(sigmay572, 2);
            break;
        case 573:
            dJydy[0] = 0.5*(-2*my573 + 2*y573)/pow(sigmay573, 2);
            break;
        case 574:
            dJydy[0] = 0.5*(-2*my574 + 2*y574)/pow(sigmay574, 2);
            break;
        case 575:
            dJydy[0] = 0.5*(-2*my575 + 2*y575)/pow(sigmay575, 2);
            break;
        case 576:
            dJydy[0] = 0.5*(-2*my576 + 2*y576)/pow(sigmay576, 2);
            break;
        case 577:
            dJydy[0] = 0.5*(-2*my577 + 2*y577)/pow(sigmay577, 2);
            break;
        case 578:
            dJydy[0] = 0.5*(-2*my578 + 2*y578)/pow(sigmay578, 2);
            break;
        case 579:
            dJydy[0] = 0.5*(-2*my579 + 2*y579)/pow(sigmay579, 2);
            break;
        case 580:
            dJydy[0] = 0.5*(-2*my580 + 2*y580)/pow(sigmay580, 2);
            break;
        case 581:
            dJydy[0] = 0.5*(-2*my581 + 2*y581)/pow(sigmay581, 2);
            break;
        case 582:
            dJydy[0] = 0.5*(-2*my582 + 2*y582)/pow(sigmay582, 2);
            break;
        case 583:
            dJydy[0] = 0.5*(-2*my583 + 2*y583)/pow(sigmay583, 2);
            break;
        case 584:
            dJydy[0] = 0.5*(-2*my584 + 2*y584)/pow(sigmay584, 2);
            break;
        case 585:
            dJydy[0] = 0.5*(-2*my585 + 2*y585)/pow(sigmay585, 2);
            break;
        case 586:
            dJydy[0] = 0.5*(-2*my586 + 2*y586)/pow(sigmay586, 2);
            break;
        case 587:
            dJydy[0] = 0.5*(-2*my587 + 2*y587)/pow(sigmay587, 2);
            break;
        case 588:
            dJydy[0] = 0.5*(-2*my588 + 2*y588)/pow(sigmay588, 2);
            break;
        case 589:
            dJydy[0] = 0.5*(-2*my589 + 2*y589)/pow(sigmay589, 2);
            break;
        case 590:
            dJydy[0] = 0.5*(-2*my590 + 2*y590)/pow(sigmay590, 2);
            break;
        case 591:
            dJydy[0] = 0.5*(-2*my591 + 2*y591)/pow(sigmay591, 2);
            break;
        case 592:
            dJydy[0] = 0.5*(-2*my592 + 2*y592)/pow(sigmay592, 2);
            break;
        case 593:
            dJydy[0] = 0.5*(-2*my593 + 2*y593)/pow(sigmay593, 2);
            break;
        case 594:
            dJydy[0] = 0.5*(-2*my594 + 2*y594)/pow(sigmay594, 2);
            break;
        case 595:
            dJydy[0] = 0.5*(-2*my595 + 2*y595)/pow(sigmay595, 2);
            break;
        case 596:
            dJydy[0] = 0.5*(-2*my596 + 2*y596)/pow(sigmay596, 2);
            break;
        case 597:
            dJydy[0] = 0.5*(-2*my597 + 2*y597)/pow(sigmay597, 2);
            break;
        case 598:
            dJydy[0] = 0.5*(-2*my598 + 2*y598)/pow(sigmay598, 2);
            break;
        case 599:
            dJydy[0] = 0.5*(-2*my599 + 2*y599)/pow(sigmay599, 2);
            break;
        case 600:
            dJydy[0] = 0.5*(-2*my600 + 2*y600)/pow(sigmay600, 2);
            break;
        case 601:
            dJydy[0] = 0.5*(-2*my601 + 2*y601)/pow(sigmay601, 2);
            break;
        case 602:
            dJydy[0] = 0.5*(-2*my602 + 2*y602)/pow(sigmay602, 2);
            break;
        case 603:
            dJydy[0] = 0.5*(-2*my603 + 2*y603)/pow(sigmay603, 2);
            break;
        case 604:
            dJydy[0] = 0.5*(-2*my604 + 2*y604)/pow(sigmay604, 2);
            break;
        case 605:
            dJydy[0] = 0.5*(-2*my605 + 2*y605)/pow(sigmay605, 2);
            break;
        case 606:
            dJydy[0] = 0.5*(-2*my606 + 2*y606)/pow(sigmay606, 2);
            break;
        case 607:
            dJydy[0] = 0.5*(-2*my607 + 2*y607)/pow(sigmay607, 2);
            break;
        case 608:
            dJydy[0] = 0.5*(-2*my608 + 2*y608)/pow(sigmay608, 2);
            break;
        case 609:
            dJydy[0] = 0.5*(-2*my609 + 2*y609)/pow(sigmay609, 2);
            break;
        case 610:
            dJydy[0] = 0.5*(-2*my610 + 2*y610)/pow(sigmay610, 2);
            break;
        case 611:
            dJydy[0] = 0.5*(-2*my611 + 2*y611)/pow(sigmay611, 2);
            break;
        case 612:
            dJydy[0] = 0.5*(-2*my612 + 2*y612)/pow(sigmay612, 2);
            break;
        case 613:
            dJydy[0] = 0.5*(-2*my613 + 2*y613)/pow(sigmay613, 2);
            break;
        case 614:
            dJydy[0] = 0.5*(-2*my614 + 2*y614)/pow(sigmay614, 2);
            break;
        case 615:
            dJydy[0] = 0.5*(-2*my615 + 2*y615)/pow(sigmay615, 2);
            break;
        case 616:
            dJydy[0] = 0.5*(-2*my616 + 2*y616)/pow(sigmay616, 2);
            break;
        case 617:
            dJydy[0] = 0.5*(-2*my617 + 2*y617)/pow(sigmay617, 2);
            break;
        case 618:
            dJydy[0] = 0.5*(-2*my618 + 2*y618)/pow(sigmay618, 2);
            break;
        case 619:
            dJydy[0] = 0.5*(-2*my619 + 2*y619)/pow(sigmay619, 2);
            break;
        case 620:
            dJydy[0] = 0.5*(-2*my620 + 2*y620)/pow(sigmay620, 2);
            break;
        case 621:
            dJydy[0] = 0.5*(-2*my621 + 2*y621)/pow(sigmay621, 2);
            break;
        case 622:
            dJydy[0] = 0.5*(-2*my622 + 2*y622)/pow(sigmay622, 2);
            break;
        case 623:
            dJydy[0] = 0.5*(-2*my623 + 2*y623)/pow(sigmay623, 2);
            break;
        case 624:
            dJydy[0] = 0.5*(-2*my624 + 2*y624)/pow(sigmay624, 2);
            break;
        case 625:
            dJydy[0] = 0.5*(-2*my625 + 2*y625)/pow(sigmay625, 2);
            break;
        case 626:
            dJydy[0] = 0.5*(-2*my626 + 2*y626)/pow(sigmay626, 2);
            break;
        case 627:
            dJydy[0] = 0.5*(-2*my627 + 2*y627)/pow(sigmay627, 2);
            break;
        case 628:
            dJydy[0] = 0.5*(-2*my628 + 2*y628)/pow(sigmay628, 2);
            break;
        case 629:
            dJydy[0] = 0.5*(-2*my629 + 2*y629)/pow(sigmay629, 2);
            break;
        case 630:
            dJydy[0] = 0.5*(-2*my630 + 2*y630)/pow(sigmay630, 2);
            break;
        case 631:
            dJydy[0] = 0.5*(-2*my631 + 2*y631)/pow(sigmay631, 2);
            break;
        case 632:
            dJydy[0] = 0.5*(-2*my632 + 2*y632)/pow(sigmay632, 2);
            break;
        case 633:
            dJydy[0] = 0.5*(-2*my633 + 2*y633)/pow(sigmay633, 2);
            break;
        case 634:
            dJydy[0] = 0.5*(-2*my634 + 2*y634)/pow(sigmay634, 2);
            break;
        case 635:
            dJydy[0] = 0.5*(-2*my635 + 2*y635)/pow(sigmay635, 2);
            break;
        case 636:
            dJydy[0] = 0.5*(-2*my636 + 2*y636)/pow(sigmay636, 2);
            break;
        case 637:
            dJydy[0] = 0.5*(-2*my637 + 2*y637)/pow(sigmay637, 2);
            break;
        case 638:
            dJydy[0] = 0.5*(-2*my638 + 2*y638)/pow(sigmay638, 2);
            break;
        case 639:
            dJydy[0] = 0.5*(-2*my639 + 2*y639)/pow(sigmay639, 2);
            break;
        case 640:
            dJydy[0] = 0.5*(-2*my640 + 2*y640)/pow(sigmay640, 2);
            break;
        case 641:
            dJydy[0] = 0.5*(-2*my641 + 2*y641)/pow(sigmay641, 2);
            break;
        case 642:
            dJydy[0] = 0.5*(-2*my642 + 2*y642)/pow(sigmay642, 2);
            break;
        case 643:
            dJydy[0] = 0.5*(-2*my643 + 2*y643)/pow(sigmay643, 2);
            break;
        case 644:
            dJydy[0] = 0.5*(-2*my644 + 2*y644)/pow(sigmay644, 2);
            break;
        case 645:
            dJydy[0] = 0.5*(-2*my645 + 2*y645)/pow(sigmay645, 2);
            break;
        case 646:
            dJydy[0] = 0.5*(-2*my646 + 2*y646)/pow(sigmay646, 2);
            break;
        case 647:
            dJydy[0] = 0.5*(-2*my647 + 2*y647)/pow(sigmay647, 2);
            break;
        case 648:
            dJydy[0] = 0.5*(-2*my648 + 2*y648)/pow(sigmay648, 2);
            break;
        case 649:
            dJydy[0] = 0.5*(-2*my649 + 2*y649)/pow(sigmay649, 2);
            break;
        case 650:
            dJydy[0] = 0.5*(-2*my650 + 2*y650)/pow(sigmay650, 2);
            break;
        case 651:
            dJydy[0] = 0.5*(-2*my651 + 2*y651)/pow(sigmay651, 2);
            break;
        case 652:
            dJydy[0] = 0.5*(-2*my652 + 2*y652)/pow(sigmay652, 2);
            break;
        case 653:
            dJydy[0] = 0.5*(-2*my653 + 2*y653)/pow(sigmay653, 2);
            break;
        case 654:
            dJydy[0] = 0.5*(-2*my654 + 2*y654)/pow(sigmay654, 2);
            break;
        case 655:
            dJydy[0] = 0.5*(-2*my655 + 2*y655)/pow(sigmay655, 2);
            break;
        case 656:
            dJydy[0] = 0.5*(-2*my656 + 2*y656)/pow(sigmay656, 2);
            break;
        case 657:
            dJydy[0] = 0.5*(-2*my657 + 2*y657)/pow(sigmay657, 2);
            break;
        case 658:
            dJydy[0] = 0.5*(-2*my658 + 2*y658)/pow(sigmay658, 2);
            break;
        case 659:
            dJydy[0] = 0.5*(-2*my659 + 2*y659)/pow(sigmay659, 2);
            break;
        case 660:
            dJydy[0] = 0.5*(-2*my660 + 2*y660)/pow(sigmay660, 2);
            break;
        case 661:
            dJydy[0] = 0.5*(-2*my661 + 2*y661)/pow(sigmay661, 2);
            break;
        case 662:
            dJydy[0] = 0.5*(-2*my662 + 2*y662)/pow(sigmay662, 2);
            break;
        case 663:
            dJydy[0] = 0.5*(-2*my663 + 2*y663)/pow(sigmay663, 2);
            break;
        case 664:
            dJydy[0] = 0.5*(-2*my664 + 2*y664)/pow(sigmay664, 2);
            break;
        case 665:
            dJydy[0] = 0.5*(-2*my665 + 2*y665)/pow(sigmay665, 2);
            break;
        case 666:
            dJydy[0] = 0.5*(-2*my666 + 2*y666)/pow(sigmay666, 2);
            break;
        case 667:
            dJydy[0] = 0.5*(-2*my667 + 2*y667)/pow(sigmay667, 2);
            break;
        case 668:
            dJydy[0] = 0.5*(-2*my668 + 2*y668)/pow(sigmay668, 2);
            break;
        case 669:
            dJydy[0] = 0.5*(-2*my669 + 2*y669)/pow(sigmay669, 2);
            break;
        case 670:
            dJydy[0] = 0.5*(-2*my670 + 2*y670)/pow(sigmay670, 2);
            break;
        case 671:
            dJydy[0] = 0.5*(-2*my671 + 2*y671)/pow(sigmay671, 2);
            break;
        case 672:
            dJydy[0] = 0.5*(-2*my672 + 2*y672)/pow(sigmay672, 2);
            break;
        case 673:
            dJydy[0] = 0.5*(-2*my673 + 2*y673)/pow(sigmay673, 2);
            break;
        case 674:
            dJydy[0] = 0.5*(-2*my674 + 2*y674)/pow(sigmay674, 2);
            break;
        case 675:
            dJydy[0] = 0.5*(-2*my675 + 2*y675)/pow(sigmay675, 2);
            break;
        case 676:
            dJydy[0] = 0.5*(-2*my676 + 2*y676)/pow(sigmay676, 2);
            break;
        case 677:
            dJydy[0] = 0.5*(-2*my677 + 2*y677)/pow(sigmay677, 2);
            break;
        case 678:
            dJydy[0] = 0.5*(-2*my678 + 2*y678)/pow(sigmay678, 2);
            break;
        case 679:
            dJydy[0] = 0.5*(-2*my679 + 2*y679)/pow(sigmay679, 2);
            break;
        case 680:
            dJydy[0] = 0.5*(-2*my680 + 2*y680)/pow(sigmay680, 2);
            break;
        case 681:
            dJydy[0] = 0.5*(-2*my681 + 2*y681)/pow(sigmay681, 2);
            break;
        case 682:
            dJydy[0] = 0.5*(-2*my682 + 2*y682)/pow(sigmay682, 2);
            break;
        case 683:
            dJydy[0] = 0.5*(-2*my683 + 2*y683)/pow(sigmay683, 2);
            break;
        case 684:
            dJydy[0] = 0.5*(-2*my684 + 2*y684)/pow(sigmay684, 2);
            break;
        case 685:
            dJydy[0] = 0.5*(-2*my685 + 2*y685)/pow(sigmay685, 2);
            break;
        case 686:
            dJydy[0] = 0.5*(-2*my686 + 2*y686)/pow(sigmay686, 2);
            break;
        case 687:
            dJydy[0] = 0.5*(-2*my687 + 2*y687)/pow(sigmay687, 2);
            break;
        case 688:
            dJydy[0] = 0.5*(-2*my688 + 2*y688)/pow(sigmay688, 2);
            break;
        case 689:
            dJydy[0] = 0.5*(-2*my689 + 2*y689)/pow(sigmay689, 2);
            break;
        case 690:
            dJydy[0] = 0.5*(-2*my690 + 2*y690)/pow(sigmay690, 2);
            break;
        case 691:
            dJydy[0] = 0.5*(-2*my691 + 2*y691)/pow(sigmay691, 2);
            break;
        case 692:
            dJydy[0] = 0.5*(-2*my692 + 2*y692)/pow(sigmay692, 2);
            break;
        case 693:
            dJydy[0] = 0.5*(-2*my693 + 2*y693)/pow(sigmay693, 2);
            break;
        case 694:
            dJydy[0] = 0.5*(-2*my694 + 2*y694)/pow(sigmay694, 2);
            break;
        case 695:
            dJydy[0] = 0.5*(-2*my695 + 2*y695)/pow(sigmay695, 2);
            break;
        case 696:
            dJydy[0] = 0.5*(-2*my696 + 2*y696)/pow(sigmay696, 2);
            break;
        case 697:
            dJydy[0] = 0.5*(-2*my697 + 2*y697)/pow(sigmay697, 2);
            break;
        case 698:
            dJydy[0] = 0.5*(-2*my698 + 2*y698)/pow(sigmay698, 2);
            break;
        case 699:
            dJydy[0] = 0.5*(-2*my699 + 2*y699)/pow(sigmay699, 2);
            break;
        case 700:
            dJydy[0] = 0.5*(-2*my700 + 2*y700)/pow(sigmay700, 2);
            break;
        case 701:
            dJydy[0] = 0.5*(-2*my701 + 2*y701)/pow(sigmay701, 2);
            break;
        case 702:
            dJydy[0] = 0.5*(-2*my702 + 2*y702)/pow(sigmay702, 2);
            break;
        case 703:
            dJydy[0] = 0.5*(-2*my703 + 2*y703)/pow(sigmay703, 2);
            break;
        case 704:
            dJydy[0] = 0.5*(-2*my704 + 2*y704)/pow(sigmay704, 2);
            break;
        case 705:
            dJydy[0] = 0.5*(-2*my705 + 2*y705)/pow(sigmay705, 2);
            break;
        case 706:
            dJydy[0] = 0.5*(-2*my706 + 2*y706)/pow(sigmay706, 2);
            break;
        case 707:
            dJydy[0] = 0.5*(-2*my707 + 2*y707)/pow(sigmay707, 2);
            break;
        case 708:
            dJydy[0] = 0.5*(-2*my708 + 2*y708)/pow(sigmay708, 2);
            break;
        case 709:
            dJydy[0] = 0.5*(-2*my709 + 2*y709)/pow(sigmay709, 2);
            break;
        case 710:
            dJydy[0] = 0.5*(-2*my710 + 2*y710)/pow(sigmay710, 2);
            break;
        case 711:
            dJydy[0] = 0.5*(-2*my711 + 2*y711)/pow(sigmay711, 2);
            break;
        case 712:
            dJydy[0] = 0.5*(-2*my712 + 2*y712)/pow(sigmay712, 2);
            break;
        case 713:
            dJydy[0] = 0.5*(-2*my713 + 2*y713)/pow(sigmay713, 2);
            break;
        case 714:
            dJydy[0] = 0.5*(-2*my714 + 2*y714)/pow(sigmay714, 2);
            break;
        case 715:
            dJydy[0] = 0.5*(-2*my715 + 2*y715)/pow(sigmay715, 2);
            break;
        case 716:
            dJydy[0] = 0.5*(-2*my716 + 2*y716)/pow(sigmay716, 2);
            break;
        case 717:
            dJydy[0] = 0.5*(-2*my717 + 2*y717)/pow(sigmay717, 2);
            break;
        case 718:
            dJydy[0] = 0.5*(-2*my718 + 2*y718)/pow(sigmay718, 2);
            break;
        case 719:
            dJydy[0] = 0.5*(-2*my719 + 2*y719)/pow(sigmay719, 2);
            break;
        case 720:
            dJydy[0] = 0.5*(-2*my720 + 2*y720)/pow(sigmay720, 2);
            break;
        case 721:
            dJydy[0] = 0.5*(-2*my721 + 2*y721)/pow(sigmay721, 2);
            break;
        case 722:
            dJydy[0] = 0.5*(-2*my722 + 2*y722)/pow(sigmay722, 2);
            break;
        case 723:
            dJydy[0] = 0.5*(-2*my723 + 2*y723)/pow(sigmay723, 2);
            break;
        case 724:
            dJydy[0] = 0.5*(-2*my724 + 2*y724)/pow(sigmay724, 2);
            break;
        case 725:
            dJydy[0] = 0.5*(-2*my725 + 2*y725)/pow(sigmay725, 2);
            break;
        case 726:
            dJydy[0] = 0.5*(-2*my726 + 2*y726)/pow(sigmay726, 2);
            break;
        case 727:
            dJydy[0] = 0.5*(-2*my727 + 2*y727)/pow(sigmay727, 2);
            break;
        case 728:
            dJydy[0] = 0.5*(-2*my728 + 2*y728)/pow(sigmay728, 2);
            break;
        case 729:
            dJydy[0] = 0.5*(-2*my729 + 2*y729)/pow(sigmay729, 2);
            break;
        case 730:
            dJydy[0] = 0.5*(-2*my730 + 2*y730)/pow(sigmay730, 2);
            break;
        case 731:
            dJydy[0] = 0.5*(-2*my731 + 2*y731)/pow(sigmay731, 2);
            break;
        case 732:
            dJydy[0] = 0.5*(-2*my732 + 2*y732)/pow(sigmay732, 2);
            break;
        case 733:
            dJydy[0] = 0.5*(-2*my733 + 2*y733)/pow(sigmay733, 2);
            break;
        case 734:
            dJydy[0] = 0.5*(-2*my734 + 2*y734)/pow(sigmay734, 2);
            break;
        case 735:
            dJydy[0] = 0.5*(-2*my735 + 2*y735)/pow(sigmay735, 2);
            break;
        case 736:
            dJydy[0] = 0.5*(-2*my736 + 2*y736)/pow(sigmay736, 2);
            break;
        case 737:
            dJydy[0] = 0.5*(-2*my737 + 2*y737)/pow(sigmay737, 2);
            break;
        case 738:
            dJydy[0] = 0.5*(-2*my738 + 2*y738)/pow(sigmay738, 2);
            break;
        case 739:
            dJydy[0] = 0.5*(-2*my739 + 2*y739)/pow(sigmay739, 2);
            break;
        case 740:
            dJydy[0] = 0.5*(-2*my740 + 2*y740)/pow(sigmay740, 2);
            break;
        case 741:
            dJydy[0] = 0.5*(-2*my741 + 2*y741)/pow(sigmay741, 2);
            break;
        case 742:
            dJydy[0] = 0.5*(-2*my742 + 2*y742)/pow(sigmay742, 2);
            break;
        case 743:
            dJydy[0] = 0.5*(-2*my743 + 2*y743)/pow(sigmay743, 2);
            break;
        case 744:
            dJydy[0] = 0.5*(-2*my744 + 2*y744)/pow(sigmay744, 2);
            break;
        case 745:
            dJydy[0] = 0.5*(-2*my745 + 2*y745)/pow(sigmay745, 2);
            break;
        case 746:
            dJydy[0] = 0.5*(-2*my746 + 2*y746)/pow(sigmay746, 2);
            break;
        case 747:
            dJydy[0] = 0.5*(-2*my747 + 2*y747)/pow(sigmay747, 2);
            break;
        case 748:
            dJydy[0] = 0.5*(-2*my748 + 2*y748)/pow(sigmay748, 2);
            break;
        case 749:
            dJydy[0] = 0.5*(-2*my749 + 2*y749)/pow(sigmay749, 2);
            break;
        case 750:
            dJydy[0] = 0.5*(-2*my750 + 2*y750)/pow(sigmay750, 2);
            break;
        case 751:
            dJydy[0] = 0.5*(-2*my751 + 2*y751)/pow(sigmay751, 2);
            break;
        case 752:
            dJydy[0] = 0.5*(-2*my752 + 2*y752)/pow(sigmay752, 2);
            break;
        case 753:
            dJydy[0] = 0.5*(-2*my753 + 2*y753)/pow(sigmay753, 2);
            break;
        case 754:
            dJydy[0] = 0.5*(-2*my754 + 2*y754)/pow(sigmay754, 2);
            break;
        case 755:
            dJydy[0] = 0.5*(-2*my755 + 2*y755)/pow(sigmay755, 2);
            break;
        case 756:
            dJydy[0] = 0.5*(-2*my756 + 2*y756)/pow(sigmay756, 2);
            break;
        case 757:
            dJydy[0] = 0.5*(-2*my757 + 2*y757)/pow(sigmay757, 2);
            break;
        case 758:
            dJydy[0] = 0.5*(-2*my758 + 2*y758)/pow(sigmay758, 2);
            break;
        case 759:
            dJydy[0] = 0.5*(-2*my759 + 2*y759)/pow(sigmay759, 2);
            break;
        case 760:
            dJydy[0] = 0.5*(-2*my760 + 2*y760)/pow(sigmay760, 2);
            break;
        case 761:
            dJydy[0] = 0.5*(-2*my761 + 2*y761)/pow(sigmay761, 2);
            break;
        case 762:
            dJydy[0] = 0.5*(-2*my762 + 2*y762)/pow(sigmay762, 2);
            break;
        case 763:
            dJydy[0] = 0.5*(-2*my763 + 2*y763)/pow(sigmay763, 2);
            break;
        case 764:
            dJydy[0] = 0.5*(-2*my764 + 2*y764)/pow(sigmay764, 2);
            break;
        case 765:
            dJydy[0] = 0.5*(-2*my765 + 2*y765)/pow(sigmay765, 2);
            break;
        case 766:
            dJydy[0] = 0.5*(-2*my766 + 2*y766)/pow(sigmay766, 2);
            break;
        case 767:
            dJydy[0] = 0.5*(-2*my767 + 2*y767)/pow(sigmay767, 2);
            break;
        case 768:
            dJydy[0] = 0.5*(-2*my768 + 2*y768)/pow(sigmay768, 2);
            break;
        case 769:
            dJydy[0] = 0.5*(-2*my769 + 2*y769)/pow(sigmay769, 2);
            break;
        case 770:
            dJydy[0] = 0.5*(-2*my770 + 2*y770)/pow(sigmay770, 2);
            break;
        case 771:
            dJydy[0] = 0.5*(-2*my771 + 2*y771)/pow(sigmay771, 2);
            break;
        case 772:
            dJydy[0] = 0.5*(-2*my772 + 2*y772)/pow(sigmay772, 2);
            break;
        case 773:
            dJydy[0] = 0.5*(-2*my773 + 2*y773)/pow(sigmay773, 2);
            break;
        case 774:
            dJydy[0] = 0.5*(-2*my774 + 2*y774)/pow(sigmay774, 2);
            break;
        case 775:
            dJydy[0] = 0.5*(-2*my775 + 2*y775)/pow(sigmay775, 2);
            break;
        case 776:
            dJydy[0] = 0.5*(-2*my776 + 2*y776)/pow(sigmay776, 2);
            break;
        case 777:
            dJydy[0] = 0.5*(-2*my777 + 2*y777)/pow(sigmay777, 2);
            break;
        case 778:
            dJydy[0] = 0.5*(-2*my778 + 2*y778)/pow(sigmay778, 2);
            break;
        case 779:
            dJydy[0] = 0.5*(-2*my779 + 2*y779)/pow(sigmay779, 2);
            break;
        case 780:
            dJydy[0] = 0.5*(-2*my780 + 2*y780)/pow(sigmay780, 2);
            break;
        case 781:
            dJydy[0] = 0.5*(-2*my781 + 2*y781)/pow(sigmay781, 2);
            break;
        case 782:
            dJydy[0] = 0.5*(-2*my782 + 2*y782)/pow(sigmay782, 2);
            break;
        case 783:
            dJydy[0] = 0.5*(-2*my783 + 2*y783)/pow(sigmay783, 2);
            break;
        case 784:
            dJydy[0] = 0.5*(-2*my784 + 2*y784)/pow(sigmay784, 2);
            break;
        case 785:
            dJydy[0] = 0.5*(-2*my785 + 2*y785)/pow(sigmay785, 2);
            break;
        case 786:
            dJydy[0] = 0.5*(-2*my786 + 2*y786)/pow(sigmay786, 2);
            break;
        case 787:
            dJydy[0] = 0.5*(-2*my787 + 2*y787)/pow(sigmay787, 2);
            break;
        case 788:
            dJydy[0] = 0.5*(-2*my788 + 2*y788)/pow(sigmay788, 2);
            break;
        case 789:
            dJydy[0] = 0.5*(-2*my789 + 2*y789)/pow(sigmay789, 2);
            break;
        case 790:
            dJydy[0] = 0.5*(-2*my790 + 2*y790)/pow(sigmay790, 2);
            break;
        case 791:
            dJydy[0] = 0.5*(-2*my791 + 2*y791)/pow(sigmay791, 2);
            break;
        case 792:
            dJydy[0] = 0.5*(-2*my792 + 2*y792)/pow(sigmay792, 2);
            break;
        case 793:
            dJydy[0] = 0.5*(-2*my793 + 2*y793)/pow(sigmay793, 2);
            break;
        case 794:
            dJydy[0] = 0.5*(-2*my794 + 2*y794)/pow(sigmay794, 2);
            break;
        case 795:
            dJydy[0] = 0.5*(-2*my795 + 2*y795)/pow(sigmay795, 2);
            break;
        case 796:
            dJydy[0] = 0.5*(-2*my796 + 2*y796)/pow(sigmay796, 2);
            break;
        case 797:
            dJydy[0] = 0.5*(-2*my797 + 2*y797)/pow(sigmay797, 2);
            break;
        case 798:
            dJydy[0] = 0.5*(-2*my798 + 2*y798)/pow(sigmay798, 2);
            break;
        case 799:
            dJydy[0] = 0.5*(-2*my799 + 2*y799)/pow(sigmay799, 2);
            break;
        case 800:
            dJydy[0] = 0.5*(-2*my800 + 2*y800)/pow(sigmay800, 2);
            break;
        case 801:
            dJydy[0] = 0.5*(-2*my801 + 2*y801)/pow(sigmay801, 2);
            break;
        case 802:
            dJydy[0] = 0.5*(-2*my802 + 2*y802)/pow(sigmay802, 2);
            break;
        case 803:
            dJydy[0] = 0.5*(-2*my803 + 2*y803)/pow(sigmay803, 2);
            break;
        case 804:
            dJydy[0] = 0.5*(-2*my804 + 2*y804)/pow(sigmay804, 2);
            break;
        case 805:
            dJydy[0] = 0.5*(-2*my805 + 2*y805)/pow(sigmay805, 2);
            break;
        case 806:
            dJydy[0] = 0.5*(-2*my806 + 2*y806)/pow(sigmay806, 2);
            break;
        case 807:
            dJydy[0] = 0.5*(-2*my807 + 2*y807)/pow(sigmay807, 2);
            break;
        case 808:
            dJydy[0] = 0.5*(-2*my808 + 2*y808)/pow(sigmay808, 2);
            break;
        case 809:
            dJydy[0] = 0.5*(-2*my809 + 2*y809)/pow(sigmay809, 2);
            break;
        case 810:
            dJydy[0] = 0.5*(-2*my810 + 2*y810)/pow(sigmay810, 2);
            break;
        case 811:
            dJydy[0] = 0.5*(-2*my811 + 2*y811)/pow(sigmay811, 2);
            break;
        case 812:
            dJydy[0] = 0.5*(-2*my812 + 2*y812)/pow(sigmay812, 2);
            break;
        case 813:
            dJydy[0] = 0.5*(-2*my813 + 2*y813)/pow(sigmay813, 2);
            break;
        case 814:
            dJydy[0] = 0.5*(-2*my814 + 2*y814)/pow(sigmay814, 2);
            break;
        case 815:
            dJydy[0] = 0.5*(-2*my815 + 2*y815)/pow(sigmay815, 2);
            break;
        case 816:
            dJydy[0] = 0.5*(-2*my816 + 2*y816)/pow(sigmay816, 2);
            break;
        case 817:
            dJydy[0] = 0.5*(-2*my817 + 2*y817)/pow(sigmay817, 2);
            break;
        case 818:
            dJydy[0] = 0.5*(-2*my818 + 2*y818)/pow(sigmay818, 2);
            break;
        case 819:
            dJydy[0] = 0.5*(-2*my819 + 2*y819)/pow(sigmay819, 2);
            break;
        case 820:
            dJydy[0] = 0.5*(-2*my820 + 2*y820)/pow(sigmay820, 2);
            break;
        case 821:
            dJydy[0] = 0.5*(-2*my821 + 2*y821)/pow(sigmay821, 2);
            break;
        case 822:
            dJydy[0] = 0.5*(-2*my822 + 2*y822)/pow(sigmay822, 2);
            break;
        case 823:
            dJydy[0] = 0.5*(-2*my823 + 2*y823)/pow(sigmay823, 2);
            break;
        case 824:
            dJydy[0] = 0.5*(-2*my824 + 2*y824)/pow(sigmay824, 2);
            break;
        case 825:
            dJydy[0] = 0.5*(-2*my825 + 2*y825)/pow(sigmay825, 2);
            break;
        case 826:
            dJydy[0] = 0.5*(-2*my826 + 2*y826)/pow(sigmay826, 2);
            break;
        case 827:
            dJydy[0] = 0.5*(-2*my827 + 2*y827)/pow(sigmay827, 2);
            break;
        case 828:
            dJydy[0] = 0.5*(-2*my828 + 2*y828)/pow(sigmay828, 2);
            break;
        case 829:
            dJydy[0] = 0.5*(-2*my829 + 2*y829)/pow(sigmay829, 2);
            break;
        case 830:
            dJydy[0] = 0.5*(-2*my830 + 2*y830)/pow(sigmay830, 2);
            break;
        case 831:
            dJydy[0] = 0.5*(-2*my831 + 2*y831)/pow(sigmay831, 2);
            break;
        case 832:
            dJydy[0] = 0.5*(-2*my832 + 2*y832)/pow(sigmay832, 2);
            break;
        case 833:
            dJydy[0] = 0.5*(-2*my833 + 2*y833)/pow(sigmay833, 2);
            break;
        case 834:
            dJydy[0] = 0.5*(-2*my834 + 2*y834)/pow(sigmay834, 2);
            break;
        case 835:
            dJydy[0] = 0.5*(-2*my835 + 2*y835)/pow(sigmay835, 2);
            break;
        case 836:
            dJydy[0] = 0.5*(-2*my836 + 2*y836)/pow(sigmay836, 2);
            break;
        case 837:
            dJydy[0] = 0.5*(-2*my837 + 2*y837)/pow(sigmay837, 2);
            break;
        case 838:
            dJydy[0] = 0.5*(-2*my838 + 2*y838)/pow(sigmay838, 2);
            break;
        case 839:
            dJydy[0] = 0.5*(-2*my839 + 2*y839)/pow(sigmay839, 2);
            break;
        case 840:
            dJydy[0] = 0.5*(-2*my840 + 2*y840)/pow(sigmay840, 2);
            break;
        case 841:
            dJydy[0] = 0.5*(-2*my841 + 2*y841)/pow(sigmay841, 2);
            break;
        case 842:
            dJydy[0] = 0.5*(-2*my842 + 2*y842)/pow(sigmay842, 2);
            break;
        case 843:
            dJydy[0] = 0.5*(-2*my843 + 2*y843)/pow(sigmay843, 2);
            break;
        case 844:
            dJydy[0] = 0.5*(-2*my844 + 2*y844)/pow(sigmay844, 2);
            break;
        case 845:
            dJydy[0] = 0.5*(-2*my845 + 2*y845)/pow(sigmay845, 2);
            break;
        case 846:
            dJydy[0] = 0.5*(-2*my846 + 2*y846)/pow(sigmay846, 2);
            break;
        case 847:
            dJydy[0] = 0.5*(-2*my847 + 2*y847)/pow(sigmay847, 2);
            break;
        case 848:
            dJydy[0] = 0.5*(-2*my848 + 2*y848)/pow(sigmay848, 2);
            break;
        case 849:
            dJydy[0] = 0.5*(-2*my849 + 2*y849)/pow(sigmay849, 2);
            break;
        case 850:
            dJydy[0] = 0.5*(-2*my850 + 2*y850)/pow(sigmay850, 2);
            break;
        case 851:
            dJydy[0] = 0.5*(-2*my851 + 2*y851)/pow(sigmay851, 2);
            break;
        case 852:
            dJydy[0] = 0.5*(-2*my852 + 2*y852)/pow(sigmay852, 2);
            break;
        case 853:
            dJydy[0] = 0.5*(-2*my853 + 2*y853)/pow(sigmay853, 2);
            break;
        case 854:
            dJydy[0] = 0.5*(-2*my854 + 2*y854)/pow(sigmay854, 2);
            break;
        case 855:
            dJydy[0] = 0.5*(-2*my855 + 2*y855)/pow(sigmay855, 2);
            break;
        case 856:
            dJydy[0] = 0.5*(-2*my856 + 2*y856)/pow(sigmay856, 2);
            break;
        case 857:
            dJydy[0] = 0.5*(-2*my857 + 2*y857)/pow(sigmay857, 2);
            break;
        case 858:
            dJydy[0] = 0.5*(-2*my858 + 2*y858)/pow(sigmay858, 2);
            break;
        case 859:
            dJydy[0] = 0.5*(-2*my859 + 2*y859)/pow(sigmay859, 2);
            break;
        case 860:
            dJydy[0] = 0.5*(-2*my860 + 2*y860)/pow(sigmay860, 2);
            break;
        case 861:
            dJydy[0] = 0.5*(-2*my861 + 2*y861)/pow(sigmay861, 2);
            break;
        case 862:
            dJydy[0] = 0.5*(-2*my862 + 2*y862)/pow(sigmay862, 2);
            break;
        case 863:
            dJydy[0] = 0.5*(-2*my863 + 2*y863)/pow(sigmay863, 2);
            break;
        case 864:
            dJydy[0] = 0.5*(-2*my864 + 2*y864)/pow(sigmay864, 2);
            break;
        case 865:
            dJydy[0] = 0.5*(-2*my865 + 2*y865)/pow(sigmay865, 2);
            break;
        case 866:
            dJydy[0] = 0.5*(-2*my866 + 2*y866)/pow(sigmay866, 2);
            break;
        case 867:
            dJydy[0] = 0.5*(-2*my867 + 2*y867)/pow(sigmay867, 2);
            break;
        case 868:
            dJydy[0] = 0.5*(-2*my868 + 2*y868)/pow(sigmay868, 2);
            break;
        case 869:
            dJydy[0] = 0.5*(-2*my869 + 2*y869)/pow(sigmay869, 2);
            break;
        case 870:
            dJydy[0] = 0.5*(-2*my870 + 2*y870)/pow(sigmay870, 2);
            break;
        case 871:
            dJydy[0] = 0.5*(-2*my871 + 2*y871)/pow(sigmay871, 2);
            break;
        case 872:
            dJydy[0] = 0.5*(-2*my872 + 2*y872)/pow(sigmay872, 2);
            break;
        case 873:
            dJydy[0] = 0.5*(-2*my873 + 2*y873)/pow(sigmay873, 2);
            break;
        case 874:
            dJydy[0] = 0.5*(-2*my874 + 2*y874)/pow(sigmay874, 2);
            break;
        case 875:
            dJydy[0] = 0.5*(-2*my875 + 2*y875)/pow(sigmay875, 2);
            break;
        case 876:
            dJydy[0] = 0.5*(-2*my876 + 2*y876)/pow(sigmay876, 2);
            break;
        case 877:
            dJydy[0] = 0.5*(-2*my877 + 2*y877)/pow(sigmay877, 2);
            break;
        case 878:
            dJydy[0] = 0.5*(-2*my878 + 2*y878)/pow(sigmay878, 2);
            break;
        case 879:
            dJydy[0] = 0.5*(-2*my879 + 2*y879)/pow(sigmay879, 2);
            break;
        case 880:
            dJydy[0] = 0.5*(-2*my880 + 2*y880)/pow(sigmay880, 2);
            break;
        case 881:
            dJydy[0] = 0.5*(-2*my881 + 2*y881)/pow(sigmay881, 2);
            break;
        case 882:
            dJydy[0] = 0.5*(-2*my882 + 2*y882)/pow(sigmay882, 2);
            break;
        case 883:
            dJydy[0] = 0.5*(-2*my883 + 2*y883)/pow(sigmay883, 2);
            break;
        case 884:
            dJydy[0] = 0.5*(-2*my884 + 2*y884)/pow(sigmay884, 2);
            break;
        case 885:
            dJydy[0] = 0.5*(-2*my885 + 2*y885)/pow(sigmay885, 2);
            break;
        case 886:
            dJydy[0] = 0.5*(-2*my886 + 2*y886)/pow(sigmay886, 2);
            break;
        case 887:
            dJydy[0] = 0.5*(-2*my887 + 2*y887)/pow(sigmay887, 2);
            break;
        case 888:
            dJydy[0] = 0.5*(-2*my888 + 2*y888)/pow(sigmay888, 2);
            break;
        case 889:
            dJydy[0] = 0.5*(-2*my889 + 2*y889)/pow(sigmay889, 2);
            break;
        case 890:
            dJydy[0] = 0.5*(-2*my890 + 2*y890)/pow(sigmay890, 2);
            break;
        case 891:
            dJydy[0] = 0.5*(-2*my891 + 2*y891)/pow(sigmay891, 2);
            break;
        case 892:
            dJydy[0] = 0.5*(-2*my892 + 2*y892)/pow(sigmay892, 2);
            break;
        case 893:
            dJydy[0] = 0.5*(-2*my893 + 2*y893)/pow(sigmay893, 2);
            break;
        case 894:
            dJydy[0] = 0.5*(-2*my894 + 2*y894)/pow(sigmay894, 2);
            break;
        case 895:
            dJydy[0] = 0.5*(-2*my895 + 2*y895)/pow(sigmay895, 2);
            break;
        case 896:
            dJydy[0] = 0.5*(-2*my896 + 2*y896)/pow(sigmay896, 2);
            break;
        case 897:
            dJydy[0] = 0.5*(-2*my897 + 2*y897)/pow(sigmay897, 2);
            break;
        case 898:
            dJydy[0] = 0.5*(-2*my898 + 2*y898)/pow(sigmay898, 2);
            break;
        case 899:
            dJydy[0] = 0.5*(-2*my899 + 2*y899)/pow(sigmay899, 2);
            break;
        case 900:
            dJydy[0] = 0.5*(-2*my900 + 2*y900)/pow(sigmay900, 2);
            break;
        case 901:
            dJydy[0] = 0.5*(-2*my901 + 2*y901)/pow(sigmay901, 2);
            break;
        case 902:
            dJydy[0] = 0.5*(-2*my902 + 2*y902)/pow(sigmay902, 2);
            break;
        case 903:
            dJydy[0] = 0.5*(-2*my903 + 2*y903)/pow(sigmay903, 2);
            break;
        case 904:
            dJydy[0] = 0.5*(-2*my904 + 2*y904)/pow(sigmay904, 2);
            break;
        case 905:
            dJydy[0] = 0.5*(-2*my905 + 2*y905)/pow(sigmay905, 2);
            break;
        case 906:
            dJydy[0] = 0.5*(-2*my906 + 2*y906)/pow(sigmay906, 2);
            break;
        case 907:
            dJydy[0] = 0.5*(-2*my907 + 2*y907)/pow(sigmay907, 2);
            break;
        case 908:
            dJydy[0] = 0.5*(-2*my908 + 2*y908)/pow(sigmay908, 2);
            break;
        case 909:
            dJydy[0] = 0.5*(-2*my909 + 2*y909)/pow(sigmay909, 2);
            break;
        case 910:
            dJydy[0] = 0.5*(-2*my910 + 2*y910)/pow(sigmay910, 2);
            break;
        case 911:
            dJydy[0] = 0.5*(-2*my911 + 2*y911)/pow(sigmay911, 2);
            break;
        case 912:
            dJydy[0] = 0.5*(-2*my912 + 2*y912)/pow(sigmay912, 2);
            break;
        case 913:
            dJydy[0] = 0.5*(-2*my913 + 2*y913)/pow(sigmay913, 2);
            break;
        case 914:
            dJydy[0] = 0.5*(-2*my914 + 2*y914)/pow(sigmay914, 2);
            break;
        case 915:
            dJydy[0] = 0.5*(-2*my915 + 2*y915)/pow(sigmay915, 2);
            break;
        case 916:
            dJydy[0] = 0.5*(-2*my916 + 2*y916)/pow(sigmay916, 2);
            break;
        case 917:
            dJydy[0] = 0.5*(-2*my917 + 2*y917)/pow(sigmay917, 2);
            break;
        case 918:
            dJydy[0] = 0.5*(-2*my918 + 2*y918)/pow(sigmay918, 2);
            break;
        case 919:
            dJydy[0] = 0.5*(-2*my919 + 2*y919)/pow(sigmay919, 2);
            break;
        case 920:
            dJydy[0] = 0.5*(-2*my920 + 2*y920)/pow(sigmay920, 2);
            break;
        case 921:
            dJydy[0] = 0.5*(-2*my921 + 2*y921)/pow(sigmay921, 2);
            break;
        case 922:
            dJydy[0] = 0.5*(-2*my922 + 2*y922)/pow(sigmay922, 2);
            break;
        case 923:
            dJydy[0] = 0.5*(-2*my923 + 2*y923)/pow(sigmay923, 2);
            break;
        case 924:
            dJydy[0] = 0.5*(-2*my924 + 2*y924)/pow(sigmay924, 2);
            break;
        case 925:
            dJydy[0] = 0.5*(-2*my925 + 2*y925)/pow(sigmay925, 2);
            break;
        case 926:
            dJydy[0] = 0.5*(-2*my926 + 2*y926)/pow(sigmay926, 2);
            break;
        case 927:
            dJydy[0] = 0.5*(-2*my927 + 2*y927)/pow(sigmay927, 2);
            break;
        case 928:
            dJydy[0] = 0.5*(-2*my928 + 2*y928)/pow(sigmay928, 2);
            break;
        case 929:
            dJydy[0] = 0.5*(-2*my929 + 2*y929)/pow(sigmay929, 2);
            break;
        case 930:
            dJydy[0] = 0.5*(-2*my930 + 2*y930)/pow(sigmay930, 2);
            break;
        case 931:
            dJydy[0] = 0.5*(-2*my931 + 2*y931)/pow(sigmay931, 2);
            break;
        case 932:
            dJydy[0] = 0.5*(-2*my932 + 2*y932)/pow(sigmay932, 2);
            break;
        case 933:
            dJydy[0] = 0.5*(-2*my933 + 2*y933)/pow(sigmay933, 2);
            break;
        case 934:
            dJydy[0] = 0.5*(-2*my934 + 2*y934)/pow(sigmay934, 2);
            break;
        case 935:
            dJydy[0] = 0.5*(-2*my935 + 2*y935)/pow(sigmay935, 2);
            break;
        case 936:
            dJydy[0] = 0.5*(-2*my936 + 2*y936)/pow(sigmay936, 2);
            break;
        case 937:
            dJydy[0] = 0.5*(-2*my937 + 2*y937)/pow(sigmay937, 2);
            break;
        case 938:
            dJydy[0] = 0.5*(-2*my938 + 2*y938)/pow(sigmay938, 2);
            break;
        case 939:
            dJydy[0] = 0.5*(-2*my939 + 2*y939)/pow(sigmay939, 2);
            break;
        case 940:
            dJydy[0] = 0.5*(-2*my940 + 2*y940)/pow(sigmay940, 2);
            break;
        case 941:
            dJydy[0] = 0.5*(-2*my941 + 2*y941)/pow(sigmay941, 2);
            break;
        case 942:
            dJydy[0] = 0.5*(-2*my942 + 2*y942)/pow(sigmay942, 2);
            break;
        case 943:
            dJydy[0] = 0.5*(-2*my943 + 2*y943)/pow(sigmay943, 2);
            break;
        case 944:
            dJydy[0] = 0.5*(-2*my944 + 2*y944)/pow(sigmay944, 2);
            break;
        case 945:
            dJydy[0] = 0.5*(-2*my945 + 2*y945)/pow(sigmay945, 2);
            break;
        case 946:
            dJydy[0] = 0.5*(-2*my946 + 2*y946)/pow(sigmay946, 2);
            break;
        case 947:
            dJydy[0] = 0.5*(-2*my947 + 2*y947)/pow(sigmay947, 2);
            break;
        case 948:
            dJydy[0] = 0.5*(-2*my948 + 2*y948)/pow(sigmay948, 2);
            break;
        case 949:
            dJydy[0] = 0.5*(-2*my949 + 2*y949)/pow(sigmay949, 2);
            break;
        case 950:
            dJydy[0] = 0.5*(-2*my950 + 2*y950)/pow(sigmay950, 2);
            break;
        case 951:
            dJydy[0] = 0.5*(-2*my951 + 2*y951)/pow(sigmay951, 2);
            break;
        case 952:
            dJydy[0] = 0.5*(-2*my952 + 2*y952)/pow(sigmay952, 2);
            break;
        case 953:
            dJydy[0] = 0.5*(-2*my953 + 2*y953)/pow(sigmay953, 2);
            break;
        case 954:
            dJydy[0] = 0.5*(-2*my954 + 2*y954)/pow(sigmay954, 2);
            break;
        case 955:
            dJydy[0] = 0.5*(-2*my955 + 2*y955)/pow(sigmay955, 2);
            break;
        case 956:
            dJydy[0] = 0.5*(-2*my956 + 2*y956)/pow(sigmay956, 2);
            break;
        case 957:
            dJydy[0] = 0.5*(-2*my957 + 2*y957)/pow(sigmay957, 2);
            break;
        case 958:
            dJydy[0] = 0.5*(-2*my958 + 2*y958)/pow(sigmay958, 2);
            break;
        case 959:
            dJydy[0] = 0.5*(-2*my959 + 2*y959)/pow(sigmay959, 2);
            break;
        case 960:
            dJydy[0] = 0.5*(-2*my960 + 2*y960)/pow(sigmay960, 2);
            break;
        case 961:
            dJydy[0] = 0.5*(-2*my961 + 2*y961)/pow(sigmay961, 2);
            break;
        case 962:
            dJydy[0] = 0.5*(-2*my962 + 2*y962)/pow(sigmay962, 2);
            break;
        case 963:
            dJydy[0] = 0.5*(-2*my963 + 2*y963)/pow(sigmay963, 2);
            break;
        case 964:
            dJydy[0] = 0.5*(-2*my964 + 2*y964)/pow(sigmay964, 2);
            break;
        case 965:
            dJydy[0] = 0.5*(-2*my965 + 2*y965)/pow(sigmay965, 2);
            break;
        case 966:
            dJydy[0] = 0.5*(-2*my966 + 2*y966)/pow(sigmay966, 2);
            break;
        case 967:
            dJydy[0] = 0.5*(-2*my967 + 2*y967)/pow(sigmay967, 2);
            break;
        case 968:
            dJydy[0] = 0.5*(-2*my968 + 2*y968)/pow(sigmay968, 2);
            break;
        case 969:
            dJydy[0] = 0.5*(-2*my969 + 2*y969)/pow(sigmay969, 2);
            break;
        case 970:
            dJydy[0] = 0.5*(-2*my970 + 2*y970)/pow(sigmay970, 2);
            break;
        case 971:
            dJydy[0] = 0.5*(-2*my971 + 2*y971)/pow(sigmay971, 2);
            break;
        case 972:
            dJydy[0] = 0.5*(-2*my972 + 2*y972)/pow(sigmay972, 2);
            break;
        case 973:
            dJydy[0] = 0.5*(-2*my973 + 2*y973)/pow(sigmay973, 2);
            break;
        case 974:
            dJydy[0] = 0.5*(-2*my974 + 2*y974)/pow(sigmay974, 2);
            break;
        case 975:
            dJydy[0] = 0.5*(-2*my975 + 2*y975)/pow(sigmay975, 2);
            break;
        case 976:
            dJydy[0] = 0.5*(-2*my976 + 2*y976)/pow(sigmay976, 2);
            break;
        case 977:
            dJydy[0] = 0.5*(-2*my977 + 2*y977)/pow(sigmay977, 2);
            break;
        case 978:
            dJydy[0] = 0.5*(-2*my978 + 2*y978)/pow(sigmay978, 2);
            break;
        case 979:
            dJydy[0] = 0.5*(-2*my979 + 2*y979)/pow(sigmay979, 2);
            break;
        case 980:
            dJydy[0] = 0.5*(-2*my980 + 2*y980)/pow(sigmay980, 2);
            break;
        case 981:
            dJydy[0] = 0.5*(-2*my981 + 2*y981)/pow(sigmay981, 2);
            break;
        case 982:
            dJydy[0] = 0.5*(-2*my982 + 2*y982)/pow(sigmay982, 2);
            break;
        case 983:
            dJydy[0] = 0.5*(-2*my983 + 2*y983)/pow(sigmay983, 2);
            break;
        case 984:
            dJydy[0] = 0.5*(-2*my984 + 2*y984)/pow(sigmay984, 2);
            break;
        case 985:
            dJydy[0] = 0.5*(-2*my985 + 2*y985)/pow(sigmay985, 2);
            break;
        case 986:
            dJydy[0] = 0.5*(-2*my986 + 2*y986)/pow(sigmay986, 2);
            break;
        case 987:
            dJydy[0] = 0.5*(-2*my987 + 2*y987)/pow(sigmay987, 2);
            break;
        case 988:
            dJydy[0] = 0.5*(-2*my988 + 2*y988)/pow(sigmay988, 2);
            break;
        case 989:
            dJydy[0] = 0.5*(-2*my989 + 2*y989)/pow(sigmay989, 2);
            break;
        case 990:
            dJydy[0] = 0.5*(-2*my990 + 2*y990)/pow(sigmay990, 2);
            break;
        case 991:
            dJydy[0] = 0.5*(-2*my991 + 2*y991)/pow(sigmay991, 2);
            break;
        case 992:
            dJydy[0] = 0.5*(-2*my992 + 2*y992)/pow(sigmay992, 2);
            break;
        case 993:
            dJydy[0] = 0.5*(-2*my993 + 2*y993)/pow(sigmay993, 2);
            break;
        case 994:
            dJydy[0] = 0.5*(-2*my994 + 2*y994)/pow(sigmay994, 2);
            break;
        case 995:
            dJydy[0] = 0.5*(-2*my995 + 2*y995)/pow(sigmay995, 2);
            break;
        case 996:
            dJydy[0] = 0.5*(-2*my996 + 2*y996)/pow(sigmay996, 2);
            break;
        case 997:
            dJydy[0] = 0.5*(-2*my997 + 2*y997)/pow(sigmay997, 2);
            break;
        case 998:
            dJydy[0] = 0.5*(-2*my998 + 2*y998)/pow(sigmay998, 2);
            break;
        case 999:
            dJydy[0] = 0.5*(-2*my999 + 2*y999)/pow(sigmay999, 2);
            break;
        case 1000:
            dJydy[0] = 0.5*(-2*my1000 + 2*y1000)/pow(sigmay1000, 2);
            break;
        case 1001:
            dJydy[0] = 0.5*(-2*my1001 + 2*y1001)/pow(sigmay1001, 2);
            break;
        case 1002:
            dJydy[0] = 0.5*(-2*my1002 + 2*y1002)/pow(sigmay1002, 2);
            break;
        case 1003:
            dJydy[0] = 0.5*(-2*my1003 + 2*y1003)/pow(sigmay1003, 2);
            break;
        case 1004:
            dJydy[0] = 0.5*(-2*my1004 + 2*y1004)/pow(sigmay1004, 2);
            break;
        case 1005:
            dJydy[0] = 0.5*(-2*my1005 + 2*y1005)/pow(sigmay1005, 2);
            break;
        case 1006:
            dJydy[0] = 0.5*(-2*my1006 + 2*y1006)/pow(sigmay1006, 2);
            break;
        case 1007:
            dJydy[0] = 0.5*(-2*my1007 + 2*y1007)/pow(sigmay1007, 2);
            break;
        case 1008:
            dJydy[0] = 0.5*(-2*my1008 + 2*y1008)/pow(sigmay1008, 2);
            break;
        case 1009:
            dJydy[0] = 0.5*(-2*my1009 + 2*y1009)/pow(sigmay1009, 2);
            break;
        case 1010:
            dJydy[0] = 0.5*(-2*my1010 + 2*y1010)/pow(sigmay1010, 2);
            break;
        case 1011:
            dJydy[0] = 0.5*(-2*my1011 + 2*y1011)/pow(sigmay1011, 2);
            break;
        case 1012:
            dJydy[0] = 0.5*(-2*my1012 + 2*y1012)/pow(sigmay1012, 2);
            break;
        case 1013:
            dJydy[0] = 0.5*(-2*my1013 + 2*y1013)/pow(sigmay1013, 2);
            break;
        case 1014:
            dJydy[0] = 0.5*(-2*my1014 + 2*y1014)/pow(sigmay1014, 2);
            break;
        case 1015:
            dJydy[0] = 0.5*(-2*my1015 + 2*y1015)/pow(sigmay1015, 2);
            break;
        case 1016:
            dJydy[0] = 0.5*(-2*my1016 + 2*y1016)/pow(sigmay1016, 2);
            break;
        case 1017:
            dJydy[0] = 0.5*(-2*my1017 + 2*y1017)/pow(sigmay1017, 2);
            break;
        case 1018:
            dJydy[0] = 0.5*(-2*my1018 + 2*y1018)/pow(sigmay1018, 2);
            break;
        case 1019:
            dJydy[0] = 0.5*(-2*my1019 + 2*y1019)/pow(sigmay1019, 2);
            break;
        case 1020:
            dJydy[0] = 0.5*(-2*my1020 + 2*y1020)/pow(sigmay1020, 2);
            break;
        case 1021:
            dJydy[0] = 0.5*(-2*my1021 + 2*y1021)/pow(sigmay1021, 2);
            break;
        case 1022:
            dJydy[0] = 0.5*(-2*my1022 + 2*y1022)/pow(sigmay1022, 2);
            break;
        case 1023:
            dJydy[0] = 0.5*(-2*my1023 + 2*y1023)/pow(sigmay1023, 2);
            break;
        case 1024:
            dJydy[0] = 0.5*(-2*my1024 + 2*y1024)/pow(sigmay1024, 2);
            break;
        case 1025:
            dJydy[0] = 0.5*(-2*my1025 + 2*y1025)/pow(sigmay1025, 2);
            break;
        case 1026:
            dJydy[0] = 0.5*(-2*my1026 + 2*y1026)/pow(sigmay1026, 2);
            break;
        case 1027:
            dJydy[0] = 0.5*(-2*my1027 + 2*y1027)/pow(sigmay1027, 2);
            break;
        case 1028:
            dJydy[0] = 0.5*(-2*my1028 + 2*y1028)/pow(sigmay1028, 2);
            break;
        case 1029:
            dJydy[0] = 0.5*(-2*my1029 + 2*y1029)/pow(sigmay1029, 2);
            break;
        case 1030:
            dJydy[0] = 0.5*(-2*my1030 + 2*y1030)/pow(sigmay1030, 2);
            break;
        case 1031:
            dJydy[0] = 0.5*(-2*my1031 + 2*y1031)/pow(sigmay1031, 2);
            break;
        case 1032:
            dJydy[0] = 0.5*(-2*my1032 + 2*y1032)/pow(sigmay1032, 2);
            break;
        case 1033:
            dJydy[0] = 0.5*(-2*my1033 + 2*y1033)/pow(sigmay1033, 2);
            break;
        case 1034:
            dJydy[0] = 0.5*(-2*my1034 + 2*y1034)/pow(sigmay1034, 2);
            break;
        case 1035:
            dJydy[0] = 0.5*(-2*my1035 + 2*y1035)/pow(sigmay1035, 2);
            break;
        case 1036:
            dJydy[0] = 0.5*(-2*my1036 + 2*y1036)/pow(sigmay1036, 2);
            break;
        case 1037:
            dJydy[0] = 0.5*(-2*my1037 + 2*y1037)/pow(sigmay1037, 2);
            break;
        case 1038:
            dJydy[0] = 0.5*(-2*my1038 + 2*y1038)/pow(sigmay1038, 2);
            break;
        case 1039:
            dJydy[0] = 0.5*(-2*my1039 + 2*y1039)/pow(sigmay1039, 2);
            break;
        case 1040:
            dJydy[0] = 0.5*(-2*my1040 + 2*y1040)/pow(sigmay1040, 2);
            break;
        case 1041:
            dJydy[0] = 0.5*(-2*my1041 + 2*y1041)/pow(sigmay1041, 2);
            break;
        case 1042:
            dJydy[0] = 0.5*(-2*my1042 + 2*y1042)/pow(sigmay1042, 2);
            break;
        case 1043:
            dJydy[0] = 0.5*(-2*my1043 + 2*y1043)/pow(sigmay1043, 2);
            break;
        case 1044:
            dJydy[0] = 0.5*(-2*my1044 + 2*y1044)/pow(sigmay1044, 2);
            break;
        case 1045:
            dJydy[0] = 0.5*(-2*my1045 + 2*y1045)/pow(sigmay1045, 2);
            break;
        case 1046:
            dJydy[0] = 0.5*(-2*my1046 + 2*y1046)/pow(sigmay1046, 2);
            break;
        case 1047:
            dJydy[0] = 0.5*(-2*my1047 + 2*y1047)/pow(sigmay1047, 2);
            break;
        case 1048:
            dJydy[0] = 0.5*(-2*my1048 + 2*y1048)/pow(sigmay1048, 2);
            break;
        case 1049:
            dJydy[0] = 0.5*(-2*my1049 + 2*y1049)/pow(sigmay1049, 2);
            break;
        case 1050:
            dJydy[0] = 0.5*(-2*my1050 + 2*y1050)/pow(sigmay1050, 2);
            break;
        case 1051:
            dJydy[0] = 0.5*(-2*my1051 + 2*y1051)/pow(sigmay1051, 2);
            break;
        case 1052:
            dJydy[0] = 0.5*(-2*my1052 + 2*y1052)/pow(sigmay1052, 2);
            break;
        case 1053:
            dJydy[0] = 0.5*(-2*my1053 + 2*y1053)/pow(sigmay1053, 2);
            break;
        case 1054:
            dJydy[0] = 0.5*(-2*my1054 + 2*y1054)/pow(sigmay1054, 2);
            break;
        case 1055:
            dJydy[0] = 0.5*(-2*my1055 + 2*y1055)/pow(sigmay1055, 2);
            break;
        case 1056:
            dJydy[0] = 0.5*(-2*my1056 + 2*y1056)/pow(sigmay1056, 2);
            break;
        case 1057:
            dJydy[0] = 0.5*(-2*my1057 + 2*y1057)/pow(sigmay1057, 2);
            break;
        case 1058:
            dJydy[0] = 0.5*(-2*my1058 + 2*y1058)/pow(sigmay1058, 2);
            break;
        case 1059:
            dJydy[0] = 0.5*(-2*my1059 + 2*y1059)/pow(sigmay1059, 2);
            break;
        case 1060:
            dJydy[0] = 0.5*(-2*my1060 + 2*y1060)/pow(sigmay1060, 2);
            break;
        case 1061:
            dJydy[0] = 0.5*(-2*my1061 + 2*y1061)/pow(sigmay1061, 2);
            break;
        case 1062:
            dJydy[0] = 0.5*(-2*my1062 + 2*y1062)/pow(sigmay1062, 2);
            break;
        case 1063:
            dJydy[0] = 0.5*(-2*my1063 + 2*y1063)/pow(sigmay1063, 2);
            break;
        case 1064:
            dJydy[0] = 0.5*(-2*my1064 + 2*y1064)/pow(sigmay1064, 2);
            break;
        case 1065:
            dJydy[0] = 0.5*(-2*my1065 + 2*y1065)/pow(sigmay1065, 2);
            break;
        case 1066:
            dJydy[0] = 0.5*(-2*my1066 + 2*y1066)/pow(sigmay1066, 2);
            break;
        case 1067:
            dJydy[0] = 0.5*(-2*my1067 + 2*y1067)/pow(sigmay1067, 2);
            break;
        case 1068:
            dJydy[0] = 0.5*(-2*my1068 + 2*y1068)/pow(sigmay1068, 2);
            break;
        case 1069:
            dJydy[0] = 0.5*(-2*my1069 + 2*y1069)/pow(sigmay1069, 2);
            break;
        case 1070:
            dJydy[0] = 0.5*(-2*my1070 + 2*y1070)/pow(sigmay1070, 2);
            break;
        case 1071:
            dJydy[0] = 0.5*(-2*my1071 + 2*y1071)/pow(sigmay1071, 2);
            break;
        case 1072:
            dJydy[0] = 0.5*(-2*my1072 + 2*y1072)/pow(sigmay1072, 2);
            break;
        case 1073:
            dJydy[0] = 0.5*(-2*my1073 + 2*y1073)/pow(sigmay1073, 2);
            break;
        case 1074:
            dJydy[0] = 0.5*(-2*my1074 + 2*y1074)/pow(sigmay1074, 2);
            break;
        case 1075:
            dJydy[0] = 0.5*(-2*my1075 + 2*y1075)/pow(sigmay1075, 2);
            break;
        case 1076:
            dJydy[0] = 0.5*(-2*my1076 + 2*y1076)/pow(sigmay1076, 2);
            break;
        case 1077:
            dJydy[0] = 0.5*(-2*my1077 + 2*y1077)/pow(sigmay1077, 2);
            break;
        case 1078:
            dJydy[0] = 0.5*(-2*my1078 + 2*y1078)/pow(sigmay1078, 2);
            break;
        case 1079:
            dJydy[0] = 0.5*(-2*my1079 + 2*y1079)/pow(sigmay1079, 2);
            break;
        case 1080:
            dJydy[0] = 0.5*(-2*my1080 + 2*y1080)/pow(sigmay1080, 2);
            break;
        case 1081:
            dJydy[0] = 0.5*(-2*my1081 + 2*y1081)/pow(sigmay1081, 2);
            break;
        case 1082:
            dJydy[0] = 0.5*(-2*my1082 + 2*y1082)/pow(sigmay1082, 2);
            break;
        case 1083:
            dJydy[0] = 0.5*(-2*my1083 + 2*y1083)/pow(sigmay1083, 2);
            break;
        case 1084:
            dJydy[0] = 0.5*(-2*my1084 + 2*y1084)/pow(sigmay1084, 2);
            break;
        case 1085:
            dJydy[0] = 0.5*(-2*my1085 + 2*y1085)/pow(sigmay1085, 2);
            break;
        case 1086:
            dJydy[0] = 0.5*(-2*my1086 + 2*y1086)/pow(sigmay1086, 2);
            break;
        case 1087:
            dJydy[0] = 0.5*(-2*my1087 + 2*y1087)/pow(sigmay1087, 2);
            break;
        case 1088:
            dJydy[0] = 0.5*(-2*my1088 + 2*y1088)/pow(sigmay1088, 2);
            break;
        case 1089:
            dJydy[0] = 0.5*(-2*my1089 + 2*y1089)/pow(sigmay1089, 2);
            break;
        case 1090:
            dJydy[0] = 0.5*(-2*my1090 + 2*y1090)/pow(sigmay1090, 2);
            break;
        case 1091:
            dJydy[0] = 0.5*(-2*my1091 + 2*y1091)/pow(sigmay1091, 2);
            break;
        case 1092:
            dJydy[0] = 0.5*(-2*my1092 + 2*y1092)/pow(sigmay1092, 2);
            break;
        case 1093:
            dJydy[0] = 0.5*(-2*my1093 + 2*y1093)/pow(sigmay1093, 2);
            break;
        case 1094:
            dJydy[0] = 0.5*(-2*my1094 + 2*y1094)/pow(sigmay1094, 2);
            break;
        case 1095:
            dJydy[0] = 0.5*(-2*my1095 + 2*y1095)/pow(sigmay1095, 2);
            break;
        case 1096:
            dJydy[0] = 0.5*(-2*my1096 + 2*y1096)/pow(sigmay1096, 2);
            break;
        case 1097:
            dJydy[0] = 0.5*(-2*my1097 + 2*y1097)/pow(sigmay1097, 2);
            break;
        case 1098:
            dJydy[0] = 0.5*(-2*my1098 + 2*y1098)/pow(sigmay1098, 2);
            break;
        case 1099:
            dJydy[0] = 0.5*(-2*my1099 + 2*y1099)/pow(sigmay1099, 2);
            break;
        case 1100:
            dJydy[0] = 0.5*(-2*my1100 + 2*y1100)/pow(sigmay1100, 2);
            break;
        case 1101:
            dJydy[0] = 0.5*(-2*my1101 + 2*y1101)/pow(sigmay1101, 2);
            break;
        case 1102:
            dJydy[0] = 0.5*(-2*my1102 + 2*y1102)/pow(sigmay1102, 2);
            break;
        case 1103:
            dJydy[0] = 0.5*(-2*my1103 + 2*y1103)/pow(sigmay1103, 2);
            break;
        case 1104:
            dJydy[0] = 0.5*(-2*my1104 + 2*y1104)/pow(sigmay1104, 2);
            break;
        case 1105:
            dJydy[0] = 0.5*(-2*my1105 + 2*y1105)/pow(sigmay1105, 2);
            break;
        case 1106:
            dJydy[0] = 0.5*(-2*my1106 + 2*y1106)/pow(sigmay1106, 2);
            break;
        case 1107:
            dJydy[0] = 0.5*(-2*my1107 + 2*y1107)/pow(sigmay1107, 2);
            break;
        case 1108:
            dJydy[0] = 0.5*(-2*my1108 + 2*y1108)/pow(sigmay1108, 2);
            break;
        case 1109:
            dJydy[0] = 0.5*(-2*my1109 + 2*y1109)/pow(sigmay1109, 2);
            break;
        case 1110:
            dJydy[0] = 0.5*(-2*my1110 + 2*y1110)/pow(sigmay1110, 2);
            break;
        case 1111:
            dJydy[0] = 0.5*(-2*my1111 + 2*y1111)/pow(sigmay1111, 2);
            break;
        case 1112:
            dJydy[0] = 0.5*(-2*my1112 + 2*y1112)/pow(sigmay1112, 2);
            break;
        case 1113:
            dJydy[0] = 0.5*(-2*my1113 + 2*y1113)/pow(sigmay1113, 2);
            break;
        case 1114:
            dJydy[0] = 0.5*(-2*my1114 + 2*y1114)/pow(sigmay1114, 2);
            break;
        case 1115:
            dJydy[0] = 0.5*(-2*my1115 + 2*y1115)/pow(sigmay1115, 2);
            break;
        case 1116:
            dJydy[0] = 0.5*(-2*my1116 + 2*y1116)/pow(sigmay1116, 2);
            break;
        case 1117:
            dJydy[0] = 0.5*(-2*my1117 + 2*y1117)/pow(sigmay1117, 2);
            break;
        case 1118:
            dJydy[0] = 0.5*(-2*my1118 + 2*y1118)/pow(sigmay1118, 2);
            break;
        case 1119:
            dJydy[0] = 0.5*(-2*my1119 + 2*y1119)/pow(sigmay1119, 2);
            break;
        case 1120:
            dJydy[0] = 0.5*(-2*my1120 + 2*y1120)/pow(sigmay1120, 2);
            break;
        case 1121:
            dJydy[0] = 0.5*(-2*my1121 + 2*y1121)/pow(sigmay1121, 2);
            break;
        case 1122:
            dJydy[0] = 0.5*(-2*my1122 + 2*y1122)/pow(sigmay1122, 2);
            break;
        case 1123:
            dJydy[0] = 0.5*(-2*my1123 + 2*y1123)/pow(sigmay1123, 2);
            break;
        case 1124:
            dJydy[0] = 0.5*(-2*my1124 + 2*y1124)/pow(sigmay1124, 2);
            break;
        case 1125:
            dJydy[0] = 0.5*(-2*my1125 + 2*y1125)/pow(sigmay1125, 2);
            break;
        case 1126:
            dJydy[0] = 0.5*(-2*my1126 + 2*y1126)/pow(sigmay1126, 2);
            break;
        case 1127:
            dJydy[0] = 0.5*(-2*my1127 + 2*y1127)/pow(sigmay1127, 2);
            break;
        case 1128:
            dJydy[0] = 0.5*(-2*my1128 + 2*y1128)/pow(sigmay1128, 2);
            break;
        case 1129:
            dJydy[0] = 0.5*(-2*my1129 + 2*y1129)/pow(sigmay1129, 2);
            break;
        case 1130:
            dJydy[0] = 0.5*(-2*my1130 + 2*y1130)/pow(sigmay1130, 2);
            break;
        case 1131:
            dJydy[0] = 0.5*(-2*my1131 + 2*y1131)/pow(sigmay1131, 2);
            break;
        case 1132:
            dJydy[0] = 0.5*(-2*my1132 + 2*y1132)/pow(sigmay1132, 2);
            break;
        case 1133:
            dJydy[0] = 0.5*(-2*my1133 + 2*y1133)/pow(sigmay1133, 2);
            break;
        case 1134:
            dJydy[0] = 0.5*(-2*my1134 + 2*y1134)/pow(sigmay1134, 2);
            break;
        case 1135:
            dJydy[0] = 0.5*(-2*my1135 + 2*y1135)/pow(sigmay1135, 2);
            break;
        case 1136:
            dJydy[0] = 0.5*(-2*my1136 + 2*y1136)/pow(sigmay1136, 2);
            break;
        case 1137:
            dJydy[0] = 0.5*(-2*my1137 + 2*y1137)/pow(sigmay1137, 2);
            break;
        case 1138:
            dJydy[0] = 0.5*(-2*my1138 + 2*y1138)/pow(sigmay1138, 2);
            break;
        case 1139:
            dJydy[0] = 0.5*(-2*my1139 + 2*y1139)/pow(sigmay1139, 2);
            break;
        case 1140:
            dJydy[0] = 0.5*(-2*my1140 + 2*y1140)/pow(sigmay1140, 2);
            break;
        case 1141:
            dJydy[0] = 0.5*(-2*my1141 + 2*y1141)/pow(sigmay1141, 2);
            break;
        case 1142:
            dJydy[0] = 0.5*(-2*my1142 + 2*y1142)/pow(sigmay1142, 2);
            break;
        case 1143:
            dJydy[0] = 0.5*(-2*my1143 + 2*y1143)/pow(sigmay1143, 2);
            break;
        case 1144:
            dJydy[0] = 0.5*(-2*my1144 + 2*y1144)/pow(sigmay1144, 2);
            break;
        case 1145:
            dJydy[0] = 0.5*(-2*my1145 + 2*y1145)/pow(sigmay1145, 2);
            break;
        case 1146:
            dJydy[0] = 0.5*(-2*my1146 + 2*y1146)/pow(sigmay1146, 2);
            break;
        case 1147:
            dJydy[0] = 0.5*(-2*my1147 + 2*y1147)/pow(sigmay1147, 2);
            break;
        case 1148:
            dJydy[0] = 0.5*(-2*my1148 + 2*y1148)/pow(sigmay1148, 2);
            break;
        case 1149:
            dJydy[0] = 0.5*(-2*my1149 + 2*y1149)/pow(sigmay1149, 2);
            break;
        case 1150:
            dJydy[0] = 0.5*(-2*my1150 + 2*y1150)/pow(sigmay1150, 2);
            break;
        case 1151:
            dJydy[0] = 0.5*(-2*my1151 + 2*y1151)/pow(sigmay1151, 2);
            break;
        case 1152:
            dJydy[0] = 0.5*(-2*my1152 + 2*y1152)/pow(sigmay1152, 2);
            break;
        case 1153:
            dJydy[0] = 0.5*(-2*my1153 + 2*y1153)/pow(sigmay1153, 2);
            break;
        case 1154:
            dJydy[0] = 0.5*(-2*my1154 + 2*y1154)/pow(sigmay1154, 2);
            break;
        case 1155:
            dJydy[0] = 0.5*(-2*my1155 + 2*y1155)/pow(sigmay1155, 2);
            break;
        case 1156:
            dJydy[0] = 0.5*(-2*my1156 + 2*y1156)/pow(sigmay1156, 2);
            break;
        case 1157:
            dJydy[0] = 0.5*(-2*my1157 + 2*y1157)/pow(sigmay1157, 2);
            break;
        case 1158:
            dJydy[0] = 0.5*(-2*my1158 + 2*y1158)/pow(sigmay1158, 2);
            break;
        case 1159:
            dJydy[0] = 0.5*(-2*my1159 + 2*y1159)/pow(sigmay1159, 2);
            break;
        case 1160:
            dJydy[0] = 0.5*(-2*my1160 + 2*y1160)/pow(sigmay1160, 2);
            break;
        case 1161:
            dJydy[0] = 0.5*(-2*my1161 + 2*y1161)/pow(sigmay1161, 2);
            break;
        case 1162:
            dJydy[0] = 0.5*(-2*my1162 + 2*y1162)/pow(sigmay1162, 2);
            break;
        case 1163:
            dJydy[0] = 0.5*(-2*my1163 + 2*y1163)/pow(sigmay1163, 2);
            break;
        case 1164:
            dJydy[0] = 0.5*(-2*my1164 + 2*y1164)/pow(sigmay1164, 2);
            break;
        case 1165:
            dJydy[0] = 0.5*(-2*my1165 + 2*y1165)/pow(sigmay1165, 2);
            break;
        case 1166:
            dJydy[0] = 0.5*(-2*my1166 + 2*y1166)/pow(sigmay1166, 2);
            break;
        case 1167:
            dJydy[0] = 0.5*(-2*my1167 + 2*y1167)/pow(sigmay1167, 2);
            break;
        case 1168:
            dJydy[0] = 0.5*(-2*my1168 + 2*y1168)/pow(sigmay1168, 2);
            break;
        case 1169:
            dJydy[0] = 0.5*(-2*my1169 + 2*y1169)/pow(sigmay1169, 2);
            break;
        case 1170:
            dJydy[0] = 0.5*(-2*my1170 + 2*y1170)/pow(sigmay1170, 2);
            break;
        case 1171:
            dJydy[0] = 0.5*(-2*my1171 + 2*y1171)/pow(sigmay1171, 2);
            break;
        case 1172:
            dJydy[0] = 0.5*(-2*my1172 + 2*y1172)/pow(sigmay1172, 2);
            break;
        case 1173:
            dJydy[0] = 0.5*(-2*my1173 + 2*y1173)/pow(sigmay1173, 2);
            break;
        case 1174:
            dJydy[0] = 0.5*(-2*my1174 + 2*y1174)/pow(sigmay1174, 2);
            break;
        case 1175:
            dJydy[0] = 0.5*(-2*my1175 + 2*y1175)/pow(sigmay1175, 2);
            break;
        case 1176:
            dJydy[0] = 0.5*(-2*my1176 + 2*y1176)/pow(sigmay1176, 2);
            break;
        case 1177:
            dJydy[0] = 0.5*(-2*my1177 + 2*y1177)/pow(sigmay1177, 2);
            break;
        case 1178:
            dJydy[0] = 0.5*(-2*my1178 + 2*y1178)/pow(sigmay1178, 2);
            break;
        case 1179:
            dJydy[0] = 0.5*(-2*my1179 + 2*y1179)/pow(sigmay1179, 2);
            break;
        case 1180:
            dJydy[0] = 0.5*(-2*my1180 + 2*y1180)/pow(sigmay1180, 2);
            break;
        case 1181:
            dJydy[0] = 0.5*(-2*my1181 + 2*y1181)/pow(sigmay1181, 2);
            break;
        case 1182:
            dJydy[0] = 0.5*(-2*my1182 + 2*y1182)/pow(sigmay1182, 2);
            break;
        case 1183:
            dJydy[0] = 0.5*(-2*my1183 + 2*y1183)/pow(sigmay1183, 2);
            break;
        case 1184:
            dJydy[0] = 0.5*(-2*my1184 + 2*y1184)/pow(sigmay1184, 2);
            break;
        case 1185:
            dJydy[0] = 0.5*(-2*my1185 + 2*y1185)/pow(sigmay1185, 2);
            break;
        case 1186:
            dJydy[0] = 0.5*(-2*my1186 + 2*y1186)/pow(sigmay1186, 2);
            break;
        case 1187:
            dJydy[0] = 0.5*(-2*my1187 + 2*y1187)/pow(sigmay1187, 2);
            break;
        case 1188:
            dJydy[0] = 0.5*(-2*my1188 + 2*y1188)/pow(sigmay1188, 2);
            break;
        case 1189:
            dJydy[0] = 0.5*(-2*my1189 + 2*y1189)/pow(sigmay1189, 2);
            break;
        case 1190:
            dJydy[0] = 0.5*(-2*my1190 + 2*y1190)/pow(sigmay1190, 2);
            break;
        case 1191:
            dJydy[0] = 0.5*(-2*my1191 + 2*y1191)/pow(sigmay1191, 2);
            break;
        case 1192:
            dJydy[0] = 0.5*(-2*my1192 + 2*y1192)/pow(sigmay1192, 2);
            break;
        case 1193:
            dJydy[0] = 0.5*(-2*my1193 + 2*y1193)/pow(sigmay1193, 2);
            break;
        case 1194:
            dJydy[0] = 0.5*(-2*my1194 + 2*y1194)/pow(sigmay1194, 2);
            break;
        case 1195:
            dJydy[0] = 0.5*(-2*my1195 + 2*y1195)/pow(sigmay1195, 2);
            break;
        case 1196:
            dJydy[0] = 0.5*(-2*my1196 + 2*y1196)/pow(sigmay1196, 2);
            break;
        case 1197:
            dJydy[0] = 0.5*(-2*my1197 + 2*y1197)/pow(sigmay1197, 2);
            break;
        case 1198:
            dJydy[0] = 0.5*(-2*my1198 + 2*y1198)/pow(sigmay1198, 2);
            break;
        case 1199:
            dJydy[0] = 0.5*(-2*my1199 + 2*y1199)/pow(sigmay1199, 2);
            break;
        case 1200:
            dJydy[0] = 0.5*(-2*my1200 + 2*y1200)/pow(sigmay1200, 2);
            break;
        case 1201:
            dJydy[0] = 0.5*(-2*my1201 + 2*y1201)/pow(sigmay1201, 2);
            break;
        case 1202:
            dJydy[0] = 0.5*(-2*my1202 + 2*y1202)/pow(sigmay1202, 2);
            break;
        case 1203:
            dJydy[0] = 0.5*(-2*my1203 + 2*y1203)/pow(sigmay1203, 2);
            break;
        case 1204:
            dJydy[0] = 0.5*(-2*my1204 + 2*y1204)/pow(sigmay1204, 2);
            break;
        case 1205:
            dJydy[0] = 0.5*(-2*my1205 + 2*y1205)/pow(sigmay1205, 2);
            break;
        case 1206:
            dJydy[0] = 0.5*(-2*my1206 + 2*y1206)/pow(sigmay1206, 2);
            break;
        case 1207:
            dJydy[0] = 0.5*(-2*my1207 + 2*y1207)/pow(sigmay1207, 2);
            break;
        case 1208:
            dJydy[0] = 0.5*(-2*my1208 + 2*y1208)/pow(sigmay1208, 2);
            break;
        case 1209:
            dJydy[0] = 0.5*(-2*my1209 + 2*y1209)/pow(sigmay1209, 2);
            break;
        case 1210:
            dJydy[0] = 0.5*(-2*my1210 + 2*y1210)/pow(sigmay1210, 2);
            break;
        case 1211:
            dJydy[0] = 0.5*(-2*my1211 + 2*y1211)/pow(sigmay1211, 2);
            break;
        case 1212:
            dJydy[0] = 0.5*(-2*my1212 + 2*y1212)/pow(sigmay1212, 2);
            break;
        case 1213:
            dJydy[0] = 0.5*(-2*my1213 + 2*y1213)/pow(sigmay1213, 2);
            break;
        case 1214:
            dJydy[0] = 0.5*(-2*my1214 + 2*y1214)/pow(sigmay1214, 2);
            break;
        case 1215:
            dJydy[0] = 0.5*(-2*my1215 + 2*y1215)/pow(sigmay1215, 2);
            break;
        case 1216:
            dJydy[0] = 0.5*(-2*my1216 + 2*y1216)/pow(sigmay1216, 2);
            break;
        case 1217:
            dJydy[0] = 0.5*(-2*my1217 + 2*y1217)/pow(sigmay1217, 2);
            break;
        case 1218:
            dJydy[0] = 0.5*(-2*my1218 + 2*y1218)/pow(sigmay1218, 2);
            break;
        case 1219:
            dJydy[0] = 0.5*(-2*my1219 + 2*y1219)/pow(sigmay1219, 2);
            break;
        case 1220:
            dJydy[0] = 0.5*(-2*my1220 + 2*y1220)/pow(sigmay1220, 2);
            break;
        case 1221:
            dJydy[0] = 0.5*(-2*my1221 + 2*y1221)/pow(sigmay1221, 2);
            break;
        case 1222:
            dJydy[0] = 0.5*(-2*my1222 + 2*y1222)/pow(sigmay1222, 2);
            break;
        case 1223:
            dJydy[0] = 0.5*(-2*my1223 + 2*y1223)/pow(sigmay1223, 2);
            break;
        case 1224:
            dJydy[0] = 0.5*(-2*my1224 + 2*y1224)/pow(sigmay1224, 2);
            break;
        case 1225:
            dJydy[0] = 0.5*(-2*my1225 + 2*y1225)/pow(sigmay1225, 2);
            break;
        case 1226:
            dJydy[0] = 0.5*(-2*my1226 + 2*y1226)/pow(sigmay1226, 2);
            break;
        case 1227:
            dJydy[0] = 0.5*(-2*my1227 + 2*y1227)/pow(sigmay1227, 2);
            break;
        case 1228:
            dJydy[0] = 0.5*(-2*my1228 + 2*y1228)/pow(sigmay1228, 2);
            break;
        case 1229:
            dJydy[0] = 0.5*(-2*my1229 + 2*y1229)/pow(sigmay1229, 2);
            break;
        case 1230:
            dJydy[0] = 0.5*(-2*my1230 + 2*y1230)/pow(sigmay1230, 2);
            break;
        case 1231:
            dJydy[0] = 0.5*(-2*my1231 + 2*y1231)/pow(sigmay1231, 2);
            break;
        case 1232:
            dJydy[0] = 0.5*(-2*my1232 + 2*y1232)/pow(sigmay1232, 2);
            break;
        case 1233:
            dJydy[0] = 0.5*(-2*my1233 + 2*y1233)/pow(sigmay1233, 2);
            break;
        case 1234:
            dJydy[0] = 0.5*(-2*my1234 + 2*y1234)/pow(sigmay1234, 2);
            break;
        case 1235:
            dJydy[0] = 0.5*(-2*my1235 + 2*y1235)/pow(sigmay1235, 2);
            break;
        case 1236:
            dJydy[0] = 0.5*(-2*my1236 + 2*y1236)/pow(sigmay1236, 2);
            break;
        case 1237:
            dJydy[0] = 0.5*(-2*my1237 + 2*y1237)/pow(sigmay1237, 2);
            break;
        case 1238:
            dJydy[0] = 0.5*(-2*my1238 + 2*y1238)/pow(sigmay1238, 2);
            break;
        case 1239:
            dJydy[0] = 0.5*(-2*my1239 + 2*y1239)/pow(sigmay1239, 2);
            break;
        case 1240:
            dJydy[0] = 0.5*(-2*my1240 + 2*y1240)/pow(sigmay1240, 2);
            break;
        case 1241:
            dJydy[0] = 0.5*(-2*my1241 + 2*y1241)/pow(sigmay1241, 2);
            break;
        case 1242:
            dJydy[0] = 0.5*(-2*my1242 + 2*y1242)/pow(sigmay1242, 2);
            break;
        case 1243:
            dJydy[0] = 0.5*(-2*my1243 + 2*y1243)/pow(sigmay1243, 2);
            break;
        case 1244:
            dJydy[0] = 0.5*(-2*my1244 + 2*y1244)/pow(sigmay1244, 2);
            break;
        case 1245:
            dJydy[0] = 0.5*(-2*my1245 + 2*y1245)/pow(sigmay1245, 2);
            break;
        case 1246:
            dJydy[0] = 0.5*(-2*my1246 + 2*y1246)/pow(sigmay1246, 2);
            break;
        case 1247:
            dJydy[0] = 0.5*(-2*my1247 + 2*y1247)/pow(sigmay1247, 2);
            break;
        case 1248:
            dJydy[0] = 0.5*(-2*my1248 + 2*y1248)/pow(sigmay1248, 2);
            break;
        case 1249:
            dJydy[0] = 0.5*(-2*my1249 + 2*y1249)/pow(sigmay1249, 2);
            break;
        case 1250:
            dJydy[0] = 0.5*(-2*my1250 + 2*y1250)/pow(sigmay1250, 2);
            break;
        case 1251:
            dJydy[0] = 0.5*(-2*my1251 + 2*y1251)/pow(sigmay1251, 2);
            break;
        case 1252:
            dJydy[0] = 0.5*(-2*my1252 + 2*y1252)/pow(sigmay1252, 2);
            break;
        case 1253:
            dJydy[0] = 0.5*(-2*my1253 + 2*y1253)/pow(sigmay1253, 2);
            break;
        case 1254:
            dJydy[0] = 0.5*(-2*my1254 + 2*y1254)/pow(sigmay1254, 2);
            break;
        case 1255:
            dJydy[0] = 0.5*(-2*my1255 + 2*y1255)/pow(sigmay1255, 2);
            break;
        case 1256:
            dJydy[0] = 0.5*(-2*my1256 + 2*y1256)/pow(sigmay1256, 2);
            break;
        case 1257:
            dJydy[0] = 0.5*(-2*my1257 + 2*y1257)/pow(sigmay1257, 2);
            break;
        case 1258:
            dJydy[0] = 0.5*(-2*my1258 + 2*y1258)/pow(sigmay1258, 2);
            break;
        case 1259:
            dJydy[0] = 0.5*(-2*my1259 + 2*y1259)/pow(sigmay1259, 2);
            break;
        case 1260:
            dJydy[0] = 0.5*(-2*my1260 + 2*y1260)/pow(sigmay1260, 2);
            break;
        case 1261:
            dJydy[0] = 0.5*(-2*my1261 + 2*y1261)/pow(sigmay1261, 2);
            break;
        case 1262:
            dJydy[0] = 0.5*(-2*my1262 + 2*y1262)/pow(sigmay1262, 2);
            break;
        case 1263:
            dJydy[0] = 0.5*(-2*my1263 + 2*y1263)/pow(sigmay1263, 2);
            break;
        case 1264:
            dJydy[0] = 0.5*(-2*my1264 + 2*y1264)/pow(sigmay1264, 2);
            break;
        case 1265:
            dJydy[0] = 0.5*(-2*my1265 + 2*y1265)/pow(sigmay1265, 2);
            break;
        case 1266:
            dJydy[0] = 0.5*(-2*my1266 + 2*y1266)/pow(sigmay1266, 2);
            break;
        case 1267:
            dJydy[0] = 0.5*(-2*my1267 + 2*y1267)/pow(sigmay1267, 2);
            break;
        case 1268:
            dJydy[0] = 0.5*(-2*my1268 + 2*y1268)/pow(sigmay1268, 2);
            break;
        case 1269:
            dJydy[0] = 0.5*(-2*my1269 + 2*y1269)/pow(sigmay1269, 2);
            break;
        case 1270:
            dJydy[0] = 0.5*(-2*my1270 + 2*y1270)/pow(sigmay1270, 2);
            break;
        case 1271:
            dJydy[0] = 0.5*(-2*my1271 + 2*y1271)/pow(sigmay1271, 2);
            break;
        case 1272:
            dJydy[0] = 0.5*(-2*my1272 + 2*y1272)/pow(sigmay1272, 2);
            break;
        case 1273:
            dJydy[0] = 0.5*(-2*my1273 + 2*y1273)/pow(sigmay1273, 2);
            break;
        case 1274:
            dJydy[0] = 0.5*(-2*my1274 + 2*y1274)/pow(sigmay1274, 2);
            break;
        case 1275:
            dJydy[0] = 0.5*(-2*my1275 + 2*y1275)/pow(sigmay1275, 2);
            break;
        case 1276:
            dJydy[0] = 0.5*(-2*my1276 + 2*y1276)/pow(sigmay1276, 2);
            break;
        case 1277:
            dJydy[0] = 0.5*(-2*my1277 + 2*y1277)/pow(sigmay1277, 2);
            break;
        case 1278:
            dJydy[0] = 0.5*(-2*my1278 + 2*y1278)/pow(sigmay1278, 2);
            break;
        case 1279:
            dJydy[0] = 0.5*(-2*my1279 + 2*y1279)/pow(sigmay1279, 2);
            break;
        case 1280:
            dJydy[0] = 0.5*(-2*my1280 + 2*y1280)/pow(sigmay1280, 2);
            break;
        case 1281:
            dJydy[0] = 0.5*(-2*my1281 + 2*y1281)/pow(sigmay1281, 2);
            break;
        case 1282:
            dJydy[0] = 0.5*(-2*my1282 + 2*y1282)/pow(sigmay1282, 2);
            break;
        case 1283:
            dJydy[0] = 0.5*(-2*my1283 + 2*y1283)/pow(sigmay1283, 2);
            break;
        case 1284:
            dJydy[0] = 0.5*(-2*my1284 + 2*y1284)/pow(sigmay1284, 2);
            break;
        case 1285:
            dJydy[0] = 0.5*(-2*my1285 + 2*y1285)/pow(sigmay1285, 2);
            break;
        case 1286:
            dJydy[0] = 0.5*(-2*my1286 + 2*y1286)/pow(sigmay1286, 2);
            break;
        case 1287:
            dJydy[0] = 0.5*(-2*my1287 + 2*y1287)/pow(sigmay1287, 2);
            break;
        case 1288:
            dJydy[0] = 0.5*(-2*my1288 + 2*y1288)/pow(sigmay1288, 2);
            break;
        case 1289:
            dJydy[0] = 0.5*(-2*my1289 + 2*y1289)/pow(sigmay1289, 2);
            break;
        case 1290:
            dJydy[0] = 0.5*(-2*my1290 + 2*y1290)/pow(sigmay1290, 2);
            break;
        case 1291:
            dJydy[0] = 0.5*(-2*my1291 + 2*y1291)/pow(sigmay1291, 2);
            break;
        case 1292:
            dJydy[0] = 0.5*(-2*my1292 + 2*y1292)/pow(sigmay1292, 2);
            break;
        case 1293:
            dJydy[0] = 0.5*(-2*my1293 + 2*y1293)/pow(sigmay1293, 2);
            break;
        case 1294:
            dJydy[0] = 0.5*(-2*my1294 + 2*y1294)/pow(sigmay1294, 2);
            break;
        case 1295:
            dJydy[0] = 0.5*(-2*my1295 + 2*y1295)/pow(sigmay1295, 2);
            break;
        case 1296:
            dJydy[0] = 0.5*(-2*my1296 + 2*y1296)/pow(sigmay1296, 2);
            break;
        case 1297:
            dJydy[0] = 0.5*(-2*my1297 + 2*y1297)/pow(sigmay1297, 2);
            break;
        case 1298:
            dJydy[0] = 0.5*(-2*my1298 + 2*y1298)/pow(sigmay1298, 2);
            break;
        case 1299:
            dJydy[0] = 0.5*(-2*my1299 + 2*y1299)/pow(sigmay1299, 2);
            break;
        case 1300:
            dJydy[0] = 0.5*(-2*my1300 + 2*y1300)/pow(sigmay1300, 2);
            break;
        case 1301:
            dJydy[0] = 0.5*(-2*my1301 + 2*y1301)/pow(sigmay1301, 2);
            break;
        case 1302:
            dJydy[0] = 0.5*(-2*my1302 + 2*y1302)/pow(sigmay1302, 2);
            break;
        case 1303:
            dJydy[0] = 0.5*(-2*my1303 + 2*y1303)/pow(sigmay1303, 2);
            break;
        case 1304:
            dJydy[0] = 0.5*(-2*my1304 + 2*y1304)/pow(sigmay1304, 2);
            break;
        case 1305:
            dJydy[0] = 0.5*(-2*my1305 + 2*y1305)/pow(sigmay1305, 2);
            break;
        case 1306:
            dJydy[0] = 0.5*(-2*my1306 + 2*y1306)/pow(sigmay1306, 2);
            break;
        case 1307:
            dJydy[0] = 0.5*(-2*my1307 + 2*y1307)/pow(sigmay1307, 2);
            break;
        case 1308:
            dJydy[0] = 0.5*(-2*my1308 + 2*y1308)/pow(sigmay1308, 2);
            break;
        case 1309:
            dJydy[0] = 0.5*(-2*my1309 + 2*y1309)/pow(sigmay1309, 2);
            break;
        case 1310:
            dJydy[0] = 0.5*(-2*my1310 + 2*y1310)/pow(sigmay1310, 2);
            break;
        case 1311:
            dJydy[0] = 0.5*(-2*my1311 + 2*y1311)/pow(sigmay1311, 2);
            break;
        case 1312:
            dJydy[0] = 0.5*(-2*my1312 + 2*y1312)/pow(sigmay1312, 2);
            break;
        case 1313:
            dJydy[0] = 0.5*(-2*my1313 + 2*y1313)/pow(sigmay1313, 2);
            break;
        case 1314:
            dJydy[0] = 0.5*(-2*my1314 + 2*y1314)/pow(sigmay1314, 2);
            break;
        case 1315:
            dJydy[0] = 0.5*(-2*my1315 + 2*y1315)/pow(sigmay1315, 2);
            break;
        case 1316:
            dJydy[0] = 0.5*(-2*my1316 + 2*y1316)/pow(sigmay1316, 2);
            break;
        case 1317:
            dJydy[0] = 0.5*(-2*my1317 + 2*y1317)/pow(sigmay1317, 2);
            break;
        case 1318:
            dJydy[0] = 0.5*(-2*my1318 + 2*y1318)/pow(sigmay1318, 2);
            break;
        case 1319:
            dJydy[0] = 0.5*(-2*my1319 + 2*y1319)/pow(sigmay1319, 2);
            break;
        case 1320:
            dJydy[0] = 0.5*(-2*my1320 + 2*y1320)/pow(sigmay1320, 2);
            break;
        case 1321:
            dJydy[0] = 0.5*(-2*my1321 + 2*y1321)/pow(sigmay1321, 2);
            break;
        case 1322:
            dJydy[0] = 0.5*(-2*my1322 + 2*y1322)/pow(sigmay1322, 2);
            break;
        case 1323:
            dJydy[0] = 0.5*(-2*my1323 + 2*y1323)/pow(sigmay1323, 2);
            break;
        case 1324:
            dJydy[0] = 0.5*(-2*my1324 + 2*y1324)/pow(sigmay1324, 2);
            break;
        case 1325:
            dJydy[0] = 0.5*(-2*my1325 + 2*y1325)/pow(sigmay1325, 2);
            break;
        case 1326:
            dJydy[0] = 0.5*(-2*my1326 + 2*y1326)/pow(sigmay1326, 2);
            break;
        case 1327:
            dJydy[0] = 0.5*(-2*my1327 + 2*y1327)/pow(sigmay1327, 2);
            break;
        case 1328:
            dJydy[0] = 0.5*(-2*my1328 + 2*y1328)/pow(sigmay1328, 2);
            break;
        case 1329:
            dJydy[0] = 0.5*(-2*my1329 + 2*y1329)/pow(sigmay1329, 2);
            break;
        case 1330:
            dJydy[0] = 0.5*(-2*my1330 + 2*y1330)/pow(sigmay1330, 2);
            break;
        case 1331:
            dJydy[0] = 0.5*(-2*my1331 + 2*y1331)/pow(sigmay1331, 2);
            break;
        case 1332:
            dJydy[0] = 0.5*(-2*my1332 + 2*y1332)/pow(sigmay1332, 2);
            break;
        case 1333:
            dJydy[0] = 0.5*(-2*my1333 + 2*y1333)/pow(sigmay1333, 2);
            break;
        case 1334:
            dJydy[0] = 0.5*(-2*my1334 + 2*y1334)/pow(sigmay1334, 2);
            break;
        case 1335:
            dJydy[0] = 0.5*(-2*my1335 + 2*y1335)/pow(sigmay1335, 2);
            break;
        case 1336:
            dJydy[0] = 0.5*(-2*my1336 + 2*y1336)/pow(sigmay1336, 2);
            break;
        case 1337:
            dJydy[0] = 0.5*(-2*my1337 + 2*y1337)/pow(sigmay1337, 2);
            break;
        case 1338:
            dJydy[0] = 0.5*(-2*my1338 + 2*y1338)/pow(sigmay1338, 2);
            break;
        case 1339:
            dJydy[0] = 0.5*(-2*my1339 + 2*y1339)/pow(sigmay1339, 2);
            break;
        case 1340:
            dJydy[0] = 0.5*(-2*my1340 + 2*y1340)/pow(sigmay1340, 2);
            break;
        case 1341:
            dJydy[0] = 0.5*(-2*my1341 + 2*y1341)/pow(sigmay1341, 2);
            break;
        case 1342:
            dJydy[0] = 0.5*(-2*my1342 + 2*y1342)/pow(sigmay1342, 2);
            break;
        case 1343:
            dJydy[0] = 0.5*(-2*my1343 + 2*y1343)/pow(sigmay1343, 2);
            break;
        case 1344:
            dJydy[0] = 0.5*(-2*my1344 + 2*y1344)/pow(sigmay1344, 2);
            break;
        case 1345:
            dJydy[0] = 0.5*(-2*my1345 + 2*y1345)/pow(sigmay1345, 2);
            break;
        case 1346:
            dJydy[0] = 0.5*(-2*my1346 + 2*y1346)/pow(sigmay1346, 2);
            break;
        case 1347:
            dJydy[0] = 0.5*(-2*my1347 + 2*y1347)/pow(sigmay1347, 2);
            break;
        case 1348:
            dJydy[0] = 0.5*(-2*my1348 + 2*y1348)/pow(sigmay1348, 2);
            break;
        case 1349:
            dJydy[0] = 0.5*(-2*my1349 + 2*y1349)/pow(sigmay1349, 2);
            break;
        case 1350:
            dJydy[0] = 0.5*(-2*my1350 + 2*y1350)/pow(sigmay1350, 2);
            break;
        case 1351:
            dJydy[0] = 0.5*(-2*my1351 + 2*y1351)/pow(sigmay1351, 2);
            break;
        case 1352:
            dJydy[0] = 0.5*(-2*my1352 + 2*y1352)/pow(sigmay1352, 2);
            break;
        case 1353:
            dJydy[0] = 0.5*(-2*my1353 + 2*y1353)/pow(sigmay1353, 2);
            break;
        case 1354:
            dJydy[0] = 0.5*(-2*my1354 + 2*y1354)/pow(sigmay1354, 2);
            break;
        case 1355:
            dJydy[0] = 0.5*(-2*my1355 + 2*y1355)/pow(sigmay1355, 2);
            break;
        case 1356:
            dJydy[0] = 0.5*(-2*my1356 + 2*y1356)/pow(sigmay1356, 2);
            break;
        case 1357:
            dJydy[0] = 0.5*(-2*my1357 + 2*y1357)/pow(sigmay1357, 2);
            break;
        case 1358:
            dJydy[0] = 0.5*(-2*my1358 + 2*y1358)/pow(sigmay1358, 2);
            break;
        case 1359:
            dJydy[0] = 0.5*(-2*my1359 + 2*y1359)/pow(sigmay1359, 2);
            break;
        case 1360:
            dJydy[0] = 0.5*(-2*my1360 + 2*y1360)/pow(sigmay1360, 2);
            break;
        case 1361:
            dJydy[0] = 0.5*(-2*my1361 + 2*y1361)/pow(sigmay1361, 2);
            break;
        case 1362:
            dJydy[0] = 0.5*(-2*my1362 + 2*y1362)/pow(sigmay1362, 2);
            break;
        case 1363:
            dJydy[0] = 0.5*(-2*my1363 + 2*y1363)/pow(sigmay1363, 2);
            break;
        case 1364:
            dJydy[0] = 0.5*(-2*my1364 + 2*y1364)/pow(sigmay1364, 2);
            break;
        case 1365:
            dJydy[0] = 0.5*(-2*my1365 + 2*y1365)/pow(sigmay1365, 2);
            break;
        case 1366:
            dJydy[0] = 0.5*(-2*my1366 + 2*y1366)/pow(sigmay1366, 2);
            break;
        case 1367:
            dJydy[0] = 0.5*(-2*my1367 + 2*y1367)/pow(sigmay1367, 2);
            break;
        case 1368:
            dJydy[0] = 0.5*(-2*my1368 + 2*y1368)/pow(sigmay1368, 2);
            break;
        case 1369:
            dJydy[0] = 0.5*(-2*my1369 + 2*y1369)/pow(sigmay1369, 2);
            break;
        case 1370:
            dJydy[0] = 0.5*(-2*my1370 + 2*y1370)/pow(sigmay1370, 2);
            break;
        case 1371:
            dJydy[0] = 0.5*(-2*my1371 + 2*y1371)/pow(sigmay1371, 2);
            break;
        case 1372:
            dJydy[0] = 0.5*(-2*my1372 + 2*y1372)/pow(sigmay1372, 2);
            break;
        case 1373:
            dJydy[0] = 0.5*(-2*my1373 + 2*y1373)/pow(sigmay1373, 2);
            break;
        case 1374:
            dJydy[0] = 0.5*(-2*my1374 + 2*y1374)/pow(sigmay1374, 2);
            break;
        case 1375:
            dJydy[0] = 0.5*(-2*my1375 + 2*y1375)/pow(sigmay1375, 2);
            break;
        case 1376:
            dJydy[0] = 0.5*(-2*my1376 + 2*y1376)/pow(sigmay1376, 2);
            break;
        case 1377:
            dJydy[0] = 0.5*(-2*my1377 + 2*y1377)/pow(sigmay1377, 2);
            break;
        case 1378:
            dJydy[0] = 0.5*(-2*my1378 + 2*y1378)/pow(sigmay1378, 2);
            break;
        case 1379:
            dJydy[0] = 0.5*(-2*my1379 + 2*y1379)/pow(sigmay1379, 2);
            break;
        case 1380:
            dJydy[0] = 0.5*(-2*my1380 + 2*y1380)/pow(sigmay1380, 2);
            break;
        case 1381:
            dJydy[0] = 0.5*(-2*my1381 + 2*y1381)/pow(sigmay1381, 2);
            break;
        case 1382:
            dJydy[0] = 0.5*(-2*my1382 + 2*y1382)/pow(sigmay1382, 2);
            break;
        case 1383:
            dJydy[0] = 0.5*(-2*my1383 + 2*y1383)/pow(sigmay1383, 2);
            break;
        case 1384:
            dJydy[0] = 0.5*(-2*my1384 + 2*y1384)/pow(sigmay1384, 2);
            break;
        case 1385:
            dJydy[0] = 0.5*(-2*my1385 + 2*y1385)/pow(sigmay1385, 2);
            break;
        case 1386:
            dJydy[0] = 0.5*(-2*my1386 + 2*y1386)/pow(sigmay1386, 2);
            break;
        case 1387:
            dJydy[0] = 0.5*(-2*my1387 + 2*y1387)/pow(sigmay1387, 2);
            break;
        case 1388:
            dJydy[0] = 0.5*(-2*my1388 + 2*y1388)/pow(sigmay1388, 2);
            break;
        case 1389:
            dJydy[0] = 0.5*(-2*my1389 + 2*y1389)/pow(sigmay1389, 2);
            break;
        case 1390:
            dJydy[0] = 0.5*(-2*my1390 + 2*y1390)/pow(sigmay1390, 2);
            break;
        case 1391:
            dJydy[0] = 0.5*(-2*my1391 + 2*y1391)/pow(sigmay1391, 2);
            break;
        case 1392:
            dJydy[0] = 0.5*(-2*my1392 + 2*y1392)/pow(sigmay1392, 2);
            break;
        case 1393:
            dJydy[0] = 0.5*(-2*my1393 + 2*y1393)/pow(sigmay1393, 2);
            break;
        case 1394:
            dJydy[0] = 0.5*(-2*my1394 + 2*y1394)/pow(sigmay1394, 2);
            break;
        case 1395:
            dJydy[0] = 0.5*(-2*my1395 + 2*y1395)/pow(sigmay1395, 2);
            break;
    }
}