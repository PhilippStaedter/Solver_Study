#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Ouzounoglou2014(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[115] = 1.0*s27;
            dwdp[116] = 1.0*s26;
            dwdp[117] = 1.0*s25;
            dwdp[118] = 1.0*s21;
            dwdp[119] = 1.0*s2;
            dwdp[120] = 1.0*s1;
            dwdp[121] = 1.0*s5;
            dwdp[122] = 1.0*s6;
            dwdp[123] = 1.0*s29;
            dwdp[124] = 1.0*s30;
            dwdp[125] = 1.0*s31;
            dwdp[126] = 1.0*s32;
            dwdp[127] = 1.0*s23;
            dwdp[128] = 1.0*s24;
            dwdp[129] = 1.0*s20;
            dwdp[130] = 1.0*s18;
            break;
        case 1:
            dwdp[2] = 1.0*pow(s17, 2);
            dwdp[31] = 1.0*pow(s7, 2);
            dwdp[32] = 1.0*s536*s7;
            dwdp[92] = 1.0*s17*s78;
            break;
        case 2:
            dwdp[84] = 1.0*s51*s7;
            break;
        case 3:
            dwdp[69] = 1.0*s501;
            dwdp[74] = 1.0*s494;
            dwdp[75] = 1.0*s495;
            dwdp[76] = 1.0*s496;
            dwdp[77] = 1.0*s498;
            dwdp[78] = 1.0*s499;
            dwdp[79] = 1.0*s500;
            break;
        case 4:
            dwdp[5] = 1.0*s18;
            dwdp[8] = 1.0*s20;
            dwdp[11] = 1.0*s24;
            dwdp[14] = 1.0*s23;
            dwdp[18] = 1.0*s32;
            dwdp[21] = 1.0*s31;
            dwdp[25] = 1.0*s30;
            dwdp[46] = 1.0*s6;
            dwdp[48] = 1.0*s5;
            dwdp[51] = 1.0*s2;
            dwdp[54] = 1.0*s1;
            dwdp[57] = 1.0*s21;
            dwdp[60] = 1.0*s25;
            dwdp[63] = 1.0*s26;
            dwdp[132] = 1.0*s17;
            dwdp[134] = 1.0*s7;
            break;
        case 5:
            dwdp[6] = 1.0*s17*s18;
            dwdp[9] = 1.0*s17*s20;
            dwdp[12] = 1.0*s17*s24;
            dwdp[15] = 1.0*s17*s23;
            dwdp[19] = 1.0*s17*s32;
            dwdp[22] = 1.0*s17*s31;
            dwdp[26] = 1.0*s17*s30;
            dwdp[33] = 1.0*s490*s7;
            dwdp[34] = 1.0*s489*s7;
            dwdp[35] = 1.0*s492*s7;
            dwdp[45] = 1.0*s17*s29;
            dwdp[47] = 1.0*s6*s7;
            dwdp[49] = 1.0*s5*s7;
            dwdp[52] = 1.0*s2*s7;
            dwdp[55] = 1.0*s1*s7;
            dwdp[58] = 1.0*s21*s7;
            dwdp[61] = 1.0*s25*s7;
            dwdp[64] = 1.0*s26*s7;
            dwdp[70] = 1.0*s482*s7;
            dwdp[71] = 1.0*s483*s7;
            dwdp[72] = 1.0*s484*s7;
            dwdp[73] = 1.0*s491*s7;
            dwdp[93] = 1.0*s17*s85;
            dwdp[94] = 1.0*s17*s494;
            dwdp[95] = 1.0*s17*s495;
            dwdp[96] = 1.0*s17*s496;
            dwdp[97] = 1.0*s17*s498;
            dwdp[98] = 1.0*s17*s499;
            dwdp[99] = 1.0*s17*s500;
            break;
        case 6:
            dwdp[10] = 1.0*s20*s35;
            dwdp[13] = 1.0*s24*s35;
            dwdp[16] = 1.0*s23*s35;
            dwdp[20] = 1.0*s32*s35;
            dwdp[23] = 1.0*s31*s35;
            dwdp[27] = 1.0*s30*s35;
            dwdp[28] = 1.0*s29*s35;
            dwdp[50] = 1.0*s35*s5;
            dwdp[53] = 1.0*s2*s35;
            dwdp[56] = 1.0*s1*s35;
            dwdp[59] = 1.0*s21*s35;
            dwdp[62] = 1.0*s25*s35;
            dwdp[65] = 1.0*s26*s35;
            dwdp[66] = 1.0*s27*s35;
            dwdp[114] = 1.0*s33*s35;
            break;
        case 7:
            dwdp[100] = 1.0*s381;
            dwdp[101] = 1.0*s383;
            dwdp[102] = 1.0*s385;
            dwdp[103] = 1.0*s387;
            dwdp[104] = 1.0*s389;
            dwdp[105] = 1.0*s391;
            dwdp[106] = 1.0*s393;
            dwdp[107] = 1.0*s473;
            dwdp[108] = 1.0*s474;
            dwdp[109] = 1.0*s475;
            dwdp[110] = 1.0*s476;
            dwdp[111] = 1.0*s477;
            dwdp[112] = 1.0*s478;
            dwdp[113] = 1.0*s479;
            break;
        case 8:
            dwdp[4] = 1.0*s17*s51;
            dwdp[7] = 1.0*s18*s51;
            break;
        case 9:
            dwdp[17] = 1.0*s23*s51;
            dwdp[24] = 1.0*s31*s51;
            dwdp[29] = 1.0*s29*s51;
            dwdp[80] = 1.0*s30*s500;
            dwdp[81] = 1.0*s20*s51;
            dwdp[82] = 1.0*s24*s51;
            dwdp[83] = 1.0*s32*s51;
            break;
        case 10:
            dwdp[38] = 1.0*s517;
            dwdp[39] = 1.0*s523;
            dwdp[40] = 1.0*s520;
            dwdp[41] = 1.0*s521;
            dwdp[42] = 1.0*s522;
            dwdp[43] = 1.0*s518;
            dwdp[44] = 1.0*s519;
            dwdp[85] = 1.0*s530;
            dwdp[86] = 1.0*s531;
            dwdp[87] = 1.0*s527;
            dwdp[88] = 1.0*s529;
            dwdp[89] = 1.0*s528;
            dwdp[90] = 1.0*s526;
            dwdp[91] = 1.0*s525;
            dwdp[133] = 1.0*s533;
            dwdp[135] = 1.0*s535;
            break;
        case 11:
            dwdp[0] = 1.0*s3;
            break;
        case 12:
            dwdp[1] = 1.0*s3;
            break;
        case 13:
            dwdp[3] = 1.0*s17*s22;
            break;
        case 14:
            dwdp[30] = 1.0*s22;
            break;
        case 15:
            dwdp[36] = 1.0*s78;
            break;
        case 16:
            dwdp[37] = 1.0*s85;
            break;
        case 17:
            dwdp[67] = 1.0*s53;
            break;
        case 18:
            dwdp[68] = 1.0*s52;
            break;
        case 19:
            dwdp[131] = 1.0*s17*s33;
            break;
    }
}