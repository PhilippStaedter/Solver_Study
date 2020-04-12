#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Ouzounoglou2014(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*re1_k1*s3;
    w[1] = 1.0*re2_k1*s3;
    w[2] = 1.0*k_2merForm*pow(s17, 2);
    w[3] = 1.0*re4_k1*s17*s22;
    w[4] = 1.0*k_WTasyn1_2merBindOnLamp*s17*s51;
    w[5] = 1.0*k_OligAutophagUptake*s18;
    w[6] = 1.0*k_OligomerForm*s17*s18;
    w[7] = 1.0*k_WTasyn1_2merBindOnLamp*s18*s51;
    w[8] = 1.0*k_OligAutophagUptake*s20;
    w[9] = 1.0*k_OligomerForm*s17*s20;
    w[10] = 1.0*k_ProteasomeBind*s20*s35;
    w[11] = 1.0*k_OligAutophagUptake*s24;
    w[12] = 1.0*k_OligomerForm*s17*s24;
    w[13] = 1.0*k_ProteasomeBind*s24*s35;
    w[14] = 1.0*k_OligAutophagUptake*s23;
    w[15] = 1.0*k_OligomerForm*s17*s23;
    w[16] = 1.0*k_ProteasomeBind*s23*s35;
    w[17] = 1.0*k_WTOligBindOnLamp*s23*s51;
    w[18] = 1.0*k_OligAutophagUptake*s32;
    w[19] = 1.0*k_OligomerForm*s17*s32;
    w[20] = 1.0*k_ProteasomeBind*s32*s35;
    w[21] = 1.0*k_OligAutophagUptake*s31;
    w[22] = 1.0*k_OligomerForm*s17*s31;
    w[23] = 1.0*k_ProteasomeBind*s31*s35;
    w[24] = 1.0*k_WTOligBindOnLamp*s31*s51;
    w[25] = 1.0*k_OligAutophagUptake*s30;
    w[26] = 1.0*k_OligomerForm*s17*s30;
    w[27] = 1.0*k_ProteasomeBind*s30*s35;
    w[28] = 1.0*k_ProteasomeBind*s29*s35;
    w[29] = 1.0*k_WTOligBindOnLamp*s29*s51;
    w[30] = 1.0*re31_k1*s22;
    w[31] = 1.0*k_2merForm*pow(s7, 2);
    w[32] = 1.0*k_2merForm*s536*s7;
    w[33] = 1.0*k_OligomerForm*s490*s7;
    w[34] = 1.0*k_OligomerForm*s489*s7;
    w[35] = 1.0*k_OligomerForm*s492*s7;
    w[36] = 1.0*re37_k1*s78;
    w[37] = 1.0*re38_k1*s85;
    w[38] = 1.0*k_M_autophagyDegr*s517;
    w[39] = 1.0*k_M_autophagyDegr*s523;
    w[40] = 1.0*k_M_autophagyDegr*s520;
    w[41] = 1.0*k_M_autophagyDegr*s521;
    w[42] = 1.0*k_M_autophagyDegr*s522;
    w[43] = 1.0*k_M_autophagyDegr*s518;
    w[44] = 1.0*k_M_autophagyDegr*s519;
    w[45] = 1.0*k_OligomerForm*s17*s29;
    w[46] = 1.0*k_OligAutophagUptake*s6;
    w[47] = 1.0*k_OligomerForm*s6*s7;
    w[48] = 1.0*k_OligAutophagUptake*s5;
    w[49] = 1.0*k_OligomerForm*s5*s7;
    w[50] = 1.0*k_ProteasomeBind*s35*s5;
    w[51] = 1.0*k_OligAutophagUptake*s2;
    w[52] = 1.0*k_OligomerForm*s2*s7;
    w[53] = 1.0*k_ProteasomeBind*s2*s35;
    w[54] = 1.0*k_OligAutophagUptake*s1;
    w[55] = 1.0*k_OligomerForm*s1*s7;
    w[56] = 1.0*k_ProteasomeBind*s1*s35;
    w[57] = 1.0*k_OligAutophagUptake*s21;
    w[58] = 1.0*k_OligomerForm*s21*s7;
    w[59] = 1.0*k_ProteasomeBind*s21*s35;
    w[60] = 1.0*k_OligAutophagUptake*s25;
    w[61] = 1.0*k_OligomerForm*s25*s7;
    w[62] = 1.0*k_ProteasomeBind*s25*s35;
    w[63] = 1.0*k_OligAutophagUptake*s26;
    w[64] = 1.0*k_OligomerForm*s26*s7;
    w[65] = 1.0*k_ProteasomeBind*s26*s35;
    w[66] = 1.0*k_ProteasomeBind*s27*s35;
    w[67] = 1.0*re69_k1*s53;
    w[68] = 1.0*re70_k1*s52;
    w[69] = 1.0*k_LampFreeWTasyn*s501;
    w[70] = 1.0*k_OligomerForm*s482*s7;
    w[71] = 1.0*k_OligomerForm*s483*s7;
    w[72] = 1.0*k_OligomerForm*s484*s7;
    w[73] = 1.0*k_OligomerForm*s491*s7;
    w[74] = 1.0*k_LampFreeWTasyn*s494;
    w[75] = 1.0*k_LampFreeWTasyn*s495;
    w[76] = 1.0*k_LampFreeWTasyn*s496;
    w[77] = 1.0*k_LampFreeWTasyn*s498;
    w[78] = 1.0*k_LampFreeWTasyn*s499;
    w[79] = 1.0*k_LampFreeWTasyn*s500;
    w[80] = 1.0*k_WTOligBindOnLamp*s30*s500;
    w[81] = 1.0*k_WTOligBindOnLamp*s20*s51;
    w[82] = 1.0*k_WTOligBindOnLamp*s24*s51;
    w[83] = 1.0*k_WTOligBindOnLamp*s32*s51;
    w[84] = 1.0*k_DopModWTasynLampBind*s51*s7;
    w[85] = 1.0*k_M_autophagyDegr*s530;
    w[86] = 1.0*k_M_autophagyDegr*s531;
    w[87] = 1.0*k_M_autophagyDegr*s527;
    w[88] = 1.0*k_M_autophagyDegr*s529;
    w[89] = 1.0*k_M_autophagyDegr*s528;
    w[90] = 1.0*k_M_autophagyDegr*s526;
    w[91] = 1.0*k_M_autophagyDegr*s525;
    w[92] = 1.0*k_2merForm*s17*s78;
    w[93] = 1.0*k_OligomerForm*s17*s85;
    w[94] = 1.0*k_OligomerForm*s17*s494;
    w[95] = 1.0*k_OligomerForm*s17*s495;
    w[96] = 1.0*k_OligomerForm*s17*s496;
    w[97] = 1.0*k_OligomerForm*s17*s498;
    w[98] = 1.0*k_OligomerForm*s17*s499;
    w[99] = 1.0*k_OligomerForm*s17*s500;
    w[100] = 1.0*k_ProtOligDegr*s381;
    w[101] = 1.0*k_ProtOligDegr*s383;
    w[102] = 1.0*k_ProtOligDegr*s385;
    w[103] = 1.0*k_ProtOligDegr*s387;
    w[104] = 1.0*k_ProtOligDegr*s389;
    w[105] = 1.0*k_ProtOligDegr*s391;
    w[106] = 1.0*k_ProtOligDegr*s393;
    w[107] = 1.0*k_ProtOligDegr*s473;
    w[108] = 1.0*k_ProtOligDegr*s474;
    w[109] = 1.0*k_ProtOligDegr*s475;
    w[110] = 1.0*k_ProtOligDegr*s476;
    w[111] = 1.0*k_ProtOligDegr*s477;
    w[112] = 1.0*k_ProtOligDegr*s478;
    w[113] = 1.0*k_ProtOligDegr*s479;
    w[114] = 1.0*k_ProteasomeBind*s33*s35;
    w[115] = 1.0*k_DisRate*s27;
    w[116] = 1.0*k_DisRate*s26;
    w[117] = 1.0*k_DisRate*s25;
    w[118] = 1.0*k_DisRate*s21;
    w[119] = 1.0*k_DisRate*s2;
    w[120] = 1.0*k_DisRate*s1;
    w[121] = 1.0*k_DisRate*s5;
    w[122] = 1.0*k_DisRate*s6;
    w[123] = 1.0*k_DisRate*s29;
    w[124] = 1.0*k_DisRate*s30;
    w[125] = 1.0*k_DisRate*s31;
    w[126] = 1.0*k_DisRate*s32;
    w[127] = 1.0*k_DisRate*s23;
    w[128] = 1.0*k_DisRate*s24;
    w[129] = 1.0*k_DisRate*s20;
    w[130] = 1.0*k_DisRate*s18;
    w[131] = 1.0*re133_k1*s17*s33;
    w[132] = 1.0*k_OligAutophagUptake*s17;
    w[133] = 1.0*k_M_autophagyDegr*s533;
    w[134] = 1.0*k_OligAutophagUptake*s7;
    w[135] = 1.0*k_M_autophagyDegr*s535;
}