#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Sengupta2015(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = s4;
    x_rdata[1] = s136;
    x_rdata[2] = s188;
    x_rdata[3] = s253;
    x_rdata[4] = s49;
    x_rdata[5] = s308;
    x_rdata[6] = s335;
    x_rdata[7] = s13;
    x_rdata[8] = s201;
    x_rdata[9] = s199;
    x_rdata[10] = s294;
    x_rdata[11] = s292;
    x_rdata[12] = s297;
    x_rdata[13] = s298;
    x_rdata[14] = s337;
    x_rdata[15] = s248;
    x_rdata[16] = s362;
    x_rdata[17] = s94;
    x_rdata[18] = s284;
    x_rdata[19] = s283;
    x_rdata[20] = s26;
    x_rdata[21] = s7;
    x_rdata[22] = s244;
    x_rdata[23] = s58;
    x_rdata[24] = s6;
    x_rdata[25] = s187;
    x_rdata[26] = s388;
    x_rdata[27] = s267;
    x_rdata[28] = s11;
    x_rdata[29] = s44;
    x_rdata[30] = s344;
    x_rdata[31] = s193;
    x_rdata[32] = s79;
    x_rdata[33] = s249;
    x_rdata[34] = s45;
    x_rdata[35] = s91;
    x_rdata[36] = s90;
    x_rdata[37] = s200;
    x_rdata[38] = s363;
    x_rdata[39] = s36;
    x_rdata[40] = s43;
    x_rdata[41] = s86;
    x_rdata[42] = s88;
    x_rdata[43] = s10;
    x_rdata[44] = s34;
    x_rdata[45] = s81;
    x_rdata[46] = s50;
    x_rdata[47] = s190;
    x_rdata[48] = s185;
    x_rdata[49] = s334;
    x_rdata[50] = s67;
    x_rdata[51] = s351;
    x_rdata[52] = s381;
    x_rdata[53] = s93;
    x_rdata[54] = s355;
    x_rdata[55] = s333;
    x_rdata[56] = s357;
    x_rdata[57] = s28;
    x_rdata[58] = s9;
    x_rdata[59] = s252;
    x_rdata[60] = s307;
    x_rdata[61] = s189;
    x_rdata[62] = s296;
    x_rdata[63] = s95;
    x_rdata[64] = s56;
    x_rdata[65] = s52;
    x_rdata[66] = s16;
    x_rdata[67] = s322;
    x_rdata[68] = s327;
    x_rdata[69] = s25;
    x_rdata[70] = s64;
    x_rdata[71] = s347;
    x_rdata[72] = s32;
    x_rdata[73] = s329;
    x_rdata[74] = s361;
    x_rdata[75] = s255;
    x_rdata[76] = s250;
    x_rdata[77] = s203;
    x_rdata[78] = s251;
    x_rdata[79] = s336;
    x_rdata[80] = s197;
    x_rdata[81] = s198;
    x_rdata[82] = s51;
    x_rdata[83] = s293;
    x_rdata[84] = s340;
    x_rdata[85] = s378;
    x_rdata[86] = s18;
    x_rdata[87] = s341;
    x_rdata[88] = s379;
    x_rdata[89] = s71;
    x_rdata[90] = s302;
    x_rdata[91] = s306;
    x_rdata[92] = s305;
    x_rdata[93] = s48;
    x_rdata[94] = s241;
    x_rdata[95] = s352;
    x_rdata[96] = s35;
    x_rdata[97] = s186;
    x_rdata[98] = s234;
    x_rdata[99] = s85;
    x_rdata[100] = s195;
    x_rdata[101] = s77;
    x_rdata[102] = s40;
    x_rdata[103] = s27;
    x_rdata[104] = s286;
    x_rdata[105] = s287;
    x_rdata[106] = s374;
    x_rdata[107] = s377;
    x_rdata[108] = s269;
    x_rdata[109] = s353;
    x_rdata[110] = s356;
    x_rdata[111] = s53;
    x_rdata[112] = s3;
    x_rdata[113] = s285;
    x_rdata[114] = s73;
    x_rdata[115] = s31;
    x_rdata[116] = s247;
    x_rdata[117] = s194;
    x_rdata[118] = s365;
    x_rdata[119] = s87;
    x_rdata[120] = s82;
    x_rdata[121] = s183;
    x_rdata[122] = s325;
    x_rdata[123] = s330;
    x_rdata[124] = s265;
    x_rdata[125] = s258;
    x_rdata[126] = s263;
    x_rdata[127] = s262;
    x_rdata[128] = s238;
    x_rdata[129] = s350;
    x_rdata[130] = s92;
    x_rdata[131] = s38;
    x_rdata[132] = s342;
    x_rdata[133] = s54;
    x_rdata[134] = s2;
    x_rdata[135] = s126;
    x_rdata[136] = s41;
    x_rdata[137] = s33;
    x_rdata[138] = s237;
    x_rdata[139] = s349;
    x_rdata[140] = s122;
    x_rdata[141] = s109;
    x_rdata[142] = s112;
    x_rdata[143] = s111;
    x_rdata[144] = s113;
    x_rdata[145] = s110;
    x_rdata[146] = s115;
    x_rdata[147] = s116;
    x_rdata[148] = s114;
    x_rdata[149] = s117;
    x_rdata[150] = s120;
    x_rdata[151] = s119;
    x_rdata[152] = s121;
    x_rdata[153] = s127;
    x_rdata[154] = s367;
    x_rdata[155] = s15;
    x_rdata[156] = s22;
    x_rdata[157] = s37;
    x_rdata[158] = s125;
    x_rdata[159] = s19;
    x_rdata[160] = s129;
    x_rdata[161] = s371;
    x_rdata[162] = s42;
    x_rdata[163] = s68;
    x_rdata[164] = s66;
    x_rdata[165] = s128;
    x_rdata[166] = s69;
    x_rdata[167] = s131;
    x_rdata[168] = s370;
    x_rdata[169] = s65;
    x_rdata[170] = s72;
    x_rdata[171] = s80;
    x_rdata[172] = s130;
    x_rdata[173] = s70;
    x_rdata[174] = s133;
    x_rdata[175] = s369;
    x_rdata[176] = s83;
    x_rdata[177] = s96;
    x_rdata[178] = s89;
    x_rdata[179] = s132;
    x_rdata[180] = s97;
    x_rdata[181] = s135;
    x_rdata[182] = s368;
    x_rdata[183] = s84;
    x_rdata[184] = s99;
    x_rdata[185] = s100;
    x_rdata[186] = s134;
    x_rdata[187] = s98;
    x_rdata[188] = s101;
    x_rdata[189] = s104;
    x_rdata[190] = s103;
    x_rdata[191] = s105;
    x_rdata[192] = s102;
    x_rdata[193] = s107;
    x_rdata[194] = s108;
    x_rdata[195] = s106;
    x_rdata[196] = s366;
    x_rdata[197] = s332;
    x_rdata[198] = s24;
    x_rdata[199] = s23;
    x_rdata[200] = s39;
    x_rdata[201] = s343;
    x_rdata[202] = s12;
    x_rdata[203] = s254;
    x_rdata[204] = s57;
    x_rdata[205] = s5;
    x_rdata[206] = s17;
    x_rdata[207] = s21;
    x_rdata[208] = s291;
    x_rdata[209] = s14;
    x_rdata[210] = s364;
    x_rdata[211] = s354;
    x_rdata[212] = s326;
    x_rdata[213] = s324;
    x_rdata[214] = s321;
    x_rdata[215] = s328;
    x_rdata[216] = s323;
    x_rdata[217] = s63;
    x_rdata[218] = s345;
    x_rdata[219] = s47;
    x_rdata[220] = s46;
    x_rdata[221] = s346;
    x_rdata[222] = s348;
    x_rdata[223] = s192;
    x_rdata[224] = s182;
    x_rdata[225] = s181;
    x_rdata[226] = s8;
    x_rdata[227] = s259;
    x_rdata[228] = s300;
    x_rdata[229] = s301;
    x_rdata[230] = s358;
    x_rdata[231] = s389;
    x_rdata[232] = s29;
    x_rdata[233] = s30;
    x_rdata[234] = s75;
    x_rdata[235] = s256;
    x_rdata[236] = s74;
    x_rdata[237] = s124;
    x_rdata[238] = s123;
    x_rdata[239] = s400;
}