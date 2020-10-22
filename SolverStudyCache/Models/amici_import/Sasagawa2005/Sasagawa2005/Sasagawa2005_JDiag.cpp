#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Sasagawa2005(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 1.0*dwdx0 - 1.0*dwdx1;
    JDiag[1] = 1.0*dwdx2 - 2.0*dwdx3;
    JDiag[2] = 1.0*dwdx4 - 1.0*dwdx5;
    JDiag[3] = -1.0*dwdx6 - 1.0*dwdx7;
    JDiag[4] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx14 - 1.0*dwdx15 + 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[5] = -1.0*dwdx16 - 1.0*dwdx17;
    JDiag[6] = 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JDiag[7] = -1.0*dwdx25 - 1.0*dwdx26;
    JDiag[8] = 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29;
    JDiag[9] = -1.0*dwdx30 - 1.0*dwdx31;
    JDiag[10] = -1.0*dwdx32;
    JDiag[11] = -1.0*dwdx33 - 1.0*dwdx34 - 1.0*dwdx35 - 1.0*dwdx36;
    JDiag[12] = -1.0*dwdx37 - 1.0*dwdx38 - 1.0*dwdx39 - 1.0*dwdx40;
    JDiag[13] = 1.0*dwdx41 - 1.0*dwdx42;
    JDiag[14] = -1.0*dwdx43;
    JDiag[15] = -1.0*dwdx44 - 1.0*dwdx45 - 1.0*dwdx46 - 1.0*dwdx47;
    JDiag[16] = -1.0*dwdx48 - 1.0*dwdx49;
    JDiag[17] = 1.0*dwdx50 - 1.0*dwdx51 - 1.0*dwdx52;
    JDiag[18] = 1.0*dwdx54 - 1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61;
    JDiag[19] = 1.0*dwdx62 - 1.0*dwdx63 - 1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66 - 1.0*dwdx67;
    JDiag[20] = 1.0*dwdx68 - 1.0*dwdx69 + 1.0*dwdx71 - 1.0*dwdx72;
    JDiag[21] = 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75 + 1.0*dwdx76;
    JDiag[22] = 1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79;
    JDiag[23] = -1.0*dwdx81 - 1.0*dwdx82 - 1.0*dwdx83 - 1.0*dwdx84 - 1.0*dwdx85 - 1.0*dwdx86 - 1.0*dwdx87;
    JDiag[24] = -1.0*dwdx88;
    JDiag[25] = -1.0*dwdx89;
    JDiag[26] = -1.0*dwdx90 - 1.0*dwdx91;
    JDiag[27] = -1.0*dwdx92 - 1.0*dwdx93 - 1.0*dwdx94;
    JDiag[29] = -1.0*dwdx99;
    JDiag[31] = -1.0*dwdx102;
    JDiag[32] = -1.0*dwdx103 + 1.0*dwdx104;
    JDiag[33] = -1.0*dwdx105 - 1.0*dwdx106 - 1.0*dwdx107 - 1.0*dwdx108 - 1.0*dwdx109 - 1.0*dwdx110;
    JDiag[34] = 1.0*dwdx112 - 1.0*dwdx113 - 1.0*dwdx114;
    JDiag[35] = -1.0*dwdx115 - 1.0*dwdx116 - 1.0*dwdx117 - 1.0*dwdx118 - 1.0*dwdx119 - 1.0*dwdx120;
    JDiag[36] = 1.0*dwdx121 - 1.0*dwdx122 - 1.0*dwdx123 - 1.0*dwdx124;
    JDiag[37] = -1.0*dwdx125 + 1.0*dwdx126 - 1.0*dwdx127 - 1.0*dwdx128 - 1.0*dwdx129;
    JDiag[38] = -1.0*dwdx130;
    JDiag[39] = 1.0*dwdx132 - 1.0*dwdx133 + 1.0*dwdx135;
    JDiag[40] = -1.0*dwdx136;
    JDiag[41] = 1.0*dwdx138 - 1.0*dwdx140;
    JDiag[42] = -1.0*dwdx141;
    JDiag[43] = 1.0*dwdx142 + 1.0*dwdx143 - 1.0*dwdx144 + 1.0*dwdx147;
    JDiag[44] = -1.0*dwdx148;
    JDiag[45] = -1.0*dwdx149;
    JDiag[47] = 1.0*dwdx150 - 1.0*dwdx151 - 1.0*dwdx152 - 1.0*dwdx153 - 1.0*dwdx154 - 1.0*dwdx155 - 1.0*dwdx156;
    JDiag[48] = 1.0*dwdx157 - 1.0*dwdx158 - 1.0*dwdx159 + 1.0*dwdx161;
    JDiag[49] = 1.0*dwdx162 + 1.0*dwdx163 - 1.0*dwdx164;
    JDiag[50] = -1.0*dwdx167 - 1.0*dwdx168 - 1.0*dwdx169 - 1.0*dwdx170 - 1.0*dwdx171;
    JDiag[51] = 1.0*dwdx173 - 1.0*dwdx174 - 1.0*dwdx175;
    JDiag[52] = 1.0*dwdx176;
    JDiag[53] = -1.0*dwdx180 - 1.0*dwdx181 - 1.0*dwdx182 - 1.0*dwdx183 - 1.0*dwdx184;
    JDiag[54] = 1.0*dwdx186 - 1.0*dwdx187 - 1.0*dwdx188 + 1.0*dwdx189;
    JDiag[55] = -1.0*dwdx190;
    JDiag[56] = -1.0*dwdx191 - 1.0*dwdx192 - 1.0*dwdx193 - 1.0*dwdx194;
    JDiag[57] = -1.0*dwdx195;
    JDiag[58] = 1.0*dwdx196 - 1.0*dwdx197 - 1.0*dwdx198 - 1.0*dwdx199 - 1.0*dwdx200 - 1.0*dwdx201;
    JDiag[59] = 1.0*dwdx202 - 1.0*dwdx203 - 1.0*dwdx204 - 1.0*dwdx205 - 1.0*dwdx206 - 1.0*dwdx207;
    JDiag[60] = -1.0*dwdx208 - 1.0*dwdx209;
    JDiag[61] = -2.0*dwdx210 - 1.0*dwdx211;
    JDiag[62] = -1.0*dwdx212 - 1.0*dwdx213 - 1.0*dwdx214 - 1.0*dwdx215 - 1.0*dwdx216 - 1.0*dwdx217 - 1.0*dwdx218;
    JDiag[63] = 1.0*dwdx220 - 1.0*dwdx221 - 1.0*dwdx222 - 1.0*dwdx223 - 1.0*dwdx224;
    JDiag[64] = -1.0*dwdx225 - 1.0*dwdx226 - 1.0*dwdx227;
    JDiag[65] = 1.0*dwdx228 - 1.0*dwdx229;
    JDiag[66] = -1.0*dwdx230 + 1.0*dwdx231 - 1.0*dwdx232;
    JDiag[67] = 1.0*dwdx235 - 1.0*dwdx237;
    JDiag[68] = 1.0*dwdx238 - 1.0*dwdx239 - 1.0*dwdx240 - 1.0*dwdx241;
    JDiag[69] = 1.0*dwdx243 - 1.0*dwdx244 - 1.0*dwdx245;
    JDiag[70] = 1.0*dwdx246 - 1.0*dwdx247 - 1.0*dwdx248 - 1.0*dwdx249;
    JDiag[71] = 1.0*dwdx251 - 1.0*dwdx252 - 1.0*dwdx253 - 1.0*dwdx254;
    JDiag[72] = 1.0*dwdx256 - 1.0*dwdx257 - 1.0*dwdx258 - 1.0*dwdx259;
    JDiag[73] = 1.0*dwdx261 - 1.0*dwdx262 - 1.0*dwdx263;
    JDiag[74] = 1.0*dwdx264 - 1.0*dwdx265 - 1.0*dwdx266;
    JDiag[75] = 1.0*dwdx267 - 1.0*dwdx268 - 1.0*dwdx269;
    JDiag[76] = 1.0*dwdx271 - 1.0*dwdx272;
    JDiag[77] = 1.0*dwdx274 - 1.0*dwdx275 + 1.0*dwdx276 - 1.0*dwdx277;
    JDiag[78] = -1.0*dwdx279 + 1.0*dwdx280 - 1.0*dwdx281;
    JDiag[79] = 1.0*dwdx283 + 1.0*dwdx284 - 1.0*dwdx285;
    JDiag[80] = 1.0*dwdx286 - 1.0*dwdx287;
    JDiag[81] = 1.0*dwdx288 - 1.0*dwdx289;
    JDiag[82] = 1.0*dwdx290 - 1.0*dwdx291;
    JDiag[83] = 1.0*dwdx292 - 1.0*dwdx293;
    JDiag[84] = 1.0*dwdx294 - 1.0*dwdx295;
    JDiag[85] = 1.0*dwdx296 - 1.0*dwdx297;
    JDiag[86] = 1.0*dwdx298 - 1.0*dwdx299;
    JDiag[87] = 1.0*dwdx300 - 1.0*dwdx301;
    JDiag[88] = 1.0*dwdx302 - 1.0*dwdx303;
    JDiag[89] = 1.0*dwdx304 - 1.0*dwdx305;
    JDiag[90] = 1.0*dwdx306 - 1.0*dwdx307;
    JDiag[91] = 1.0*dwdx308 - 1.0*dwdx309;
    JDiag[92] = 1.0*dwdx310 - 1.0*dwdx311;
    JDiag[93] = 1.0*dwdx312 - 1.0*dwdx313;
}