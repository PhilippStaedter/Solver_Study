#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Qi2013a(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx10 - 1.0*dwdx2 - 1.0*dwdx3 - 2.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[2] = -1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13;
    JDiag[3] = -1.0*dwdx14;
    JDiag[4] = -1.0*dwdx15;
    JDiag[5] = -1.0*dwdx16;
    JDiag[6] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19;
    JDiag[7] = -1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx38 - 1.0*dwdx39 - 1.0*dwdx40 - 1.0*dwdx41 - 1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44 - 1.0*dwdx45;
    JDiag[8] = -1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54 - 1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63 - 1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66 - 1.0*dwdx67 - 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79 - 1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx82 - 1.0*dwdx83 - 1.0*dwdx84 - 1.0*dwdx85 - 1.0*dwdx86 - 1.0*dwdx87 - 1.0*dwdx88 - 1.0*dwdx89 - 1.0*dwdx90 - 1.0*dwdx91 - 1.0*dwdx92 - 1.0*dwdx93 - 1.0*dwdx94;
    JDiag[9] = -1.0*dwdx99;
    JDiag[10] = -1.0*dwdx100;
    JDiag[11] = -1.0*dwdx101 - 1.0*dwdx102 - 1.0*dwdx103 - 1.0*dwdx104 - 1.0*dwdx105 - 1.0*dwdx106 - 1.0*dwdx107 - 1.0*dwdx108 - 1.0*dwdx109 - 1.0*dwdx110 - 1.0*dwdx111 - 1.0*dwdx112 - 1.0*dwdx113 - 1.0*dwdx114 - 1.0*dwdx115 - 1.0*dwdx116 - 1.0*dwdx117 - 1.0*dwdx118 - 1.0*dwdx119 - 1.0*dwdx120 - 1.0*dwdx121 - 1.0*dwdx122 - 1.0*dwdx123 - 1.0*dwdx124;
    JDiag[12] = -1.0*dwdx125 - 1.0*dwdx126 - 1.0*dwdx127;
    JDiag[13] = -1.0*dwdx128 - 1.0*dwdx129 - 1.0*dwdx130;
    JDiag[14] = -1.0*dwdx131 - 1.0*dwdx132 - 1.0*dwdx133;
    JDiag[15] = -1.0*dwdx134 - 1.0*dwdx135 - 1.0*dwdx136 - 1.0*dwdx137;
    JDiag[16] = -1.0*dwdx138 - 1.0*dwdx139 - 1.0*dwdx140 - 1.0*dwdx141;
    JDiag[17] = -1.0*dwdx142 - 1.0*dwdx143 - 1.0*dwdx144 - 1.0*dwdx145;
    JDiag[18] = -1.0*dwdx146 - 1.0*dwdx147 - 1.0*dwdx148 - 1.0*dwdx149;
    JDiag[19] = -1.0*dwdx150 - 1.0*dwdx151 - 1.0*dwdx152;
    JDiag[20] = -1.0*dwdx153 - 1.0*dwdx154;
    JDiag[21] = -1.0*dwdx155 - 1.0*dwdx156;
    JDiag[22] = -1.0*dwdx157 - 1.0*dwdx158;
    JDiag[23] = -1.0*dwdx159 - 1.0*dwdx160;
    JDiag[24] = -1.0*dwdx161 - 1.0*dwdx162;
    JDiag[25] = -1.0*dwdx163 - 1.0*dwdx164;
    JDiag[26] = -1.0*dwdx165 - 1.0*dwdx166;
    JDiag[27] = -1.0*dwdx167 - 1.0*dwdx168;
    JDiag[28] = -1.0*dwdx169 - 1.0*dwdx170;
    JDiag[29] = -1.0*dwdx171 - 1.0*dwdx172;
    JDiag[30] = -1.0*dwdx173 - 1.0*dwdx174;
    JDiag[31] = -1.0*dwdx175 - 1.0*dwdx176;
    JDiag[32] = -1.0*dwdx177 - 1.0*dwdx178;
    JDiag[33] = -1.0*dwdx179 - 1.0*dwdx180 - 1.0*dwdx181;
    JDiag[34] = -1.0*dwdx183 - 1.0*dwdx184 - 1.0*dwdx185;
    JDiag[35] = -1.0*dwdx187 - 1.0*dwdx188 - 1.0*dwdx189;
    JDiag[36] = -1.0*dwdx191 - 1.0*dwdx192 - 1.0*dwdx193;
    JDiag[37] = -1.0*dwdx195 - 1.0*dwdx196 - 1.0*dwdx197;
    JDiag[38] = 1.0*dwdx199 + 1.0*dwdx200 + 1.0*dwdx201 + 1.0*dwdx202 + 1.0*dwdx203 + 1.0*dwdx204 + 1.0*dwdx205 + 1.0*dwdx206 + 1.0*dwdx207 + 1.0*dwdx208 + 1.0*dwdx209 + 1.0*dwdx210 + 1.0*dwdx211 + 1.0*dwdx212 + 1.0*dwdx213 + 1.0*dwdx214 + 1.0*dwdx215 + 1.0*dwdx216 + 1.0*dwdx217 + 1.0*dwdx218 + 1.0*dwdx219 + 1.0*dwdx220 + 1.0*dwdx221 + 1.0*dwdx222 + 1.0*dwdx223 + 1.0*dwdx224 + 1.0*dwdx225 + 1.0*dwdx226 + 1.0*dwdx227 + 1.0*dwdx228 + 1.0*dwdx229 + 1.0*dwdx230 + 1.0*dwdx231 + 1.0*dwdx232 + 1.0*dwdx233 + 1.0*dwdx234 + 1.0*dwdx235 + 1.0*dwdx236 + 1.0*dwdx237 + 1.0*dwdx238 + 1.0*dwdx239 + 1.0*dwdx240 + 1.0*dwdx241 + 1.0*dwdx242 + 1.0*dwdx243 + 1.0*dwdx244 + 1.0*dwdx245 + 1.0*dwdx246 + 1.0*dwdx247 + 1.0*dwdx248 + 1.0*dwdx249 + 1.0*dwdx250 + 1.0*dwdx251 + 1.0*dwdx252 + 1.0*dwdx253 + 1.0*dwdx254 + 1.0*dwdx255;
    JDiag[43] = -1.0*dwdx272 - 1.0*dwdx273 - 1.0*dwdx274 - 1.0*dwdx275 - 1.0*dwdx276 - 1.0*dwdx277 - 1.0*dwdx278 - 1.0*dwdx279 - 1.0*dwdx280 - 1.0*dwdx281 - 1.0*dwdx282 - 1.0*dwdx283;
    JDiag[44] = -1.0*dwdx284;
    JDiag[45] = -1.0*dwdx285;
    JDiag[47] = -1.0*dwdx288 - 1.0*dwdx289 - 1.0*dwdx290 - 2.0*dwdx291 - 1.0*dwdx292 - 1.0*dwdx293 - 1.0*dwdx294 - 1.0*dwdx295 - 1.0*dwdx296 - 1.0*dwdx297;
    JDiag[48] = -1.0*dwdx298 - 1.0*dwdx299;
    JDiag[49] = -1.0*dwdx300 - 1.0*dwdx301;
    JDiag[50] = -1.0*dwdx302;
    JDiag[51] = -1.0*dwdx303;
    JDiag[52] = -1.0*dwdx304 - 1.0*dwdx305 - 2.0*dwdx306 - 1.0*dwdx307 - 1.0*dwdx308 - 1.0*dwdx309 - 1.0*dwdx310 - 1.0*dwdx311 - 1.0*dwdx312;
    JDiag[53] = -1.0*dwdx313;
    JDiag[54] = -1.0*dwdx314 - 1.0*dwdx315 - 1.0*dwdx316;
    JDiag[55] = -1.0*dwdx317 - 1.0*dwdx318 - 1.0*dwdx319;
    JDiag[56] = -1.0*dwdx320 - 1.0*dwdx321 - 1.0*dwdx322;
    JDiag[57] = -1.0*dwdx323 - 1.0*dwdx324 - 1.0*dwdx325;
    JDiag[58] = -1.0*dwdx326 - 1.0*dwdx327 - 1.0*dwdx328 - 1.0*dwdx329;
    JDiag[59] = -1.0*dwdx330 - 1.0*dwdx331 - 1.0*dwdx332 - 1.0*dwdx333;
    JDiag[60] = -1.0*dwdx334 - 1.0*dwdx335 - 1.0*dwdx336 - 1.0*dwdx337;
    JDiag[61] = -1.0*dwdx338 - 1.0*dwdx339 - 1.0*dwdx340 - 1.0*dwdx341;
    JDiag[62] = -1.0*dwdx342 - 1.0*dwdx343 - 1.0*dwdx344;
    JDiag[63] = -1.0*dwdx345 - 1.0*dwdx346;
    JDiag[64] = -1.0*dwdx347 - 1.0*dwdx348;
    JDiag[65] = -1.0*dwdx349 - 1.0*dwdx350;
    JDiag[66] = -1.0*dwdx351 - 1.0*dwdx352;
    JDiag[67] = -1.0*dwdx353 - 1.0*dwdx354;
    JDiag[68] = -1.0*dwdx355 - 1.0*dwdx356;
    JDiag[69] = -1.0*dwdx357 - 1.0*dwdx358;
    JDiag[70] = -1.0*dwdx359 - 1.0*dwdx360;
    JDiag[71] = -1.0*dwdx361 - 1.0*dwdx362;
    JDiag[72] = -1.0*dwdx363 - 1.0*dwdx364;
    JDiag[73] = -1.0*dwdx365 - 1.0*dwdx366;
    JDiag[74] = -1.0*dwdx367 - 1.0*dwdx368;
    JDiag[75] = -1.0*dwdx369 - 1.0*dwdx370;
    JDiag[76] = -1.0*dwdx371 - 1.0*dwdx372 - 1.0*dwdx373 - 2.0*dwdx374 - 1.0*dwdx375 - 1.0*dwdx376 - 1.0*dwdx377 - 1.0*dwdx378 - 1.0*dwdx379 - 1.0*dwdx380;
    JDiag[77] = -1.0*dwdx381;
    JDiag[78] = -1.0*dwdx382;
    JDiag[79] = -1.0*dwdx383 - 2.0*dwdx384 - 1.0*dwdx385 - 1.0*dwdx386 - 1.0*dwdx387 - 1.0*dwdx388 - 1.0*dwdx389 - 1.0*dwdx390;
    JDiag[80] = -1.0*dwdx391;
    JDiag[81] = -1.0*dwdx392 - 1.0*dwdx393 - 1.0*dwdx394;
    JDiag[82] = -1.0*dwdx395 - 1.0*dwdx396 - 1.0*dwdx397;
    JDiag[83] = -1.0*dwdx398 - 1.0*dwdx399 - 1.0*dwdx400;
    JDiag[84] = -1.0*dwdx401 - 1.0*dwdx402 - 1.0*dwdx403;
    JDiag[85] = -1.0*dwdx404 - 1.0*dwdx405 - 1.0*dwdx406 - 1.0*dwdx407;
    JDiag[86] = -1.0*dwdx408 - 1.0*dwdx409 - 1.0*dwdx410 - 1.0*dwdx411;
    JDiag[87] = -1.0*dwdx412 - 1.0*dwdx413 - 1.0*dwdx414 - 1.0*dwdx415;
    JDiag[88] = -1.0*dwdx416 - 1.0*dwdx417 - 1.0*dwdx418 - 1.0*dwdx419;
    JDiag[89] = -1.0*dwdx420 - 1.0*dwdx421 - 1.0*dwdx422;
    JDiag[90] = -1.0*dwdx423 - 1.0*dwdx424;
    JDiag[91] = -1.0*dwdx425 - 1.0*dwdx426;
    JDiag[92] = -1.0*dwdx427 - 1.0*dwdx428;
    JDiag[93] = -1.0*dwdx429 - 1.0*dwdx430;
    JDiag[94] = -1.0*dwdx431 - 1.0*dwdx432;
    JDiag[95] = -1.0*dwdx433 - 1.0*dwdx434;
    JDiag[96] = -1.0*dwdx435 - 1.0*dwdx436;
    JDiag[97] = -1.0*dwdx437 - 1.0*dwdx438;
    JDiag[98] = -1.0*dwdx439 - 1.0*dwdx440;
    JDiag[99] = -1.0*dwdx441 - 1.0*dwdx442;
    JDiag[100] = -1.0*dwdx443 - 1.0*dwdx444;
    JDiag[101] = -1.0*dwdx445 - 1.0*dwdx446;
    JDiag[102] = -1.0*dwdx447 - 1.0*dwdx448;
    JDiag[103] = -1.0*dwdx449 - 1.0*dwdx450 - 1.0*dwdx451;
    JDiag[104] = -1.0*dwdx453 - 1.0*dwdx454 - 1.0*dwdx455;
    JDiag[105] = -1.0*dwdx457 - 1.0*dwdx458 - 1.0*dwdx459;
    JDiag[106] = -1.0*dwdx461 - 1.0*dwdx462 - 1.0*dwdx463;
    JDiag[107] = -1.0*dwdx465 - 1.0*dwdx466 - 1.0*dwdx467;
    JDiag[108] = -1.0*dwdx469 - 1.0*dwdx470 - 1.0*dwdx471;
    JDiag[109] = -1.0*dwdx473 - 1.0*dwdx474 - 1.0*dwdx475;
    JDiag[110] = -1.0*dwdx477 - 1.0*dwdx478 - 1.0*dwdx479;
    JDiag[111] = -1.0*dwdx481 - 1.0*dwdx482 - 1.0*dwdx483;
    JDiag[112] = -1.0*dwdx485 - 1.0*dwdx486 - 1.0*dwdx487;
    JDiag[113] = -1.0*dwdx489 - 1.0*dwdx490 - 1.0*dwdx491;
    JDiag[114] = -1.0*dwdx493 - 1.0*dwdx494 - 1.0*dwdx495;
    JDiag[115] = -1.0*dwdx497 - 1.0*dwdx498 - 1.0*dwdx499;
    JDiag[116] = -1.0*dwdx501 - 1.0*dwdx502 - 1.0*dwdx503;
    JDiag[117] = -1.0*dwdx505 - 1.0*dwdx506 - 1.0*dwdx507;
    JDiag[118] = -1.0*dwdx509 - 1.0*dwdx510 - 1.0*dwdx511;
    JDiag[119] = -1.0*dwdx513 - 1.0*dwdx514 - 1.0*dwdx515;
    JDiag[120] = -1.0*dwdx517 - 1.0*dwdx518 - 1.0*dwdx519;
    JDiag[121] = -1.0*dwdx521 - 1.0*dwdx522 - 1.0*dwdx523;
    JDiag[122] = -1.0*dwdx525 - 1.0*dwdx526 - 1.0*dwdx527;
}