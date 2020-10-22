#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Levchenko2000a(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1;
    JDiag[1] = -1.0*dwdx2 - 1.0*dwdx3;
    JDiag[2] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15 - 1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[3] = -1.0*dwdx21;
    JDiag[4] = -1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33 - 1.0*dwdx34;
    JDiag[5] = -1.0*dwdx35 - 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx38 - 1.0*dwdx39 - 1.0*dwdx40 - 1.0*dwdx41 - 1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44 - 1.0*dwdx45 - 1.0*dwdx46 - 1.0*dwdx47 - 1.0*dwdx48;
    JDiag[6] = -1.0*dwdx49 - 1.0*dwdx50 - 1.0*dwdx51 - 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54 - 1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61;
    JDiag[7] = -1.0*dwdx62 - 1.0*dwdx63 - 1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66 - 1.0*dwdx67 - 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73 - 1.0*dwdx74;
    JDiag[8] = -1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79 - 1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx82 - 1.0*dwdx83 - 1.0*dwdx84 - 1.0*dwdx85 - 1.0*dwdx86 - 1.0*dwdx87 - 1.0*dwdx88;
    JDiag[9] = -1.0*dwdx100 - 1.0*dwdx101 - 1.0*dwdx102 - 1.0*dwdx103 - 1.0*dwdx89 - 1.0*dwdx90 - 1.0*dwdx91 - 1.0*dwdx92 - 1.0*dwdx93 - 1.0*dwdx94 - 1.0*dwdx95 - 1.0*dwdx96 - 1.0*dwdx97 - 1.0*dwdx98 - 1.0*dwdx99;
    JDiag[10] = -1.0*dwdx104 - 1.0*dwdx105 - 1.0*dwdx106 - 1.0*dwdx107 - 1.0*dwdx108 - 1.0*dwdx109 - 1.0*dwdx110 - 1.0*dwdx111 - 1.0*dwdx112 - 1.0*dwdx113 - 1.0*dwdx114 - 1.0*dwdx115 - 1.0*dwdx116 - 1.0*dwdx117 - 1.0*dwdx118 - 1.0*dwdx119 - 1.0*dwdx120;
    JDiag[11] = -1.0*dwdx121 - 1.0*dwdx122 - 1.0*dwdx123 - 1.0*dwdx124 - 1.0*dwdx125 - 1.0*dwdx126 - 1.0*dwdx127 - 1.0*dwdx128 - 1.0*dwdx129 - 1.0*dwdx130 - 1.0*dwdx131 - 1.0*dwdx132 - 1.0*dwdx133 - 1.0*dwdx134 - 1.0*dwdx135 - 1.0*dwdx136 - 1.0*dwdx137 - 1.0*dwdx138 - 1.0*dwdx139;
    JDiag[12] = -1.0*dwdx140 - 1.0*dwdx141;
    JDiag[13] = -1.0*dwdx142 - 1.0*dwdx143;
    JDiag[14] = -1.0*dwdx144 - 1.0*dwdx145;
    JDiag[15] = -1.0*dwdx146 - 1.0*dwdx147;
    JDiag[16] = -1.0*dwdx148 - 1.0*dwdx149;
    JDiag[17] = -1.0*dwdx150 - 1.0*dwdx151;
    JDiag[18] = -1.0*dwdx152 - 1.0*dwdx153;
    JDiag[19] = -1.0*dwdx154 - 1.0*dwdx155;
    JDiag[20] = -1.0*dwdx156 - 1.0*dwdx157;
    JDiag[21] = -1.0*dwdx158 - 1.0*dwdx159;
    JDiag[22] = -1.0*dwdx160 - 1.0*dwdx161 - 1.0*dwdx162 - 1.0*dwdx163 - 1.0*dwdx164 - 1.0*dwdx165 - 1.0*dwdx166 - 1.0*dwdx167;
    JDiag[23] = -1.0*dwdx168 - 1.0*dwdx169 - 1.0*dwdx170 - 1.0*dwdx171 - 1.0*dwdx172 - 1.0*dwdx173 - 1.0*dwdx174 - 1.0*dwdx175;
    JDiag[24] = -1.0*dwdx176 - 1.0*dwdx177 - 1.0*dwdx178 - 1.0*dwdx179 - 1.0*dwdx180 - 1.0*dwdx181 - 1.0*dwdx182;
    JDiag[25] = -1.0*dwdx183 - 1.0*dwdx184 - 1.0*dwdx185 - 1.0*dwdx186 - 1.0*dwdx187 - 1.0*dwdx188;
    JDiag[26] = -1.0*dwdx189 - 1.0*dwdx190 - 1.0*dwdx191 - 1.0*dwdx192 - 1.0*dwdx193 - 1.0*dwdx194;
    JDiag[27] = -1.0*dwdx195 - 1.0*dwdx196 - 1.0*dwdx197 - 1.0*dwdx198 - 1.0*dwdx199 - 1.0*dwdx200;
    JDiag[28] = -1.0*dwdx201 - 1.0*dwdx202 - 1.0*dwdx203 - 1.0*dwdx204 - 1.0*dwdx205 - 1.0*dwdx206;
    JDiag[29] = -1.0*dwdx207 - 1.0*dwdx208 - 1.0*dwdx209 - 1.0*dwdx210 - 1.0*dwdx211 - 1.0*dwdx212;
    JDiag[30] = -1.0*dwdx213 - 1.0*dwdx214 - 1.0*dwdx215 - 1.0*dwdx216 - 1.0*dwdx217 - 1.0*dwdx218;
    JDiag[31] = -1.0*dwdx219 - 1.0*dwdx220 - 1.0*dwdx221 - 1.0*dwdx222 - 1.0*dwdx223 - 1.0*dwdx224;
    JDiag[32] = -1.0*dwdx225 - 1.0*dwdx226 - 1.0*dwdx227 - 1.0*dwdx228 - 1.0*dwdx229 - 1.0*dwdx230;
    JDiag[33] = -1.0*dwdx231 - 1.0*dwdx232 - 1.0*dwdx233 - 1.0*dwdx234 - 1.0*dwdx235;
    JDiag[34] = -1.0*dwdx236 - 1.0*dwdx237 - 1.0*dwdx238 - 1.0*dwdx239 - 1.0*dwdx240 - 1.0*dwdx241;
    JDiag[35] = -1.0*dwdx242 - 1.0*dwdx243 - 1.0*dwdx244 - 1.0*dwdx245 - 1.0*dwdx246 - 1.0*dwdx247;
    JDiag[36] = -1.0*dwdx248 - 1.0*dwdx249 - 1.0*dwdx250 - 1.0*dwdx251 - 1.0*dwdx252;
    JDiag[37] = -1.0*dwdx253 - 1.0*dwdx254 - 1.0*dwdx255 - 1.0*dwdx256;
    JDiag[38] = -1.0*dwdx257 - 1.0*dwdx258 - 1.0*dwdx259 - 1.0*dwdx260;
    JDiag[39] = -1.0*dwdx261 - 1.0*dwdx262 - 1.0*dwdx263 - 1.0*dwdx264;
    JDiag[40] = -1.0*dwdx265 - 1.0*dwdx266 - 1.0*dwdx267 - 1.0*dwdx268;
    JDiag[41] = -1.0*dwdx269 - 1.0*dwdx270 - 1.0*dwdx271 - 1.0*dwdx272;
    JDiag[42] = -1.0*dwdx273 - 1.0*dwdx274 - 1.0*dwdx275 - 1.0*dwdx276;
    JDiag[43] = -1.0*dwdx277 - 1.0*dwdx278 - 1.0*dwdx279 - 1.0*dwdx280 - 1.0*dwdx281;
    JDiag[44] = -1.0*dwdx282 - 1.0*dwdx283 - 1.0*dwdx284 - 1.0*dwdx285 - 1.0*dwdx286;
    JDiag[45] = -1.0*dwdx287 - 1.0*dwdx288 - 1.0*dwdx289 - 1.0*dwdx290;
    JDiag[46] = -1.0*dwdx291 - 1.0*dwdx292 - 1.0*dwdx293 - 1.0*dwdx294 - 1.0*dwdx295 - 1.0*dwdx296;
    JDiag[47] = -1.0*dwdx297 - 1.0*dwdx298 - 1.0*dwdx299 - 1.0*dwdx300 - 1.0*dwdx301 - 1.0*dwdx302;
    JDiag[48] = -1.0*dwdx303 - 1.0*dwdx304 - 1.0*dwdx305 - 1.0*dwdx306 - 1.0*dwdx307;
    JDiag[49] = -1.0*dwdx308 - 1.0*dwdx309 - 1.0*dwdx310 - 1.0*dwdx311;
    JDiag[50] = -1.0*dwdx312 - 1.0*dwdx313 - 1.0*dwdx314 - 1.0*dwdx315;
    JDiag[51] = -1.0*dwdx316 - 1.0*dwdx317 - 1.0*dwdx318 - 1.0*dwdx319;
    JDiag[52] = -1.0*dwdx320 - 1.0*dwdx321 - 1.0*dwdx322 - 1.0*dwdx323;
    JDiag[53] = -1.0*dwdx324 - 1.0*dwdx325 - 1.0*dwdx326 - 1.0*dwdx327;
    JDiag[54] = -1.0*dwdx328 - 1.0*dwdx329 - 1.0*dwdx330 - 1.0*dwdx331;
    JDiag[55] = -1.0*dwdx332 - 1.0*dwdx333 - 1.0*dwdx334 - 1.0*dwdx335 - 1.0*dwdx336;
    JDiag[56] = -1.0*dwdx337 - 1.0*dwdx338 - 1.0*dwdx339 - 1.0*dwdx340 - 1.0*dwdx341;
    JDiag[57] = -1.0*dwdx342 - 1.0*dwdx343 - 1.0*dwdx344 - 1.0*dwdx345;
    JDiag[58] = -1.0*dwdx346 - 1.0*dwdx347 - 1.0*dwdx348 - 1.0*dwdx349 - 1.0*dwdx350 - 1.0*dwdx351;
    JDiag[59] = -1.0*dwdx352 - 1.0*dwdx353 - 1.0*dwdx354 - 1.0*dwdx355 - 1.0*dwdx356 - 1.0*dwdx357;
    JDiag[60] = -1.0*dwdx358 - 1.0*dwdx359 - 1.0*dwdx360 - 1.0*dwdx361 - 1.0*dwdx362;
    JDiag[61] = -1.0*dwdx363 - 1.0*dwdx364 - 1.0*dwdx365 - 1.0*dwdx366;
    JDiag[62] = -1.0*dwdx367 - 1.0*dwdx368 - 1.0*dwdx369 - 1.0*dwdx370;
    JDiag[63] = -1.0*dwdx371 - 1.0*dwdx372 - 1.0*dwdx373 - 1.0*dwdx374;
    JDiag[64] = -1.0*dwdx375 - 1.0*dwdx376 - 1.0*dwdx377 - 1.0*dwdx378;
    JDiag[65] = -1.0*dwdx379 - 1.0*dwdx380 - 1.0*dwdx381 - 1.0*dwdx382;
    JDiag[66] = -1.0*dwdx383 - 1.0*dwdx384 - 1.0*dwdx385 - 1.0*dwdx386;
    JDiag[67] = -1.0*dwdx387 - 1.0*dwdx388 - 1.0*dwdx389 - 1.0*dwdx390;
    JDiag[68] = -1.0*dwdx391 - 1.0*dwdx392 - 1.0*dwdx393 - 1.0*dwdx394;
    JDiag[69] = -1.0*dwdx395 - 1.0*dwdx396 - 1.0*dwdx397;
    JDiag[70] = -1.0*dwdx398 - 1.0*dwdx399;
    JDiag[71] = -1.0*dwdx400 - 1.0*dwdx401;
    JDiag[72] = -1.0*dwdx402 - 1.0*dwdx403;
    JDiag[73] = -1.0*dwdx404 - 1.0*dwdx405;
    JDiag[74] = -1.0*dwdx406 - 1.0*dwdx407;
    JDiag[75] = -1.0*dwdx408 - 1.0*dwdx409;
    JDiag[76] = -1.0*dwdx410 - 1.0*dwdx411;
    JDiag[77] = -1.0*dwdx412 - 1.0*dwdx413;
    JDiag[78] = -1.0*dwdx414 - 1.0*dwdx415;
    JDiag[79] = -1.0*dwdx416 - 1.0*dwdx417;
    JDiag[80] = -1.0*dwdx418 - 1.0*dwdx419;
    JDiag[81] = -1.0*dwdx420 - 1.0*dwdx421;
    JDiag[82] = -1.0*dwdx422 - 1.0*dwdx423;
    JDiag[83] = -1.0*dwdx424 - 1.0*dwdx425;
    JDiag[84] = -1.0*dwdx426 - 1.0*dwdx427;
    JDiag[85] = -1.0*dwdx428 - 1.0*dwdx429;
}