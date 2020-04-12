#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_bachmann2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 2.5*dwdx1;
    JB[20] = 2.5*dwdx0;
    JB[25] = -2.5*dwdx0;
    JB[26] = -2.5*dwdx3;
    JB[27] = 2.5*dwdx2;
    JB[54] = 3.6363636363636362*dwdx4;
    JB[55] = -3.6363636363636362*dwdx4;
    JB[81] = 3.6363636363636362*dwdx5;
    JB[82] = -3.6363636363636362*dwdx5;
    JB[108] = 3.6363636363636362*dwdx6;
    JB[109] = -3.6363636363636362*dwdx6;
    JB[135] = 3.6363636363636362*dwdx7;
    JB[136] = -3.6363636363636362*dwdx7;
    JB[157] = -2.5*dwdx8;
    JB[162] = 3.6363636363636362*dwdx8;
    JB[190] = 2.5*dwdx9;
    JB[192] = -2.5*dwdx9;
    JB[216] = 2.5*dwdx10;
    JB[218] = -2.5*dwdx10;
    JB[243] = 2.5*dwdx11;
    JB[244] = 2.5*dwdx12;
    JB[256] = -2.5*dwdx13;
    JB[257] = 2.5*dwdx13;
    JB[258] = -2.5*dwdx12;
    JB[268] = -2.5*dwdx16;
    JB[270] = 2.5*dwdx16 + 2.5*dwdx17 + 2.5*dwdx18;
    JB[271] = 2.5*dwdx14;
    JB[272] = -2.5*dwdx14;
    JB[280] = 2.5*dwdx15;
    JB[283] = -2.5*dwdx17;
    JB[284] = -2.5*dwdx18;
    JB[285] = -2.5*dwdx15;
    JB[297] = 2.5*dwdx19;
    JB[298] = -2.5*dwdx19;
    JB[320] = -2.5*dwdx21 - 2.5*dwdx22 - 2.5*dwdx23 - 2.5*dwdx24;
    JB[322] = 2.5*dwdx21;
    JB[323] = -2.5*dwdx20;
    JB[324] = 2.5*dwdx20;
    JB[334] = 2.5*dwdx24;
    JB[335] = 2.5*dwdx22;
    JB[336] = 2.5*dwdx23;
    JB[346] = 2.5*dwdx25;
    JB[348] = -2.5*dwdx25 + 2.5*dwdx28 + 2.5*dwdx30;
    JB[351] = 2.5*dwdx29;
    JB[358] = 2.5*dwdx26 + 2.5*dwdx27;
    JB[360] = -2.5*dwdx31 - 2.5*dwdx32;
    JB[361] = -2.5*dwdx28 + 2.5*dwdx31;
    JB[362] = -2.5*dwdx30 + 2.5*dwdx32;
    JB[363] = -2.5*dwdx26 - 2.5*dwdx27;
    JB[377] = -2.5*dwdx34;
    JB[378] = 2.5*dwdx33;
    JB[405] = 3.6363636363636362*dwdx35;
    JB[406] = -3.6363636363636362*dwdx35;
    JB[432] = 3.6363636363636362*dwdx36;
    JB[433] = -3.6363636363636362*dwdx36;
    JB[459] = 3.6363636363636362*dwdx37;
    JB[460] = -3.6363636363636362*dwdx37;
    JB[486] = 3.6363636363636362*dwdx38;
    JB[487] = -3.6363636363636362*dwdx38;
    JB[508] = -2.5*dwdx39;
    JB[513] = 3.6363636363636362*dwdx39;
    JB[540] = 2.5*dwdx40 + 2.5*dwdx41;
    JB[545] = -2.5*dwdx40 - 2.5*dwdx41;
    JB[548] = -3.6363636363636362*dwdx43;
    JB[561] = -3.6363636363636362*dwdx44;
    JB[566] = -2.5*dwdx42;
    JB[567] = 3.6363636363636362*dwdx42;
    JB[580] = -2.5*dwdx49;
    JB[581] = 2.5*dwdx45;
    JB[583] = 2.5*dwdx46;
    JB[584] = -2.5*dwdx46;
    JB[592] = 2.5*dwdx47 + 2.5*dwdx48;
    JB[594] = 2.5*dwdx49;
    JB[597] = -2.5*dwdx47 - 2.5*dwdx48;
    JB[606] = -2.5*dwdx55;
    JB[607] = 2.5*dwdx50;
    JB[609] = 2.5*dwdx51;
    JB[610] = -2.5*dwdx51;
    JB[618] = 2.5*dwdx52 + 2.5*dwdx53;
    JB[620] = -2.5*dwdx54;
    JB[621] = 2.5*dwdx54 + 2.5*dwdx55;
    JB[623] = -2.5*dwdx52 - 2.5*dwdx53;
    JB[632] = -2.5*dwdx59;
    JB[635] = 2.5*dwdx56;
    JB[636] = -2.5*dwdx56;
    JB[644] = 2.5*dwdx57;
    JB[646] = -2.5*dwdx58;
    JB[648] = 2.5*dwdx58 + 2.5*dwdx59;
    JB[649] = -2.5*dwdx57;
    JB[671] = -3.6363636363636362*dwdx60;
    JB[675] = 2.5*dwdx60;
}