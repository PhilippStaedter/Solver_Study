#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_bachmann2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -2.5*dwdx1;
    J[1] = 2.5*dwdx3;
    J[27] = -2.5*dwdx2;
    J[32] = 2.5*dwdx8;
    J[54] = -3.6363636363636362*dwdx4;
    J[73] = 3.6363636363636362*dwdx43;
    J[80] = 3.6363636363636362*dwdx4;
    J[81] = -3.6363636363636362*dwdx5;
    J[107] = 3.6363636363636362*dwdx5;
    J[108] = -3.6363636363636362*dwdx6;
    J[134] = 3.6363636363636362*dwdx6;
    J[135] = -3.6363636363636362*dwdx7;
    J[161] = 3.6363636363636362*dwdx7;
    J[162] = -3.6363636363636362*dwdx8;
    J[215] = -2.5*dwdx9;
    J[216] = -2.5*dwdx10;
    J[218] = 2.5*dwdx16;
    J[220] = 2.5*dwdx21 + 2.5*dwdx22 + 2.5*dwdx23 + 2.5*dwdx24;
    J[221] = -2.5*dwdx25;
    J[230] = 2.5*dwdx49;
    J[231] = 2.5*dwdx55;
    J[232] = 2.5*dwdx59;
    J[243] = -2.5*dwdx11;
    J[256] = -2.5*dwdx45;
    J[257] = -2.5*dwdx50;
    J[267] = 2.5*dwdx9;
    J[268] = 2.5*dwdx10;
    J[269] = -2.5*dwdx12;
    J[270] = -2.5*dwdx16 - 2.5*dwdx17 - 2.5*dwdx18;
    J[272] = -2.5*dwdx21;
    J[273] = 2.5*dwdx25 - 2.5*dwdx28 - 2.5*dwdx30;
    J[296] = -2.5*dwdx14;
    J[297] = -2.5*dwdx19;
    J[298] = 2.5*dwdx20;
    J[308] = -2.5*dwdx46;
    J[309] = -2.5*dwdx51;
    J[310] = -2.5*dwdx56;
    J[322] = 2.5*dwdx14;
    J[323] = 2.5*dwdx19;
    J[324] = -2.5*dwdx20;
    J[334] = 2.5*dwdx46;
    J[335] = 2.5*dwdx51;
    J[336] = 2.5*dwdx56;
    J[351] = -2.5*dwdx29;
    J[352] = 2.5*dwdx34;
    J[378] = -2.5*dwdx33;
    J[383] = 2.5*dwdx39;
    J[405] = -3.6363636363636362*dwdx35;
    J[411] = 3.6363636363636362*dwdx44;
    J[431] = 3.6363636363636362*dwdx35;
    J[432] = -3.6363636363636362*dwdx36;
    J[458] = 3.6363636363636362*dwdx36;
    J[459] = -3.6363636363636362*dwdx37;
    J[485] = 3.6363636363636362*dwdx37;
    J[486] = -3.6363636363636362*dwdx38;
    J[512] = 3.6363636363636362*dwdx38;
    J[513] = -3.6363636363636362*dwdx39;
    J[520] = -2.5*dwdx0;
    J[530] = -2.5*dwdx15;
    J[533] = -2.5*dwdx26 - 2.5*dwdx27;
    J[540] = -2.5*dwdx40 - 2.5*dwdx41;
    J[541] = 2.5*dwdx42;
    J[542] = -2.5*dwdx47 - 2.5*dwdx48;
    J[543] = -2.5*dwdx52 - 2.5*dwdx53;
    J[544] = -2.5*dwdx57;
    J[567] = -3.6363636363636362*dwdx42;
    J[571] = 3.6363636363636362*dwdx60;
    J[581] = 2.5*dwdx13;
    J[584] = -2.5*dwdx24;
    J[585] = 2.5*dwdx31 + 2.5*dwdx32;
    J[594] = -2.5*dwdx49;
    J[595] = 2.5*dwdx54;
    J[596] = 2.5*dwdx58;
    J[607] = -2.5*dwdx13;
    J[608] = 2.5*dwdx17;
    J[610] = -2.5*dwdx22;
    J[611] = 2.5*dwdx28 - 2.5*dwdx31;
    J[621] = -2.5*dwdx54 - 2.5*dwdx55;
    J[633] = 2.5*dwdx12;
    J[634] = 2.5*dwdx18;
    J[636] = -2.5*dwdx23;
    J[637] = 2.5*dwdx30 - 2.5*dwdx32;
    J[648] = -2.5*dwdx58 - 2.5*dwdx59;
    J[650] = 2.5*dwdx0;
    J[660] = 2.5*dwdx15;
    J[663] = 2.5*dwdx26 + 2.5*dwdx27;
    J[670] = 2.5*dwdx40 + 2.5*dwdx41;
    J[672] = 2.5*dwdx47 + 2.5*dwdx48;
    J[673] = 2.5*dwdx52 + 2.5*dwdx53;
    J[674] = 2.5*dwdx57;
    J[675] = -2.5*dwdx60;
}