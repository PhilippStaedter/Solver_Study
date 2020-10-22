#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model6_levchenko1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[1] = -1.0*dwdx0;
    JB[14] = 1.0*dwdx0;
    JB[22] = -1.0*dwdx1;
    JB[23] = 1.0*dwdx1 + 1.0*dwdx2;
    JB[25] = -1.0*dwdx2;
    JB[36] = -1.0*dwdx1 - 1.0*dwdx2;
    JB[46] = 1.0*dwdx3 + 1.0*dwdx4;
    JB[47] = 1.0*dwdx3;
    JB[48] = -1.0*dwdx3;
    JB[50] = 1.0*dwdx4;
    JB[51] = -1.0*dwdx4;
    JB[68] = 1.0*dwdx5;
    JB[69] = 1.0*dwdx5 + 1.0*dwdx6;
    JB[70] = -1.0*dwdx5;
    JB[71] = -1.0*dwdx6;
    JB[80] = 1.0*dwdx6;
    JB[88] = -1.0*dwdx8;
    JB[90] = -1.0*dwdx7 - 1.0*dwdx8;
    JB[91] = -1.0*dwdx7;
    JB[92] = 1.0*dwdx7 + 1.0*dwdx8;
    JB[113] = -1.0*dwdx9;
    JB[115] = 1.0*dwdx10 + 1.0*dwdx9;
    JB[116] = -1.0*dwdx10;
    JB[124] = -1.0*dwdx10 - 1.0*dwdx9;
    JB[134] = 1.0*dwdx11;
    JB[138] = 1.0*dwdx11;
    JB[139] = -1.0*dwdx11;
    JB[156] = -1.0*dwdx12 - 1.0*dwdx13;
    JB[157] = -1.0*dwdx13;
    JB[160] = -1.0*dwdx12;
    JB[161] = 1.0*dwdx12 + 1.0*dwdx13;
    JB[184] = 1.0*dwdx14;
    JB[186] = -1.0*dwdx14;
    JB[196] = 1.0*dwdx14;
    JB[207] = 1.0*dwdx15 + 1.0*dwdx16;
    JB[209] = 1.0*dwdx15;
    JB[210] = -1.0*dwdx15;
    JB[212] = 1.0*dwdx16;
    JB[213] = -1.0*dwdx16;
    JB[228] = -1.0*dwdx17;
    JB[230] = 1.0*dwdx17 + 1.0*dwdx18;
    JB[231] = -1.0*dwdx18;
    JB[240] = -1.0*dwdx17 - 1.0*dwdx18;
    JB[251] = 1.0*dwdx19;
    JB[253] = 1.0*dwdx19 + 1.0*dwdx20;
    JB[254] = -1.0*dwdx19;
    JB[255] = -1.0*dwdx20;
    JB[262] = 1.0*dwdx20;
    JB[272] = -1.0*dwdx22;
    JB[273] = -1.0*dwdx21 - 1.0*dwdx22;
    JB[275] = -1.0*dwdx21;
    JB[276] = 1.0*dwdx21 + 1.0*dwdx22;
    JB[297] = -1.0*dwdx23;
    JB[299] = 1.0*dwdx23 + 1.0*dwdx24;
    JB[300] = -1.0*dwdx24;
    JB[306] = -1.0*dwdx23 - 1.0*dwdx24;
    JB[308] = 1.0*dwdx26;
    JB[309] = -1.0*dwdx26;
    JB[311] = 1.0*dwdx27;
    JB[313] = -1.0*dwdx27;
    JB[317] = 1.0*dwdx25;
    JB[322] = 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27;
    JB[323] = -1.0*dwdx25;
    JB[339] = -1.0*dwdx28 - 1.0*dwdx29;
    JB[341] = -1.0*dwdx29;
    JB[344] = -1.0*dwdx28;
    JB[345] = 1.0*dwdx28 + 1.0*dwdx29;
    JB[368] = 1.0*dwdx30;
    JB[369] = 1.0*dwdx30;
    JB[371] = -1.0*dwdx30;
    JB[390] = 1.0*dwdx31;
    JB[391] = 1.0*dwdx31;
    JB[393] = -1.0*dwdx31;
    JB[414] = 1.0*dwdx32;
    JB[416] = 1.0*dwdx32;
    JB[417] = -1.0*dwdx32;
    JB[434] = -1.0*dwdx33;
    JB[435] = -1.0*dwdx33 - 1.0*dwdx34;
    JB[437] = 1.0*dwdx33 + 1.0*dwdx34;
    JB[438] = -1.0*dwdx34;
    JB[448] = 1.0*dwdx37;
    JB[450] = -1.0*dwdx37;
    JB[451] = 1.0*dwdx35;
    JB[453] = -1.0*dwdx35;
    JB[458] = 1.0*dwdx36;
    JB[460] = 1.0*dwdx35 + 1.0*dwdx36 + 1.0*dwdx37;
    JB[461] = -1.0*dwdx36;
    JB[478] = -1.0*dwdx39;
    JB[480] = -1.0*dwdx38 - 1.0*dwdx39;
    JB[482] = -1.0*dwdx38;
    JB[483] = 1.0*dwdx38 + 1.0*dwdx39;
}