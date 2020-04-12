#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model6_levchenko1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[1] = 1.0*dwdx1;
    J[4] = 1.0*dwdx8;
    J[14] = -1.0*dwdx26;
    J[22] = 1.0*dwdx0;
    J[23] = -1.0*dwdx1 - 1.0*dwdx2;
    J[36] = 1.0*dwdx26;
    J[46] = -1.0*dwdx3 - 1.0*dwdx4;
    J[47] = -1.0*dwdx5;
    J[48] = 1.0*dwdx7 + 1.0*dwdx8;
    J[50] = -1.0*dwdx11;
    J[51] = 1.0*dwdx12 + 1.0*dwdx13;
    J[67] = 1.0*dwdx2;
    J[68] = -1.0*dwdx3;
    J[69] = -1.0*dwdx5 - 1.0*dwdx6;
    J[70] = 1.0*dwdx7;
    J[71] = 1.0*dwdx9;
    J[73] = 1.0*dwdx13;
    J[80] = -1.0*dwdx27;
    J[90] = 1.0*dwdx3;
    J[91] = 1.0*dwdx5;
    J[92] = -1.0*dwdx7 - 1.0*dwdx8;
    J[113] = 1.0*dwdx6;
    J[115] = -1.0*dwdx10 - 1.0*dwdx9;
    J[124] = 1.0*dwdx27;
    J[134] = -1.0*dwdx4;
    J[137] = 1.0*dwdx10;
    J[138] = -1.0*dwdx11;
    J[139] = 1.0*dwdx12;
    J[156] = 1.0*dwdx4;
    J[160] = 1.0*dwdx11;
    J[161] = -1.0*dwdx12 - 1.0*dwdx13;
    J[184] = -1.0*dwdx14;
    J[186] = 1.0*dwdx17;
    J[188] = 1.0*dwdx22;
    J[196] = -1.0*dwdx37;
    J[207] = -1.0*dwdx15 - 1.0*dwdx16;
    J[209] = -1.0*dwdx19;
    J[210] = 1.0*dwdx21 + 1.0*dwdx22;
    J[212] = -1.0*dwdx25;
    J[213] = 1.0*dwdx28 + 1.0*dwdx29;
    J[228] = 1.0*dwdx14;
    J[230] = -1.0*dwdx17 - 1.0*dwdx18;
    J[240] = 1.0*dwdx37;
    J[251] = -1.0*dwdx15;
    J[252] = 1.0*dwdx18;
    J[253] = -1.0*dwdx19 - 1.0*dwdx20;
    J[254] = 1.0*dwdx21;
    J[255] = 1.0*dwdx23;
    J[257] = 1.0*dwdx29;
    J[262] = -1.0*dwdx35;
    J[273] = 1.0*dwdx15;
    J[275] = 1.0*dwdx19;
    J[276] = -1.0*dwdx21 - 1.0*dwdx22;
    J[297] = 1.0*dwdx20;
    J[299] = -1.0*dwdx23 - 1.0*dwdx24;
    J[306] = 1.0*dwdx35;
    J[308] = -1.0*dwdx0;
    J[309] = 1.0*dwdx1 + 1.0*dwdx2;
    J[311] = -1.0*dwdx6;
    J[313] = 1.0*dwdx10 + 1.0*dwdx9;
    J[317] = -1.0*dwdx16;
    J[321] = 1.0*dwdx24;
    J[322] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27;
    J[323] = 1.0*dwdx28;
    J[339] = 1.0*dwdx16;
    J[344] = 1.0*dwdx25;
    J[345] = -1.0*dwdx28 - 1.0*dwdx29;
    J[368] = -1.0*dwdx30;
    J[369] = -1.0*dwdx31;
    J[371] = 1.0*dwdx33;
    J[373] = 1.0*dwdx39;
    J[390] = -1.0*dwdx30;
    J[391] = -1.0*dwdx31;
    J[393] = 1.0*dwdx33 + 1.0*dwdx34;
    J[414] = -1.0*dwdx32;
    J[416] = -1.0*dwdx36;
    J[417] = 1.0*dwdx38 + 1.0*dwdx39;
    J[434] = 1.0*dwdx30;
    J[435] = 1.0*dwdx31;
    J[437] = -1.0*dwdx33 - 1.0*dwdx34;
    J[448] = -1.0*dwdx14;
    J[450] = 1.0*dwdx17 + 1.0*dwdx18;
    J[451] = -1.0*dwdx20;
    J[453] = 1.0*dwdx23 + 1.0*dwdx24;
    J[458] = -1.0*dwdx32;
    J[459] = 1.0*dwdx34;
    J[460] = -1.0*dwdx35 - 1.0*dwdx36 - 1.0*dwdx37;
    J[461] = 1.0*dwdx38;
    J[480] = 1.0*dwdx32;
    J[482] = 1.0*dwdx36;
    J[483] = -1.0*dwdx38 - 1.0*dwdx39;
}