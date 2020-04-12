#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_Yang2007(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx5 - 1.0*dwdx6;
    J[1] = -1.0*dwdx7;
    J[2] = 1.0*dwdx11 - 1.0*dwdx12;
    J[4] = 1.0*dwdx22;
    J[6] = 1.0*dwdx28;
    J[7] = -1.0*dwdx29;
    J[8] = -1.0*dwdx31;
    J[9] = -1.0*dwdx33;
    J[11] = 1.0*dwdx35 - 1.0*dwdx36;
    J[13] = -1.0*dwdx45;
    J[18] = -1.0*dwdx58 - 1.0*dwdx60;
    J[19] = 1.0*dwdx63 - 1.0*dwdx64;
    J[20] = -1.0*dwdx67;
    J[21] = -1.0*dwdx70;
    J[22] = -1.0*dwdx74 - 1.0*dwdx76;
    J[25] = 1.0*dwdx5;
    J[26] = 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    J[27] = 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14;
    J[28] = -1.0*dwdx17;
    J[38] = 1.0*dwdx45 - 1.0*dwdx46;
    J[41] = -1.0*dwdx56;
    J[43] = 1.0*dwdx60 - 1.0*dwdx61;
    J[45] = 1.0*dwdx67 - 1.0*dwdx68;
    J[47] = 1.0*dwdx76 - 1.0*dwdx77;
    J[51] = 1.0*dwdx8;
    J[52] = 1.0*dwdx13 - 1.0*dwdx16;
    J[66] = 1.0*dwdx56;
    J[76] = 1.0*dwdx9;
    J[77] = 1.0*dwdx14;
    J[78] = 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx20;
    J[79] = -1.0*dwdx23;
    J[88] = 1.0*dwdx46;
    J[89] = -1.0*dwdx51;
    J[93] = 1.0*dwdx61;
    J[95] = 1.0*dwdx68;
    J[97] = 1.0*dwdx77;
    J[102] = -1.0*dwdx15;
    J[103] = 1.0*dwdx18;
    J[104] = 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx26;
    J[105] = -1.0*dwdx27;
    J[114] = 1.0*dwdx51;
    J[115] = -1.0*dwdx53;
    J[120] = -1.0*dwdx69;
    J[127] = 1.0*dwdx15;
    J[129] = 1.0*dwdx24;
    J[130] = 1.0*dwdx27;
    J[140] = 1.0*dwdx53;
    J[145] = 1.0*dwdx69;
    J[182] = -1.0*dwdx30;
    J[197] = 1.0*dwdx78;
    J[208] = -1.0*dwdx32;
    J[211] = -1.0*dwdx38;
    J[275] = 1.0*dwdx1;
    J[282] = 1.0*dwdx29;
    J[286] = 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx40;
    J[291] = -1.0*dwdx54;
    J[293] = -1.0*dwdx57;
    J[311] = -1.0*dwdx39;
    J[312] = -1.0*dwdx43 - 1.0*dwdx44;
    J[321] = -1.0*dwdx73;
    J[326] = -1.0*dwdx10;
    J[328] = -1.0*dwdx21;
    J[329] = 1.0*dwdx25;
    J[336] = -1.0*dwdx41;
    J[338] = 1.0*dwdx47 - 1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50;
    J[353] = -1.0*dwdx19;
    J[364] = -1.0*dwdx52;
    J[461] = 1.0*dwdx37;
    J[466] = 1.0*dwdx54;
    J[468] = 1.0*dwdx57 - 1.0*dwdx62;
    J[475] = 1.0*dwdx2;
    J[483] = 1.0*dwdx31;
    J[491] = -1.0*dwdx55;
    J[493] = 1.0*dwdx58;
    J[494] = 1.0*dwdx64 - 1.0*dwdx65;
    J[495] = -1.0*dwdx66;
    J[516] = 1.0*dwdx55;
    J[519] = 1.0*dwdx65;
    J[520] = 1.0*dwdx66;
    J[525] = 1.0*dwdx3 - 1.0*dwdx4;
    J[534] = 1.0*dwdx33;
    J[535] = -1.0*dwdx34;
    J[537] = -1.0*dwdx42;
    J[543] = -1.0*dwdx59;
    J[546] = 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72;
    J[547] = 1.0*dwdx74 - 1.0*dwdx75;
    J[548] = -1.0*dwdx79;
    J[550] = 1.0*dwdx4;
    J[560] = 1.0*dwdx34;
    J[568] = 1.0*dwdx59;
    J[571] = 1.0*dwdx71;
    J[572] = 1.0*dwdx75;
    J[587] = 1.0*dwdx42;
    J[596] = 1.0*dwdx72;
    J[598] = 1.0*dwdx79 - 1.0*dwdx81;
    J[623] = 1.0*dwdx80;
    J[624] = -1.0*dwdx82;
}