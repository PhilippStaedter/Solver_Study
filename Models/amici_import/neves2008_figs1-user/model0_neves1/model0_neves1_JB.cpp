#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_neves1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[38] = -5.0*dwdx0;
    JB[39] = 5.0*dwdx0;
    JB[41] = 1.0*dwdx1;
    JB[49] = 1.0*dwdx0;
    JB[68] = -1.0*dwdx1;
    JB[75] = -5.0*dwdx2;
    JB[76] = 5.0*dwdx2;
    JB[78] = 1.0*dwdx3;
    JB[86] = 1.0*dwdx2;
    JB[105] = -1.0*dwdx3;
    JB[152] = 1.0*dwdx4 + 1.0*dwdx5;
    JB[179] = -1.0*dwdx4 - 1.0*dwdx5;
    JB[190] = -5.0*dwdx6 + 5.0*dwdx7;
    JB[191] = 5.0*dwdx6;
    JB[198] = 1.0*dwdx6;
    JB[218] = -5.0*dwdx7;
    JB[221] = 9.0000000000000089*dwdx7;
    JB[227] = -5.0*dwdx8;
    JB[228] = 5.0*dwdx8 + 5.0*dwdx9;
    JB[235] = 1.0*dwdx8;
    JB[256] = -5.0*dwdx9;
    JB[258] = 9.0000000000000089*dwdx9;
    JB[266] = 1.0*dwdx10;
    JB[267] = -1.0*dwdx10;
    JB[275] = -1.0*dwdx11;
    JB[276] = 1.0*dwdx11;
    JB[303] = -1.0*dwdx12;
    JB[304] = 1.0*dwdx12;
    JB[342] = -1.0*dwdx14;
    JB[343] = 1.0*dwdx14;
    JB[361] = 1.0*dwdx14;
    JB[367] = 5.0*dwdx13;
    JB[368] = -5.0*dwdx13;
    JB[379] = -1.0*dwdx16;
    JB[380] = 1.0*dwdx16;
    JB[398] = 1.0*dwdx16;
    JB[404] = 5.0*dwdx15;
    JB[405] = -5.0*dwdx15;
    JB[418] = -1.0*dwdx17 + 1.0*dwdx18;
    JB[419] = 1.0*dwdx17;
    JB[420] = -1.0*dwdx18;
    JB[435] = 1.0*dwdx18;
    JB[445] = -5.0*dwdx19;
    JB[446] = 5.0*dwdx19;
    JB[455] = -1.0*dwdx20;
    JB[456] = 1.0*dwdx19 + 1.0*dwdx20 - 1.0*dwdx21;
    JB[472] = -1.0*dwdx21;
    JB[477] = 5.0*dwdx21;
    JB[478] = -5.0*dwdx21;
    JB[486] = -5.0*dwdx22;
    JB[487] = 5.0*dwdx22;
    JB[492] = 1.0*dwdx24;
    JB[494] = 1.0*dwdx22 + 1.0*dwdx23 - 1.0*dwdx24;
    JB[509] = 1.0*dwdx24;
    JB[514] = -5.0*dwdx23;
    JB[515] = 5.0*dwdx23;
    JB[532] = 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27;
    JB[533] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27;
    JB[569] = -1.0*dwdx28;
    JB[570] = 1.0*dwdx28;
    JB[606] = -1.0*dwdx29;
    JB[607] = 1.0*dwdx29;
    JB[608] = 1.0*dwdx30;
    JB[609] = -1.0*dwdx30;
    JB[645] = -1.0*dwdx31;
    JB[646] = 1.0*dwdx31;
    JB[669] = -1.0*dwdx32;
    JB[684] = 1.0*dwdx33;
    JB[685] = -1.0*dwdx33;
    JB[697] = 1.0*dwdx32;
    JB[706] = -1.0*dwdx34;
    JB[721] = -1.0*dwdx35;
    JB[722] = 1.0*dwdx35;
    JB[734] = 1.0*dwdx34;
    JB[743] = -1.0*dwdx36;
    JB[771] = 1.0*dwdx36;
    JB[784] = -1.0*dwdx40;
    JB[785] = 1.0*dwdx40;
    JB[795] = -1.0*dwdx38;
    JB[796] = 1.0*dwdx38;
    JB[798] = -1.0*dwdx37;
    JB[801] = -1.0*dwdx39;
    JB[803] = 1.0*dwdx39;
    JB[807] = 1.0*dwdx37;
    JB[808] = 1.0*dwdx37;
    JB[821] = 1.0*dwdx42;
    JB[822] = -1.0*dwdx42;
    JB[828] = 1.0*dwdx41;
    JB[829] = -1.0*dwdx41;
    JB[830] = 1.0*dwdx43;
    JB[831] = -1.0*dwdx43;
    JB[869] = 1.0*dwdx44;
    JB[870] = -1.0*dwdx44;
    JB[902] = 1.0*dwdx45;
    JB[903] = -1.0*dwdx45;
    JB[912] = 1.0*dwdx46;
    JB[914] = -1.0*dwdx46;
    JB[949] = 1.0*dwdx47;
    JB[951] = -1.0*dwdx47;
    JB[976] = 1.0*dwdx49;
    JB[977] = -1.0*dwdx49;
    JB[986] = -1.0*dwdx48;
    JB[988] = 1.0*dwdx48;
    JB[1026] = 1.0*dwdx50;
    JB[1030] = 1.0*dwdx50;
    JB[1031] = -1.0*dwdx50;
    JB[1045] = -1.0*dwdx52;
    JB[1046] = 1.0*dwdx52;
    JB[1047] = 1.0*dwdx53;
    JB[1048] = -1.0*dwdx51;
    JB[1049] = -1.0*dwdx53;
    JB[1064] = -1.0*dwdx51 + 1.0*dwdx52 + 1.0*dwdx53;
    JB[1069] = 5.0*dwdx51;
    JB[1070] = -5.0*dwdx51;
    JB[1102] = 1.0*dwdx54 - 1.0*dwdx55;
    JB[1103] = -1.0*dwdx54;
    JB[1104] = 1.0*dwdx54 + 1.0*dwdx55;
    JB[1105] = 1.0*dwdx55;
    JB[1131] = -1.0*dwdx57;
    JB[1139] = 1.0*dwdx56;
    JB[1140] = -1.0*dwdx56 + 1.0*dwdx57;
    JB[1141] = 1.0*dwdx56 + 1.0*dwdx57;
    JB[1150] = -1.0*dwdx62 - 1.0*dwdx63 - 1.0*dwdx64;
    JB[1168] = -1.0*dwdx59;
    JB[1174] = 1.0*dwdx60;
    JB[1176] = 1.0*dwdx58 - 1.0*dwdx61;
    JB[1177] = -1.0*dwdx58 + 1.0*dwdx59;
    JB[1178] = 1.0*dwdx58 + 1.0*dwdx59 + 1.0*dwdx60 + 1.0*dwdx61 + 1.0*dwdx62 + 1.0*dwdx63 + 1.0*dwdx64;
    JB[1179] = -1.0*dwdx60 + 1.0*dwdx61;
    JB[1211] = 1.0*dwdx65;
    JB[1213] = -1.0*dwdx66;
    JB[1215] = 1.0*dwdx65 + 1.0*dwdx66;
    JB[1216] = -1.0*dwdx65 + 1.0*dwdx66;
    JB[1226] = 5.0*dwdx69;
    JB[1233] = -1.0*dwdx68;
    JB[1234] = 1.0*dwdx67;
    JB[1249] = -1.0*dwdx68;
    JB[1254] = -5.0*dwdx67 + 5.0*dwdx68 - 5.0*dwdx69;
    JB[1255] = 5.0*dwdx67 - 5.0*dwdx68;
    JB[1257] = 9.0000000000000089*dwdx69;
    JB[1264] = 5.0*dwdx74;
    JB[1270] = -1.0*dwdx73;
    JB[1271] = 1.0*dwdx72;
    JB[1286] = -1.0*dwdx73;
    JB[1291] = -5.0*dwdx72 + 5.0*dwdx73;
    JB[1292] = 5.0*dwdx70 + 5.0*dwdx71 + 5.0*dwdx72 - 5.0*dwdx73 - 5.0*dwdx74;
    JB[1293] = -5.0*dwdx70 - 5.0*dwdx71;
    JB[1294] = 9.0000000000000089*dwdx74;
    JB[1337] = 5.0*dwdx76;
    JB[1338] = 5.0*dwdx75;
    JB[1365] = -5.0*dwdx76;
    JB[1366] = -5.0*dwdx75;
    JB[1368] = 9.0000000000000089*dwdx75 + 9.0000000000000089*dwdx76;
}