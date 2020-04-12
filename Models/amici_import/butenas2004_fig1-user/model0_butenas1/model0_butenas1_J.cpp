#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_butenas1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4;
    J[2] = -1.0*dwdx9;
    J[5] = -1.0*dwdx13;
    J[12] = -1.0*dwdx28;
    J[27] = -1.0*dwdx58;
    J[32] = -1.0*dwdx70;
    J[35] = -1.0*dwdx5 - 1.0*dwdx6;
    J[61] = -1.0*dwdx61;
    J[64] = -1.0*dwdx65;
    J[65] = -1.0*dwdx67;
    J[68] = -1.0*dwdx3;
    J[69] = 1.0*dwdx6;
    J[70] = -1.0*dwdx9;
    J[95] = 1.0*dwdx61;
    J[98] = 1.0*dwdx66;
    J[100] = 1.0*dwdx69;
    J[102] = 1.0*dwdx3;
    J[104] = 1.0*dwdx9;
    J[140] = -1.0*dwdx11;
    J[148] = -1.0*dwdx32;
    J[150] = -1.0*dwdx33;
    J[170] = -1.0*dwdx2;
    J[175] = -1.0*dwdx12 - 1.0*dwdx13;
    J[177] = -1.0*dwdx15 + 1.0*dwdx17;
    J[178] = 1.0*dwdx20;
    J[184] = 1.0*dwdx34;
    J[191] = -1.0*dwdx47;
    J[204] = 1.0*dwdx2;
    J[209] = 1.0*dwdx13;
    J[243] = 1.0*dwdx12;
    J[245] = 1.0*dwdx15 - 1.0*dwdx16 - 1.0*dwdx17;
    J[246] = -1.0*dwdx18 + 1.0*dwdx19;
    J[259] = 1.0*dwdx47;
    J[264] = -1.0*dwdx53;
    J[279] = 1.0*dwdx16;
    J[280] = 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20;
    J[298] = 1.0*dwdx53;
    J[315] = -1.0*dwdx21 - 1.0*dwdx22;
    J[317] = -1.0*dwdx25;
    J[318] = -1.0*dwdx26;
    J[325] = -1.0*dwdx42;
    J[330] = -1.0*dwdx51;
    J[350] = -1.0*dwdx23 - 1.0*dwdx24;
    J[356] = -1.0*dwdx37;
    J[357] = -1.0*dwdx39;
    J[367] = -1.0*dwdx57;
    J[369] = -1.0*dwdx62;
    J[383] = 1.0*dwdx21;
    J[385] = 1.0*dwdx25;
    J[393] = 1.0*dwdx42;
    J[408] = -1.0*dwdx4;
    J[412] = -1.0*dwdx11;
    J[417] = 1.0*dwdx22;
    J[420] = 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx30 - 1.0*dwdx31 - 1.0*dwdx32;
    J[422] = -1.0*dwdx33 + 1.0*dwdx34;
    J[423] = -1.0*dwdx35;
    J[424] = -1.0*dwdx38;
    J[432] = 1.0*dwdx51;
    J[434] = -1.0*dwdx55;
    J[435] = -1.0*dwdx60;
    J[437] = -1.0*dwdx63;
    J[442] = 1.0*dwdx4;
    J[454] = 1.0*dwdx28;
    J[480] = 1.0*dwdx11;
    J[488] = 1.0*dwdx32;
    J[490] = 1.0*dwdx33 - 1.0*dwdx34;
    J[522] = 1.0*dwdx30;
    J[525] = 1.0*dwdx35 - 1.0*dwdx36;
    J[536] = 1.0*dwdx55;
    J[554] = -1.0*dwdx24;
    J[556] = 1.0*dwdx31;
    J[559] = 1.0*dwdx36;
    J[560] = -1.0*dwdx37 + 1.0*dwdx38;
    J[561] = -1.0*dwdx39;
    J[571] = 1.0*dwdx60;
    J[588] = 1.0*dwdx24;
    J[590] = 1.0*dwdx27;
    J[594] = 1.0*dwdx37;
    J[595] = 1.0*dwdx39;
    J[607] = 1.0*dwdx63;
    J[614] = -1.0*dwdx8;
    J[630] = -1.0*dwdx40 - 1.0*dwdx41;
    J[644] = -1.0*dwdx71;
    J[648] = -1.0*dwdx10;
    J[655] = -1.0*dwdx21;
    J[657] = -1.0*dwdx25;
    J[658] = -1.0*dwdx29;
    J[665] = -1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44 - 1.0*dwdx45;
    J[673] = -1.0*dwdx59;
    J[682] = -1.0*dwdx7;
    J[700] = -1.0*dwdx46;
    J[716] = 1.0*dwdx7;
    J[719] = -1.0*dwdx12;
    J[721] = -1.0*dwdx15;
    J[734] = 1.0*dwdx46;
    J[735] = -1.0*dwdx47 - 1.0*dwdx48;
    J[736] = -1.0*dwdx49;
    J[737] = -1.0*dwdx50;
    J[755] = 1.0*dwdx17;
    J[756] = 1.0*dwdx20;
    J[769] = 1.0*dwdx48;
    J[770] = 1.0*dwdx49;
    J[771] = 1.0*dwdx50;
    J[789] = 1.0*dwdx17;
    J[790] = 1.0*dwdx20;
    J[803] = 1.0*dwdx48;
    J[804] = 1.0*dwdx49;
    J[805] = 1.0*dwdx50;
    J[818] = 1.0*dwdx10;
    J[825] = -1.0*dwdx22;
    J[828] = -1.0*dwdx26 + 1.0*dwdx29;
    J[835] = 1.0*dwdx43 + 1.0*dwdx44 + 1.0*dwdx45;
    J[840] = -1.0*dwdx51;
    J[843] = 1.0*dwdx59;
    J[852] = 1.0*dwdx8;
    J[868] = 1.0*dwdx40 + 1.0*dwdx41;
    J[875] = -1.0*dwdx52;
    J[877] = -1.0*dwdx56;
    J[880] = -1.0*dwdx64;
    J[882] = 1.0*dwdx71;
    J[889] = -1.0*dwdx14;
    J[891] = -1.0*dwdx16;
    J[892] = -1.0*dwdx18 + 1.0*dwdx20;
    J[896] = -1.0*dwdx30;
    J[899] = -1.0*dwdx35;
    J[910] = -1.0*dwdx53 - 1.0*dwdx54 - 1.0*dwdx55;
    J[918] = -1.0*dwdx0;
    J[923] = 1.0*dwdx14;
    J[926] = 1.0*dwdx19;
    J[928] = -1.0*dwdx23;
    J[930] = -1.0*dwdx31;
    J[934] = -1.0*dwdx38;
    J[943] = -1.0*dwdx52;
    J[944] = 1.0*dwdx54;
    J[945] = -1.0*dwdx56 - 1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx60;
    J[947] = -1.0*dwdx62;
    J[948] = -1.0*dwdx64;
    J[952] = 1.0*dwdx0;
    J[979] = 1.0*dwdx58;
    J[996] = 1.0*dwdx23;
    J[998] = -1.0*dwdx27;
    J[1013] = 1.0*dwdx57;
    J[1015] = 1.0*dwdx62 - 1.0*dwdx63;
    J[1021] = -1.0*dwdx5;
    J[1045] = 1.0*dwdx52;
    J[1047] = 1.0*dwdx56;
    J[1050] = 1.0*dwdx64 - 1.0*dwdx65;
    J[1051] = -1.0*dwdx67 + 1.0*dwdx68;
    J[1055] = 1.0*dwdx5;
    J[1084] = 1.0*dwdx65;
    J[1085] = 1.0*dwdx67 - 1.0*dwdx68;
    J[1088] = -1.0*dwdx1;
    J[1118] = -1.0*dwdx66;
    J[1119] = 1.0*dwdx68;
    J[1120] = -1.0*dwdx69 - 1.0*dwdx70;
    J[1122] = 1.0*dwdx1;
    J[1154] = 1.0*dwdx70;
}