#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model5_levchenko2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[1] = 1.0*dwdx5;
    J[2] = 1.0*dwdx11;
    J[3] = 1.0*dwdx16;
    J[4] = 1.0*dwdx21;
    J[9] = -1.0*dwdx50;
    J[17] = -1.0*dwdx71;
    J[31] = 1.0*dwdx1;
    J[32] = -1.0*dwdx6;
    J[36] = 1.0*dwdx26;
    J[39] = 1.0*dwdx43;
    J[40] = -1.0*dwdx51;
    J[48] = 1.0*dwdx72;
    J[60] = -1.0*dwdx101;
    J[63] = 1.0*dwdx7;
    J[64] = -1.0*dwdx12;
    J[68] = 1.0*dwdx32;
    J[69] = 1.0*dwdx38;
    J[71] = -1.0*dwdx52;
    J[91] = 1.0*dwdx102;
    J[93] = 1.0*dwdx2;
    J[96] = -1.0*dwdx17;
    J[98] = 1.0*dwdx27;
    J[99] = 1.0*dwdx33;
    J[102] = 1.0*dwdx53;
    J[110] = -1.0*dwdx73;
    J[128] = -1.0*dwdx22;
    J[131] = 1.0*dwdx39;
    J[132] = 1.0*dwdx44;
    J[141] = -1.0*dwdx74;
    J[156] = 1.0*dwdx8;
    J[158] = 1.0*dwdx18;
    J[160] = -1.0*dwdx28;
    J[164] = 1.0*dwdx54;
    J[172] = 1.0*dwdx75;
    J[184] = -1.0*dwdx103;
    J[188] = 1.0*dwdx13;
    J[191] = 1.0*dwdx29;
    J[192] = -1.0*dwdx34;
    J[195] = 1.0*dwdx55;
    J[215] = 1.0*dwdx104;
    J[223] = 1.0*dwdx35;
    J[224] = -1.0*dwdx40;
    J[225] = 1.0*dwdx45;
    J[246] = 1.0*dwdx105;
    J[252] = 1.0*dwdx23;
    J[256] = -1.0*dwdx46;
    J[265] = 1.0*dwdx76;
    J[277] = -1.0*dwdx106;
    J[279] = -1.0*dwdx3;
    J[280] = -1.0*dwdx9;
    J[281] = -1.0*dwdx14;
    J[282] = 1.0*dwdx19;
    J[284] = 1.0*dwdx30;
    J[285] = 1.0*dwdx36;
    J[288] = -1.0*dwdx49 - 1.0*dwdx56;
    J[289] = 1.0*dwdx57;
    J[292] = 1.0*dwdx64;
    J[302] = -1.0*dwdx89;
    J[319] = 1.0*dwdx49;
    J[320] = -1.0*dwdx57 - 1.0*dwdx58;
    J[333] = 1.0*dwdx89;
    J[352] = -1.0*dwdx59 - 1.0*dwdx60;
    J[353] = -1.0*dwdx61;
    J[354] = 1.0*dwdx63 + 1.0*dwdx64;
    J[356] = -1.0*dwdx67;
    J[357] = 1.0*dwdx68 + 1.0*dwdx69;
    J[382] = 1.0*dwdx58;
    J[383] = -1.0*dwdx59;
    J[384] = -1.0*dwdx61 - 1.0*dwdx62;
    J[385] = 1.0*dwdx63;
    J[386] = 1.0*dwdx65;
    J[388] = 1.0*dwdx69;
    J[395] = -1.0*dwdx90;
    J[414] = 1.0*dwdx59;
    J[415] = 1.0*dwdx61;
    J[416] = -1.0*dwdx63 - 1.0*dwdx64;
    J[446] = 1.0*dwdx62;
    J[448] = -1.0*dwdx65 - 1.0*dwdx66;
    J[457] = 1.0*dwdx90;
    J[469] = 1.0*dwdx24;
    J[472] = 1.0*dwdx41;
    J[473] = 1.0*dwdx47;
    J[476] = -1.0*dwdx60;
    J[479] = 1.0*dwdx66;
    J[480] = -1.0*dwdx67;
    J[481] = 1.0*dwdx68;
    J[507] = 1.0*dwdx60;
    J[511] = 1.0*dwdx67;
    J[512] = -1.0*dwdx68 - 1.0*dwdx69;
    J[527] = -1.0*dwdx4;
    J[528] = 1.0*dwdx10;
    J[530] = -1.0*dwdx20;
    J[531] = -1.0*dwdx25;
    J[532] = 1.0*dwdx31;
    J[535] = 1.0*dwdx48;
    J[544] = -1.0*dwdx70 - 1.0*dwdx77;
    J[546] = 1.0*dwdx80;
    J[548] = 1.0*dwdx85;
    J[556] = -1.0*dwdx100;
    J[576] = -1.0*dwdx78 - 1.0*dwdx79;
    J[578] = -1.0*dwdx82;
    J[579] = 1.0*dwdx84 + 1.0*dwdx85;
    J[581] = -1.0*dwdx88;
    J[582] = 1.0*dwdx91 + 1.0*dwdx92;
    J[606] = 1.0*dwdx70;
    J[608] = -1.0*dwdx80 - 1.0*dwdx81;
    J[618] = 1.0*dwdx100;
    J[638] = -1.0*dwdx78;
    J[639] = 1.0*dwdx81;
    J[640] = -1.0*dwdx82 - 1.0*dwdx83;
    J[641] = 1.0*dwdx84;
    J[642] = 1.0*dwdx86;
    J[644] = 1.0*dwdx92;
    J[649] = -1.0*dwdx98;
    J[669] = 1.0*dwdx78;
    J[671] = 1.0*dwdx82;
    J[672] = -1.0*dwdx84 - 1.0*dwdx85;
    J[702] = 1.0*dwdx83;
    J[704] = -1.0*dwdx86 - 1.0*dwdx87;
    J[711] = 1.0*dwdx98;
    J[715] = 1.0*dwdx15;
    J[719] = 1.0*dwdx37;
    J[720] = 1.0*dwdx42;
    J[722] = -1.0*dwdx49;
    J[723] = 1.0*dwdx57 + 1.0*dwdx58;
    J[725] = -1.0*dwdx62;
    J[727] = 1.0*dwdx65 + 1.0*dwdx66;
    J[731] = -1.0*dwdx79;
    J[735] = 1.0*dwdx87;
    J[736] = -1.0*dwdx88 - 1.0*dwdx89 - 1.0*dwdx90;
    J[737] = 1.0*dwdx91;
    J[762] = 1.0*dwdx79;
    J[767] = 1.0*dwdx88;
    J[768] = -1.0*dwdx91 - 1.0*dwdx92;
    J[800] = -1.0*dwdx93;
    J[801] = -1.0*dwdx94;
    J[803] = 1.0*dwdx96;
    J[805] = 1.0*dwdx108;
    J[831] = -1.0*dwdx93;
    J[832] = -1.0*dwdx94;
    J[834] = 1.0*dwdx96 + 1.0*dwdx97;
    J[864] = -1.0*dwdx95;
    J[866] = -1.0*dwdx99;
    J[867] = 1.0*dwdx107 + 1.0*dwdx108;
    J[893] = 1.0*dwdx93;
    J[894] = 1.0*dwdx94;
    J[896] = -1.0*dwdx96 - 1.0*dwdx97;
    J[916] = -1.0*dwdx70;
    J[918] = 1.0*dwdx80 + 1.0*dwdx81;
    J[919] = -1.0*dwdx83;
    J[921] = 1.0*dwdx86 + 1.0*dwdx87;
    J[926] = -1.0*dwdx95;
    J[927] = 1.0*dwdx97;
    J[928] = -1.0*dwdx100 - 1.0*dwdx98 - 1.0*dwdx99;
    J[929] = 1.0*dwdx107;
    J[957] = 1.0*dwdx95;
    J[959] = 1.0*dwdx99;
    J[960] = -1.0*dwdx107 - 1.0*dwdx108;
}