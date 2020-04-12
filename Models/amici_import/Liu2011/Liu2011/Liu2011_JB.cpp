#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_Liu2011(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JB[1] = 1.0*dwdx0;
    JB[2] = -1.0*dwdx0;
    JB[19] = 1.0*dwdx1;
    JB[24] = -1.0*dwdx1;
    JB[42] = 1.0*dwdx2;
    JB[43] = 1.0*dwdx2;
    JB[44] = -1.0*dwdx2;
    JB[84] = 1.0*dwdx3;
    JB[85] = 1.0*dwdx3;
    JB[86] = -1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6;
    JB[91] = 1.0*dwdx4;
    JB[92] = -1.0*dwdx4;
    JB[104] = 1.0*dwdx5;
    JB[106] = -1.0*dwdx5;
    JB[110] = 1.0*dwdx6;
    JB[111] = -1.0*dwdx6;
    JB[129] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[130] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JB[131] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JB[215] = 1.0*dwdx15 + 1.0*dwdx16;
    JB[219] = 1.0*dwdx15;
    JB[221] = -1.0*dwdx15;
    JB[236] = 1.0*dwdx16;
    JB[240] = -1.0*dwdx16;
    JB[258] = 1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24;
    JB[261] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JB[262] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JB[296] = 1.0*dwdx25;
    JB[301] = 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27 + 1.0*dwdx28;
    JB[302] = -1.0*dwdx25;
    JB[316] = 1.0*dwdx27;
    JB[317] = 1.0*dwdx28;
    JB[318] = 1.0*dwdx26;
    JB[319] = -1.0*dwdx26;
    JB[327] = -1.0*dwdx27;
    JB[330] = -1.0*dwdx28;
    JB[338] = 1.0*dwdx29;
    JB[339] = 1.0*dwdx30;
    JB[340] = -1.0*dwdx30;
    JB[341] = -1.0*dwdx30;
    JB[342] = 1.0*dwdx31;
    JB[343] = 1.0*dwdx29;
    JB[344] = -1.0*dwdx29;
    JB[345] = -1.0*dwdx31;
    JB[346] = -1.0*dwdx31;
    JB[383] = 1.0*dwdx32;
    JB[387] = 1.0*dwdx32;
    JB[389] = -1.0*dwdx32;
    JB[467] = 1.0*dwdx33 - 1.0*dwdx37;
    JB[471] = 1.0*dwdx33 - 1.0*dwdx37;
    JB[473] = -1.0*dwdx33 + 1.0*dwdx35 + 1.0*dwdx36 + 1.0*dwdx37 + 1.0*dwdx38 + 1.0*dwdx39;
    JB[474] = 1.0*dwdx34;
    JB[475] = -1.0*dwdx34;
    JB[476] = -1.0*dwdx34;
    JB[479] = -1.0*dwdx35;
    JB[488] = 1.0*dwdx38;
    JB[491] = -1.0*dwdx36;
    JB[493] = -1.0*dwdx38;
    JB[516] = 1.0*dwdx40;
    JB[517] = -1.0*dwdx40;
    JB[518] = -1.0*dwdx40;
    JB[602] = 1.0*dwdx41 + 1.0*dwdx42;
    JB[603] = -1.0*dwdx41;
    JB[644] = 1.0*dwdx43;
    JB[645] = -1.0*dwdx43;
    JB[688] = 1.0*dwdx44 + 1.0*dwdx45 + 1.0*dwdx46 + 1.0*dwdx47 + 1.0*dwdx48;
    JB[691] = 1.0*dwdx45;
    JB[693] = -1.0*dwdx45;
    JB[694] = 1.0*dwdx44;
    JB[695] = -1.0*dwdx44;
    JB[696] = 1.0*dwdx46;
    JB[705] = 1.0*dwdx47;
    JB[707] = -1.0*dwdx46;
    JB[708] = -1.0*dwdx47;
    JB[710] = 1.0*dwdx48;
    JB[712] = -1.0*dwdx48;
    JB[725] = 1.0*dwdx49;
    JB[731] = -1.0*dwdx49 + 1.0*dwdx50;
    JB[740] = 1.0*dwdx50;
    JB[746] = -1.0*dwdx50;
    JB[774] = 1.0*dwdx51;
    JB[775] = -1.0*dwdx51;
    JB[776] = 1.0*dwdx51;
    JB[798] = 1.0*dwdx54;
    JB[814] = 1.0*dwdx53;
    JB[816] = 1.0*dwdx52;
    JB[817] = -1.0*dwdx52 + 1.0*dwdx53 + 1.0*dwdx54;
    JB[818] = 1.0*dwdx52;
    JB[819] = -1.0*dwdx53;
    JB[822] = -1.0*dwdx54;
    JB[842] = 1.0*dwdx55;
    JB[858] = 1.0*dwdx56;
    JB[859] = -1.0*dwdx56;
    JB[860] = 1.0*dwdx55 + 1.0*dwdx56;
    JB[862] = -1.0*dwdx55;
    JB[885] = 1.0*dwdx58;
    JB[886] = -1.0*dwdx58;
    JB[887] = -1.0*dwdx58;
    JB[888] = 1.0*dwdx59;
    JB[891] = -1.0*dwdx59;
    JB[892] = -1.0*dwdx59;
    JB[898] = 1.0*dwdx57;
    JB[901] = 1.0*dwdx57;
    JB[903] = -1.0*dwdx57;
    JB[926] = 1.0*dwdx60;
    JB[931] = 1.0*dwdx62;
    JB[940] = 1.0*dwdx61;
    JB[944] = 1.0*dwdx60;
    JB[946] = -1.0*dwdx60 + 1.0*dwdx61 + 1.0*dwdx62 + 1.0*dwdx63;
    JB[947] = -1.0*dwdx61;
    JB[950] = 1.0*dwdx63;
    JB[957] = -1.0*dwdx62;
    JB[958] = -1.0*dwdx63;
    JB[969] = 1.0*dwdx65;
    JB[970] = -1.0*dwdx65;
    JB[971] = -1.0*dwdx65;
    JB[972] = 1.0*dwdx66;
    JB[973] = 1.0*dwdx67;
    JB[975] = -1.0*dwdx66;
    JB[976] = -1.0*dwdx66;
    JB[982] = 1.0*dwdx64;
    JB[988] = 1.0*dwdx64;
    JB[989] = -1.0*dwdx64 + 1.0*dwdx67;
    JB[1002] = -1.0*dwdx67;
    JB[1008] = 1.0*dwdx68;
    JB[1015] = 1.0*dwdx69;
    JB[1024] = 1.0*dwdx71;
    JB[1027] = 1.0*dwdx68;
    JB[1032] = -1.0*dwdx68 + 1.0*dwdx69 + 1.0*dwdx70 + 1.0*dwdx71;
    JB[1033] = -1.0*dwdx69;
    JB[1034] = 1.0*dwdx70;
    JB[1036] = -1.0*dwdx70;
    JB[1043] = -1.0*dwdx71;
    JB[1053] = 1.0*dwdx73;
    JB[1054] = -1.0*dwdx73;
    JB[1055] = -1.0*dwdx73;
    JB[1056] = 1.0*dwdx74;
    JB[1057] = 1.0*dwdx72;
    JB[1059] = -1.0*dwdx74;
    JB[1060] = -1.0*dwdx74;
    JB[1074] = 1.0*dwdx72;
    JB[1075] = -1.0*dwdx72;
    JB[1094] = 1.0*dwdx75;
    JB[1097] = 1.0*dwdx78 - 1.0*dwdx79;
    JB[1101] = -1.0*dwdx79;
    JB[1103] = 1.0*dwdx77 + 1.0*dwdx79 + 1.0*dwdx80;
    JB[1109] = 1.0*dwdx81;
    JB[1114] = 1.0*dwdx83;
    JB[1116] = 1.0*dwdx76;
    JB[1118] = 1.0*dwdx75 + 1.0*dwdx76 + 1.0*dwdx78 + 1.0*dwdx80 + 1.0*dwdx81 + 1.0*dwdx82 + 1.0*dwdx83;
    JB[1119] = -1.0*dwdx75;
    JB[1120] = -1.0*dwdx76;
    JB[1121] = -1.0*dwdx77;
    JB[1122] = -1.0*dwdx78;
    JB[1123] = -1.0*dwdx80;
    JB[1124] = -1.0*dwdx81;
    JB[1126] = -1.0*dwdx83;
    JB[1136] = 1.0*dwdx84;
    JB[1160] = 1.0*dwdx84;
    JB[1161] = -1.0*dwdx84;
    JB[1200] = 1.0*dwdx85;
    JB[1202] = 1.0*dwdx85;
    JB[1204] = -1.0*dwdx85;
    JB[1265] = 1.0*dwdx86;
    JB[1286] = 1.0*dwdx86;
    JB[1290] = -1.0*dwdx86;
    JB[1313] = 1.0*dwdx87;
    JB[1328] = 1.0*dwdx87;
    JB[1333] = -1.0*dwdx87;
    JB[1361] = 1.0*dwdx88;
    JB[1370] = 1.0*dwdx88;
    JB[1376] = -1.0*dwdx88;
    JB[1389] = 1.0*dwdx90;
    JB[1390] = -1.0*dwdx90;
    JB[1391] = -1.0*dwdx90;
    JB[1392] = 1.0*dwdx91;
    JB[1393] = 1.0*dwdx89;
    JB[1395] = -1.0*dwdx91;
    JB[1396] = -1.0*dwdx91;
    JB[1402] = 1.0*dwdx92;
    JB[1408] = 1.0*dwdx89;
    JB[1419] = -1.0*dwdx89 + 1.0*dwdx92;
    JB[1422] = -1.0*dwdx92;
    JB[1450] = 1.0*dwdx93;
    JB[1454] = 1.0*dwdx93;
    JB[1462] = -1.0*dwdx93;
    JB[1473] = 1.0*dwdx95;
    JB[1474] = -1.0*dwdx95;
    JB[1475] = -1.0*dwdx95;
    JB[1476] = 1.0*dwdx96;
    JB[1479] = -1.0*dwdx96;
    JB[1480] = -1.0*dwdx96;
    JB[1486] = 1.0*dwdx94;
    JB[1494] = 1.0*dwdx94;
    JB[1505] = -1.0*dwdx94;
    JB[1515] = 1.0*dwdx99;
    JB[1516] = -1.0*dwdx99;
    JB[1517] = -1.0*dwdx99;
    JB[1518] = 1.0*dwdx100;
    JB[1519] = 1.0*dwdx98;
    JB[1521] = -1.0*dwdx100;
    JB[1522] = -1.0*dwdx100;
    JB[1528] = 1.0*dwdx97;
    JB[1535] = 1.0*dwdx98;
    JB[1545] = 1.0*dwdx97;
    JB[1548] = -1.0*dwdx97 - 1.0*dwdx98;
    JB[1612] = 1.0*dwdx102;
    JB[1634] = -1.0*dwdx101 + 1.0*dwdx102;
    JB[1635] = 1.0*dwdx101;
    JB[1636] = -1.0*dwdx102;
    JB[1637] = 1.0*dwdx101;
    JB[1676] = -1.0*dwdx103;
    JB[1677] = 1.0*dwdx103;
    JB[1679] = 1.0*dwdx103;
    JB[1683] = 1.0*dwdx105;
    JB[1684] = -1.0*dwdx105;
    JB[1685] = -1.0*dwdx105;
    JB[1686] = 1.0*dwdx106;
    JB[1689] = -1.0*dwdx106;
    JB[1690] = -1.0*dwdx106;
    JB[1696] = 1.0*dwdx104;
    JB[1718] = 1.0*dwdx104;
    JB[1720] = -1.0*dwdx104;
    JB[1760] = -1.0*dwdx107;
    JB[1761] = 1.0*dwdx107;
    JB[1763] = 1.0*dwdx107;
}