#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_Sivakumar2011c(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 + 1.0*dwdx1;
    J[1] = -1.0*dwdx2;
    J[2] = -1.0*dwdx4;
    J[24] = 1.0*dwdx32;
    J[50] = 1.0*dwdx0;
    J[51] = 1.0*dwdx2 - 1.0*dwdx3;
    J[52] = 1.0*dwdx4;
    J[53] = -1.0*dwdx5;
    J[54] = -1.0*dwdx6;
    J[100] = -1.0*dwdx0;
    J[101] = -1.0*dwdx2;
    J[102] = -1.0*dwdx4;
    J[151] = 1.0*dwdx3;
    J[153] = 1.0*dwdx5;
    J[154] = 1.0*dwdx6;
    J[201] = -1.0*dwdx3;
    J[203] = -1.0*dwdx5;
    J[204] = -1.0*dwdx6;
    J[255] = 1.0*dwdx7;
    J[291] = 1.0*dwdx67;
    J[306] = -1.0*dwdx10;
    J[341] = -1.0*dwdx66;
    J[342] = -1.0*dwdx69;
    J[355] = -1.0*dwdx8 - 1.0*dwdx9;
    J[357] = -1.0*dwdx11 - 1.0*dwdx12 - 3.0*dwdx13;
    J[358] = -1.0*dwdx14 - 1.0*dwdx15 - 3.0*dwdx16;
    J[359] = -3.0*dwdx18;
    J[364] = -1.0*dwdx23;
    J[377] = -1.0*dwdx39;
    J[384] = -1.0*dwdx52 - 3.0*dwdx53;
    J[385] = -3.0*dwdx55;
    J[391] = -1.0*dwdx68;
    J[398] = -3.0*dwdx81;
    J[405] = 1.0*dwdx8 + 1.0*dwdx9;
    J[407] = 1.0*dwdx11 + 1.0*dwdx12 + 3.0*dwdx13;
    J[408] = 1.0*dwdx14 + 1.0*dwdx15 + 3.0*dwdx16;
    J[409] = 3.0*dwdx18;
    J[414] = 1.0*dwdx23;
    J[427] = 1.0*dwdx39;
    J[434] = 1.0*dwdx52 + 3.0*dwdx53;
    J[435] = 3.0*dwdx55;
    J[441] = 1.0*dwdx68;
    J[448] = 3.0*dwdx81;
    J[459] = 1.0*dwdx17;
    J[493] = 1.0*dwdx72;
    J[510] = -1.0*dwdx19;
    J[526] = -1.0*dwdx36;
    J[527] = -1.0*dwdx38;
    J[561] = -1.0*dwdx20;
    J[576] = -1.0*dwdx37;
    J[578] = -1.0*dwdx40;
    J[612] = 1.0*dwdx21;
    J[639] = 1.0*dwdx63;
    J[640] = 1.0*dwdx64;
    J[663] = -1.0*dwdx22;
    J[692] = -1.0*dwdx70;
    J[693] = -1.0*dwdx71;
    J[705] = -1.0*dwdx8;
    J[707] = -1.0*dwdx11;
    J[708] = -1.0*dwdx14;
    J[714] = -1.0*dwdx23;
    J[741] = -1.0*dwdx68;
    J[765] = -1.0*dwdx24;
    J[775] = -1.0*dwdx33;
    J[776] = -1.0*dwdx35;
    J[816] = 1.0*dwdx25;
    J[846] = 1.0*dwdx78;
    J[847] = 1.0*dwdx79;
    J[867] = 1.0*dwdx26;
    J[888] = 1.0*dwdx60;
    J[890] = 1.0*dwdx65;
    J[918] = -1.0*dwdx27;
    J[931] = -1.0*dwdx46;
    J[932] = -1.0*dwdx49;
    J[969] = -1.0*dwdx28;
    J[983] = -1.0*dwdx51;
    J[999] = -1.0*dwdx83;
    J[1020] = -1.0*dwdx29;
    J[1032] = -1.0*dwdx48;
    J[1033] = -1.0*dwdx50;
    J[1071] = -1.0*dwdx30;
    J[1080] = -1.0*dwdx45;
    J[1081] = -1.0*dwdx47;
    J[1149] = 1.0*dwdx82;
    J[1173] = -1.0*dwdx31;
    J[1179] = -1.0*dwdx43;
    J[1180] = -1.0*dwdx44;
    J[1200] = -1.0*dwdx1;
    J[1224] = -1.0*dwdx32;
    J[1265] = -1.0*dwdx24;
    J[1275] = -1.0*dwdx33 - 1.0*dwdx34;
    J[1276] = -1.0*dwdx35;
    J[1279] = -1.0*dwdx42;
    J[1310] = -1.0*dwdx19;
    J[1311] = -1.0*dwdx20;
    J[1315] = 1.0*dwdx24;
    J[1325] = 1.0*dwdx33;
    J[1326] = 1.0*dwdx35 - 1.0*dwdx36 - 1.0*dwdx37;
    J[1327] = -1.0*dwdx38;
    J[1328] = -1.0*dwdx40;
    J[1355] = -1.0*dwdx9;
    J[1357] = -1.0*dwdx12;
    J[1358] = -1.0*dwdx15;
    J[1360] = 1.0*dwdx19;
    J[1376] = 1.0*dwdx36;
    J[1377] = 1.0*dwdx38 - 1.0*dwdx39;
    J[1384] = -1.0*dwdx52;
    J[1411] = 1.0*dwdx20;
    J[1426] = 1.0*dwdx37;
    J[1428] = 1.0*dwdx40 - 1.0*dwdx41;
    J[1436] = -1.0*dwdx57;
    J[1439] = -1.0*dwdx62;
    J[1473] = -1.0*dwdx31;
    J[1475] = 1.0*dwdx34;
    J[1479] = 1.0*dwdx42 - 1.0*dwdx43;
    J[1480] = -1.0*dwdx44;
    J[1487] = 1.0*dwdx59;
    J[1488] = 1.0*dwdx61;
    J[1521] = -1.0*dwdx30;
    J[1523] = 1.0*dwdx31;
    J[1529] = 1.0*dwdx43;
    J[1530] = 1.0*dwdx44 - 1.0*dwdx45;
    J[1531] = -1.0*dwdx47;
    J[1568] = -1.0*dwdx27;
    J[1571] = 1.0*dwdx30;
    J[1580] = 1.0*dwdx45;
    J[1581] = -1.0*dwdx46 + 1.0*dwdx47;
    J[1582] = -1.0*dwdx49;
    J[1618] = 1.0*dwdx27;
    J[1620] = -1.0*dwdx29;
    J[1631] = 1.0*dwdx46;
    J[1632] = -1.0*dwdx48 + 1.0*dwdx49;
    J[1633] = -1.0*dwdx50;
    J[1669] = -1.0*dwdx28;
    J[1670] = 1.0*dwdx29;
    J[1682] = 1.0*dwdx48;
    J[1683] = 1.0*dwdx50 - 1.0*dwdx51;
    J[1699] = -1.0*dwdx83;
    J[1705] = 1.0*dwdx9;
    J[1707] = 1.0*dwdx12 - 1.0*dwdx13;
    J[1708] = 1.0*dwdx15 - 1.0*dwdx16;
    J[1709] = -1.0*dwdx18;
    J[1727] = 1.0*dwdx39;
    J[1734] = 1.0*dwdx52 - 1.0*dwdx53;
    J[1735] = -1.0*dwdx55;
    J[1748] = -1.0*dwdx81;
    J[1757] = 1.0*dwdx13;
    J[1758] = 1.0*dwdx16;
    J[1759] = 1.0*dwdx18;
    J[1784] = 1.0*dwdx53;
    J[1785] = -1.0*dwdx54 + 1.0*dwdx55;
    J[1794] = -1.0*dwdx73;
    J[1795] = -1.0*dwdx75;
    J[1798] = 1.0*dwdx81;
    J[1828] = -1.0*dwdx41;
    J[1836] = -1.0*dwdx56 - 1.0*dwdx57;
    J[1839] = -1.0*dwdx62;
    J[1845] = -1.0*dwdx76;
    J[1846] = -1.0*dwdx77;
    J[1887] = 1.0*dwdx58 - 1.0*dwdx59;
    J[1894] = 1.0*dwdx74;
    J[1897] = 1.0*dwdx80;
    J[1917] = 1.0*dwdx26;
    J[1938] = 1.0*dwdx60 - 1.0*dwdx61;
    J[1940] = 1.0*dwdx65;
    J[1962] = -1.0*dwdx21;
    J[1978] = 1.0*dwdx41;
    J[1986] = 1.0*dwdx57;
    J[1989] = 1.0*dwdx62 - 1.0*dwdx63;
    J[1990] = -1.0*dwdx64;
    J[2012] = 1.0*dwdx21;
    J[2017] = -1.0*dwdx26;
    J[2038] = -1.0*dwdx60;
    J[2039] = 1.0*dwdx63;
    J[2040] = 1.0*dwdx64 - 1.0*dwdx65;
    J[2055] = -1.0*dwdx7 + 1.0*dwdx8;
    J[2056] = -1.0*dwdx10;
    J[2057] = 1.0*dwdx11;
    J[2058] = 1.0*dwdx14;
    J[2064] = 1.0*dwdx23;
    J[2091] = -1.0*dwdx66 - 1.0*dwdx67 + 1.0*dwdx68;
    J[2092] = -1.0*dwdx69;
    J[2106] = 1.0*dwdx10;
    J[2113] = -1.0*dwdx22;
    J[2141] = 1.0*dwdx66;
    J[2142] = 1.0*dwdx69 - 1.0*dwdx70;
    J[2143] = -1.0*dwdx71;
    J[2159] = -1.0*dwdx17;
    J[2163] = 1.0*dwdx22;
    J[2192] = 1.0*dwdx70;
    J[2193] = 1.0*dwdx71 - 1.0*dwdx72;
    J[2235] = -1.0*dwdx54;
    J[2237] = 1.0*dwdx58;
    J[2244] = -1.0*dwdx73 + 1.0*dwdx74;
    J[2245] = -1.0*dwdx75;
    J[2247] = 1.0*dwdx80;
    J[2285] = 1.0*dwdx54;
    J[2286] = -1.0*dwdx56;
    J[2294] = 1.0*dwdx73;
    J[2295] = 1.0*dwdx75 - 1.0*dwdx76;
    J[2296] = -1.0*dwdx77;
    J[2316] = -1.0*dwdx25;
    J[2336] = 1.0*dwdx56;
    J[2345] = 1.0*dwdx76;
    J[2346] = 1.0*dwdx77 - 1.0*dwdx78;
    J[2347] = -1.0*dwdx79;
    J[2366] = 1.0*dwdx25;
    J[2387] = -1.0*dwdx58;
    J[2394] = -1.0*dwdx74;
    J[2396] = 1.0*dwdx78;
    J[2397] = 1.0*dwdx79 - 1.0*dwdx80;
    J[2469] = 1.0*dwdx28;
    J[2483] = 1.0*dwdx51;
    J[2499] = -1.0*dwdx82 + 1.0*dwdx83;
}