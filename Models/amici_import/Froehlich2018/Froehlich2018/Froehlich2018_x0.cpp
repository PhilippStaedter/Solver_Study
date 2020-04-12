#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Froehlich2018(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 1.0;
    x0[3] = 1.0;
    x0[11] = 1.0;
    x0[21] = 1.0;
    x0[31] = 1.0;
    x0[41] = 1.0;
    x0[52] = 1.0;
    x0[62] = 1.0;
    x0[72] = 1.0;
    x0[82] = 1.0;
    x0[92] = 1.0;
    x0[102] = 1.0;
    x0[113] = 1.0;
    x0[119] = 1.0;
    x0[130] = 1.0;
    x0[141] = 1.0;
    x0[185] = 1.0;
    x0[212] = 1.0;
    x0[213] = 1.0;
    x0[214] = 1.0;
    x0[215] = 1.0;
    x0[216] = 1.0;
    x0[230] = 1.0;
    x0[238] = 1.0;
    x0[534] = 1.0;
    x0[560] = 1.0;
    x0[584] = 1.0;
    x0[608] = 1.0;
    x0[610] = 1.0;
    x0[622] = 1.0;
    x0[625] = 1.0;
    x0[628] = 1.0;
    x0[636] = 1.0;
    x0[651] = 1.0;
    x0[660] = 1.0;
    x0[661] = 1.0;
    x0[662] = 1.0;
    x0[675] = 1.0;
    x0[678] = 1.0;
    x0[699] = 1.0;
    x0[703] = 1.0;
    x0[724] = 1.0;
    x0[736] = 1.0;
    x0[764] = 1.0;
    x0[770] = 1.0;
    x0[780] = 1.0;
    x0[797] = 1.0;
    x0[831] = 1.0;
    x0[835] = 1.0;
    x0[846] = 1.0;
    x0[923] = 1.0;
    x0[934] = 1.0;
    x0[935] = 1.0;
    x0[940] = 1.0;
    x0[945] = 1.0;
    x0[947] = 1.0;
    x0[948] = 1.0;
    x0[968] = 1.0;
    x0[1002] = 1.0;
    x0[1013] = 1.0;
    x0[1036] = 1.0;
    x0[1044] = 1.0;
    x0[1048] = 1.0;
    x0[1060] = 1.0;
    x0[1066] = 1.0;
    x0[1084] = 1.0;
    x0[1089] = 1.0;
    x0[1109] = 1.0;
    x0[1114] = 1.0;
    x0[1134] = 1.0;
    x0[1139] = 1.0;
    x0[1150] = 1.0;
    x0[1158] = 1.0;
    x0[1170] = 1.0;
    x0[1194] = 1.0;
    x0[1195] = 1.0;
    x0[1236] = 1.0;
    x0[1238] = 1.0;
    x0[1245] = 1.0;
    x0[1278] = 1.0;
    x0[1280] = 1.0;
    x0[1282] = 1.0;
    x0[1287] = 1.0;
    x0[1288] = 1.0;
    x0[1299] = 1.0;
    x0[1304] = 1.0;
    x0[1307] = 1.0;
    x0[1313] = 1.0;
    x0[1319] = 1.0;
    x0[1332] = 1.0;
    x0[1337] = 1.0;
    x0[1338] = 1.0;
    x0[1339] = 1.0;
    x0[1340] = 1.0;
    x0[1341] = 1.0;
    x0[1342] = 1.0;
    x0[1343] = 1.0;
    x0[1344] = 1.0;
    x0[1345] = 1.0;
    x0[1347] = 1.0;
    x0[1348] = 1.0;
    x0[1349] = 1.0;
    x0[1350] = 1.0;
    x0[1351] = 1.0;
    x0[1352] = 1.0;
    x0[1353] = 1.0;
    x0[1355] = 1.0;
    x0[1356] = 1.0;
    x0[1357] = 1.0;
    x0[1358] = 1.0;
    x0[1359] = 1.0;
    x0[1361] = 1.0;
    x0[1362] = 1.0;
    x0[1370] = 1.0;
    x0[1371] = 1.0;
    x0[1372] = 1.0;
    x0[1373] = 1.0;
    x0[1375] = 1.0;
    x0[1380] = 1.0;
    x0[1384] = 1.0;
    x0[1385] = 1.0;
    x0[1386] = 1.0;
    x0[1394] = 1.0;
}