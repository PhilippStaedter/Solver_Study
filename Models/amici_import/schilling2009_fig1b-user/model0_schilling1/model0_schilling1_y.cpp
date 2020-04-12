#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_schilling1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Delay01_mSHP1;
    y[1] = Delay02_mSHP1;
    y[2] = Delay03_mSHP1;
    y[3] = Delay04_mSHP1;
    y[4] = Delay05_mSHP1;
    y[5] = Delay06_mSHP1;
    y[6] = Delay07_mSHP1;
    y[7] = Delay08_mSHP1;
    y[8] = ERK1;
    y[9] = ERK2;
    y[10] = Epo;
    y[11] = EpoR;
    y[12] = JAK2;
    y[13] = MEK1;
    y[14] = MEK2;
    y[15] = Raf;
    y[16] = SHP1;
    y[17] = SOS;
    y[18] = actSHP1;
    y[19] = mSHP1;
    y[20] = mSOS;
    y[21] = pERK1;
    y[22] = pERK2;
    y[23] = pEpoR;
    y[24] = pJAK2;
    y[25] = pMEK1;
    y[26] = pMEK2;
    y[27] = pRaf;
    y[28] = pSOS;
    y[29] = ppERK1;
    y[30] = ppERK2;
    y[31] = ppMEK1;
    y[32] = ppMEK2;
}