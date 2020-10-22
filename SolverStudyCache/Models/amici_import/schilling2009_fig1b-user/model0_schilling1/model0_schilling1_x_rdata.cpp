#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_schilling1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Delay01_mSHP1;
    x_rdata[1] = Delay02_mSHP1;
    x_rdata[2] = Delay03_mSHP1;
    x_rdata[3] = Delay04_mSHP1;
    x_rdata[4] = Delay05_mSHP1;
    x_rdata[5] = Delay06_mSHP1;
    x_rdata[6] = Delay07_mSHP1;
    x_rdata[7] = Delay08_mSHP1;
    x_rdata[8] = ERK1;
    x_rdata[9] = ERK2;
    x_rdata[10] = Epo;
    x_rdata[11] = EpoR;
    x_rdata[12] = JAK2;
    x_rdata[13] = MEK1;
    x_rdata[14] = MEK2;
    x_rdata[15] = Raf;
    x_rdata[16] = SHP1;
    x_rdata[17] = SOS;
    x_rdata[18] = actSHP1;
    x_rdata[19] = mSHP1;
    x_rdata[20] = mSOS;
    x_rdata[21] = pERK1;
    x_rdata[22] = pERK2;
    x_rdata[23] = pEpoR;
    x_rdata[24] = pJAK2;
    x_rdata[25] = pMEK1;
    x_rdata[26] = pMEK2;
    x_rdata[27] = pRaf;
    x_rdata[28] = pSOS;
    x_rdata[29] = ppERK1;
    x_rdata[30] = ppERK2;
    x_rdata[31] = ppMEK1;
    x_rdata[32] = ppMEK2;
}