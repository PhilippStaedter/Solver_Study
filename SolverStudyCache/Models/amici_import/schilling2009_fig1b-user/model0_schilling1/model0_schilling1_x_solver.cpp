#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_schilling1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Delay01_mSHP1;
    x_solver[1] = Delay02_mSHP1;
    x_solver[2] = Delay03_mSHP1;
    x_solver[3] = Delay04_mSHP1;
    x_solver[4] = Delay05_mSHP1;
    x_solver[5] = Delay06_mSHP1;
    x_solver[6] = Delay07_mSHP1;
    x_solver[7] = Delay08_mSHP1;
    x_solver[8] = ERK1;
    x_solver[9] = ERK2;
    x_solver[10] = Epo;
    x_solver[11] = EpoR;
    x_solver[12] = JAK2;
    x_solver[13] = MEK1;
    x_solver[14] = MEK2;
    x_solver[15] = Raf;
    x_solver[16] = SHP1;
    x_solver[17] = SOS;
    x_solver[18] = actSHP1;
    x_solver[19] = mSHP1;
    x_solver[20] = mSOS;
    x_solver[21] = pERK1;
    x_solver[22] = pERK2;
    x_solver[23] = pEpoR;
    x_solver[24] = pJAK2;
    x_solver[25] = pMEK1;
    x_solver[26] = pMEK2;
    x_solver[27] = pRaf;
    x_solver[28] = pSOS;
    x_solver[29] = ppERK1;
    x_solver[30] = ppERK2;
    x_solver[31] = ppMEK1;
    x_solver[32] = ppMEK2;
}