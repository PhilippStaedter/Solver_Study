#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_Morris2002_CellCycle_CDK2Cyclin(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
}