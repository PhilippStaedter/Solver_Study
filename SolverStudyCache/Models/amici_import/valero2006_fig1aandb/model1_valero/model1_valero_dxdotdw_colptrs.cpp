#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model1_valero(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 10;
}