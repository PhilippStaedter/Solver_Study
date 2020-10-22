#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model1_kouril4(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 8;
}