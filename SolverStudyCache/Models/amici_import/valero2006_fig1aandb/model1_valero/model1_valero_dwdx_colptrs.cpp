#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_valero(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 4;
    colptrs[4] = 4;
    colptrs[5] = 4;
    colptrs[6] = 5;
}