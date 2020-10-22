#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model8_panteleev3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 11;
    colptrs[5] = 14;
}