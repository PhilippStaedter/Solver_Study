#include "sundials/sundials_types.h"

void dxdotdw_colptrs_dupont2_Fig4B(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 3;
    colptrs[4] = 5;
    colptrs[5] = 7;
    colptrs[6] = 8;
    colptrs[7] = 10;
}