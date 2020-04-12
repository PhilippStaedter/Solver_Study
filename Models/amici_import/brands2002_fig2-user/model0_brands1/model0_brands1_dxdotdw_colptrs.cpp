#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_brands1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 7;
    colptrs[4] = 9;
    colptrs[5] = 12;
    colptrs[6] = 15;
    colptrs[7] = 17;
    colptrs[8] = 20;
    colptrs[9] = 23;
    colptrs[10] = 26;
    colptrs[11] = 28;
}