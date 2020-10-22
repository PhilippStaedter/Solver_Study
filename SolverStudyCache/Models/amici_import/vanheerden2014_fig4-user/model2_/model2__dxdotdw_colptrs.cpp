#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model2_(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 7;
    colptrs[4] = 8;
    colptrs[5] = 10;
}