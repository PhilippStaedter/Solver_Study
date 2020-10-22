#include "sundials/sundials_types.h"

void JSparse_colptrs_model1_alexander2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 3;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 10;
}