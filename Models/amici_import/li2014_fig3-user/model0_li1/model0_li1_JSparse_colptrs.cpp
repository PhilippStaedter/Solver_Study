#include "sundials/sundials_types.h"

void JSparse_colptrs_model0_li1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 8;
}