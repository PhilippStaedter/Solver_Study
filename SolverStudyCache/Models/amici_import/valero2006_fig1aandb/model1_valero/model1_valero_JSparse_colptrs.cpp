#include "sundials/sundials_types.h"

void JSparse_colptrs_model1_valero(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 9;
    colptrs[4] = 9;
    colptrs[5] = 9;
    colptrs[6] = 11;
}