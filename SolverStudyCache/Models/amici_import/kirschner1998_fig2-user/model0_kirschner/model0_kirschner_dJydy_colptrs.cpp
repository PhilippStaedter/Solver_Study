#include "sundials/sundials_types.h"

void dJydy_colptrs_model0_kirschner(sunindextype *colptrs, int index){
    switch(index) {
        case 0:
                colptrs[0] = 0;
                colptrs[1] = 1;
                colptrs[2] = 1;
            break;
        case 1:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 1;
            break;
    }
}