#include "sundials/sundials_types.h"

void dJydy_colptrs_model0_montanez1(sunindextype *colptrs, int index){
    switch(index) {
        case 0:
                colptrs[0] = 0;
                colptrs[1] = 1;
                colptrs[2] = 1;
                colptrs[3] = 1;
            break;
        case 1:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 1;
                colptrs[3] = 1;
            break;
        case 2:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 1;
            break;
    }
}