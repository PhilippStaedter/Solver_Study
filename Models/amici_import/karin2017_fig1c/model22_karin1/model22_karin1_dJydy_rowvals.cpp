#include "sundials/sundials_types.h"

void dJydy_rowvals_model22_karin1(sunindextype *rowvals, int index){
    switch(index) {
        case 0:
                rowvals[0] = 0;
            break;
        case 1:
                rowvals[0] = 0;
            break;
    }
}