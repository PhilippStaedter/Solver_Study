#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_kouril7(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 0;
    rowvals[4] = 2;
    rowvals[5] = 3;
    rowvals[6] = 2;
    rowvals[7] = 3;
    rowvals[8] = 2;
    rowvals[9] = 0;
    rowvals[10] = 2;
    rowvals[11] = 1;
}