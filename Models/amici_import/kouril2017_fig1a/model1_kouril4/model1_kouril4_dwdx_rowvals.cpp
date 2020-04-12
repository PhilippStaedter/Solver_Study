#include "sundials/sundials_types.h"

void dwdx_rowvals_model1_kouril4(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 0;
    rowvals[4] = 0;
    rowvals[5] = 1;
}