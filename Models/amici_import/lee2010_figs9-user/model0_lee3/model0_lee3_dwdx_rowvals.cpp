#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_lee3(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 2;
    rowvals[2] = 1;
    rowvals[3] = 3;
}