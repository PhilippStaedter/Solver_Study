#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_bier2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 3;
    rowvals[4] = 5;
    rowvals[5] = 1;
    rowvals[6] = 4;
    rowvals[7] = 5;
}