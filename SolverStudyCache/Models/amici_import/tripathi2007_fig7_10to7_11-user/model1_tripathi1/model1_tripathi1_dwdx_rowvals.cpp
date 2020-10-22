#include "sundials/sundials_types.h"

void dwdx_rowvals_model1_tripathi1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 9;
    rowvals[3] = 2;
    rowvals[4] = 4;
    rowvals[5] = 5;
    rowvals[6] = 6;
    rowvals[7] = 2;
    rowvals[8] = 7;
    rowvals[9] = 8;
    rowvals[10] = 2;
    rowvals[11] = 3;
}