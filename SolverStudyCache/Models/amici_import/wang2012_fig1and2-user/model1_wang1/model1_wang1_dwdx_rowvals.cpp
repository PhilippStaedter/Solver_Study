#include "sundials/sundials_types.h"

void dwdx_rowvals_model1_wang1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 3;
    rowvals[3] = 4;
    rowvals[4] = 5;
    rowvals[5] = 3;
    rowvals[6] = 6;
    rowvals[7] = 7;
    rowvals[8] = 4;
    rowvals[9] = 5;
    rowvals[10] = 8;
}