#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_montanez1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 4;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 4;
    rowvals[6] = 0;
    rowvals[7] = 1;
    rowvals[8] = 3;
    rowvals[9] = 4;
}