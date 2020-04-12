#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_ouattara1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 3;
    rowvals[3] = 4;
    rowvals[4] = 2;
    rowvals[5] = 5;
}