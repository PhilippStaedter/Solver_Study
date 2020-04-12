#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_perelson1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 4;
    rowvals[2] = 0;
    rowvals[3] = 2;
    rowvals[4] = 3;
}