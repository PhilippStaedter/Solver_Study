#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_becker1(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 4;
    rowvals[5] = 5;
    rowvals[6] = 6;
    rowvals[7] = 7;
}