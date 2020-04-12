#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_kirschner(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 3;
    rowvals[2] = 5;
    rowvals[3] = 1;
    rowvals[4] = 3;
    rowvals[5] = 4;
    rowvals[6] = 5;
}