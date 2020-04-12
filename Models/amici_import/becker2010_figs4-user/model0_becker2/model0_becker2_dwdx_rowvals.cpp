#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_becker2(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 4;
    rowvals[5] = 5;
    rowvals[6] = 6;
    rowvals[7] = 7;
}