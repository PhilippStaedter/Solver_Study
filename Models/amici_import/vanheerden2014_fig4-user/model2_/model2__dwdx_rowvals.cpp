#include "sundials/sundials_types.h"

void dwdx_rowvals_model2_(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 1;
    rowvals[4] = 4;
    rowvals[5] = 1;
    rowvals[6] = 3;
}