#include "sundials/sundials_types.h"

void dwdx_rowvals_model1_valero(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 3;
}