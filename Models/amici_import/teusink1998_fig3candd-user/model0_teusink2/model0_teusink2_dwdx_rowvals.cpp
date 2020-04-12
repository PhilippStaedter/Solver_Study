#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_teusink2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 2;
    rowvals[5] = 0;
    rowvals[6] = 0;
    rowvals[7] = 1;
}