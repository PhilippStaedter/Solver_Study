#include "sundials/sundials_types.h"

void dwdx_rowvals_model24_karin3(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 1;
}