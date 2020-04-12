#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_Morris2002_CellCycle_CDK2Cyclin(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 1;
}