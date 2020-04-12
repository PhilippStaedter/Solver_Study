#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_bachmann(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = CIS;
    x_solver[1] = CISRNA;
    x_solver[2] = CISnRNA1;
    x_solver[3] = CISnRNA2;
    x_solver[4] = CISnRNA3;
    x_solver[5] = CISnRNA4;
    x_solver[6] = CISnRNA5;
    x_solver[7] = Epo;
    x_solver[8] = EpoRJAK2;
    x_solver[9] = EpoRJAK2_CIS;
    x_solver[10] = EpoRpJAK2;
    x_solver[11] = SHP1;
    x_solver[12] = SHP1Act;
    x_solver[13] = SOCS3;
    x_solver[14] = SOCS3RNA;
    x_solver[15] = SOCS3nRNA1;
    x_solver[16] = SOCS3nRNA2;
    x_solver[17] = SOCS3nRNA3;
    x_solver[18] = SOCS3nRNA4;
    x_solver[19] = SOCS3nRNA5;
    x_solver[20] = STAT5;
    x_solver[21] = npSTAT5;
    x_solver[22] = p12EpoRpJAK2;
    x_solver[23] = p1EpoRpJAK2;
    x_solver[24] = p2EpoRpJAK2;
    x_solver[25] = pSTAT5;
}