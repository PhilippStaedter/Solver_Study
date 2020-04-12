#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_bachmann(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = CIS;
    x_rdata[1] = CISRNA;
    x_rdata[2] = CISnRNA1;
    x_rdata[3] = CISnRNA2;
    x_rdata[4] = CISnRNA3;
    x_rdata[5] = CISnRNA4;
    x_rdata[6] = CISnRNA5;
    x_rdata[7] = Epo;
    x_rdata[8] = EpoRJAK2;
    x_rdata[9] = EpoRJAK2_CIS;
    x_rdata[10] = EpoRpJAK2;
    x_rdata[11] = SHP1;
    x_rdata[12] = SHP1Act;
    x_rdata[13] = SOCS3;
    x_rdata[14] = SOCS3RNA;
    x_rdata[15] = SOCS3nRNA1;
    x_rdata[16] = SOCS3nRNA2;
    x_rdata[17] = SOCS3nRNA3;
    x_rdata[18] = SOCS3nRNA4;
    x_rdata[19] = SOCS3nRNA5;
    x_rdata[20] = STAT5;
    x_rdata[21] = npSTAT5;
    x_rdata[22] = p12EpoRpJAK2;
    x_rdata[23] = p1EpoRpJAK2;
    x_rdata[24] = p2EpoRpJAK2;
    x_rdata[25] = pSTAT5;
}