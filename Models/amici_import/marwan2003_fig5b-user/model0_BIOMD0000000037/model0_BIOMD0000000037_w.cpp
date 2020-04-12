#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_BIOMD0000000037(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Gluc*Glucose_sensor_inactivation_kG*Ya;
    w[1] = 1.0*Pfr*Photoreceptor_activation_IfrSfrPfr;
    w[2] = 1.0*Photoreceptor_decay_kd*Pr;
    w[3] = 1.0*Photoreceptor_inactivation_IrSrPr*Pr;
    w[4] = 1.0*S*S_degradation_kd_s;
    w[5] = 1.0*S_formation_alpha1/(pow(V, 3) + 1);
    w[6] = 1.0*S_generation_ky*Ya*preS;
    w[7] = 1.0*Pr*Transducer_activation_kia*Xi;
    w[8] = 1.0*Transducer_inactivation_kai*Xa;
    w[9] = 1.0*V*V_degradation_kd_v;
    w[10] = 1.0*V_formation_alpha2/(pow(S, 3) + 1);
    w[11] = 1.0*Xa*preS_formation_kx*prepreS;
}