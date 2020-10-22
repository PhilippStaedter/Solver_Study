#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model1_hald(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(-EtOH*NAD*V2adh1/Kiadh1ACA - EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) - Kadh1NADH*V1adh1/Keq - Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) - NADH*V1adh1/Keq)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2) - 2*NADH*V1adh1*V2adh1*sc/(Keq*(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1));
    dwdx[1] = NAD*kOAc*sc;
    dwdx[2] = dACA*sc/Vratio;
    dwdx[3] = 0.5*kVap;
    dwdx[4] = HCN0*kacaf;
    dwdx[5] = -dACA*sc/Vratio;
    dwdx[6] = -2*ADP*kakr*sc;
    dwdx[7] = -ADP*PEP*Vpkm*sc/(pow(ADP + KpkADP, 2)*(KpkPEP + PEP)) + PEP*Vpkm*sc/((ADP + KpkADP)*(KpkPEP + PEP));
    dwdx[8] = DPG*klppepf*sc;
    dwdx[9] = ATP*kakf*sc;
    dwdx[10] = 2*pow(ATP, 3)*pow(F6P, 2)*Kpfk*Vpfkm*kappapfk*sc/(pow(AMP, 3)*(ATP + KpfkmATP)*pow(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2)), 2));
    dwdx[11] = AMP*kakf*sc;
    dwdx[12] = -0.83199999999999996*ATP*Vatpconm*sc/pow(ATP + KatpconATP, 2) + 0.83199999999999996*Vatpconm*sc/(ATP + KatpconATP);
    dwdx[13] = ATP*Glc*Vhkm*sc*(-Glc - KhkGlc)/pow(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc, 2) + Glc*Vhkm*sc/(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc);
    dwdx[14] = -ATP*pow(F6P, 2)*Vpfkm*sc/(pow(ATP + KpfkmATP, 2)*(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2)))) + pow(F6P, 2)*Vpfkm*sc/((ATP + KpfkmATP)*(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2)))) - 2*pow(ATP, 2)*pow(F6P, 2)*Kpfk*Vpfkm*kappapfk*sc/(pow(AMP, 2)*(ATP + KpfkmATP)*pow(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2)), 2));
    dwdx[15] = G6P*kg6pst*sc;
    dwdx[16] = -PEP*klppepr*sc;
    dwdx[17] = sc*(-DHAP*GAP*Valdf*(-GAP*Valdf/(Kaldeq*Valdr) - KaldGAP*Valdf/(Kaldeq*Valdr))/(Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) + FBP*Valdf*(-GAP*Valdf/(Kaldeq*Valdr) - KaldGAP*Valdf/(Kaldeq*Valdr))/pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2) - GAP*Valdf/(Kaldeq*(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP)));
    dwdx[18] = HCN*kdhapf;
    dwdx[19] = sc*(-DHAP*Vtimm/pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2) + GAP*Vtimm/(Ktimeq*pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2)) + Vtimm/(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP));
    dwdx[20] = DHAP*Vlpglycm*sc*(-KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH - 1)/pow(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1), 2) + Vlpglycm*sc/(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1));
    dwdx[21] = -kdhapb;
    dwdx[22] = sc*(DPG*NADH*Vgapdh/(KgapdhDPG*KgapdhGAP*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)) - GAP*NAD*Vgapdh/(KgapdhDPG*KgapdhGAP*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)) - NADH*Vgapdh/(KgapdhGAP*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
    dwdx[23] = ADP*klppepf*sc;
    dwdx[24] = 2*NAD*V1adh1*V2adh1*sc/(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1) + 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(-ACA*NAD*V2adh1/Kiadh1ACA - ACA*NADH*V1adh1/(Keq*Kiadh1EtOH) - Kadh1NAD*V2adh1 - Kadh1NAD*NADH*V2adh1/Kiadh1NADH - NAD*V2adh1)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
    dwdx[25] = dEtOH*sc/Vratio;
    dwdx[26] = -dEtOH*sc/Vratio;
    dwdx[27] = -2*ATP*pow(F6P, 3)*Vpfkm*sc/((ATP + KpfkmATP)*pow(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2)), 2)) + 2*ATP*F6P*Vpfkm*sc/((ATP + KpfkmATP)*(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2))));
    dwdx[28] = sc*(F6P*KpgiG6P*Vpgi/(KpgiF6P*Kpgieq*pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2)) - G6P*KpgiG6P*Vpgi/(KpgiF6P*pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2)) - Vpgi/(Kpgieq*(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P)));
    dwdx[29] = sc*(-DHAP*GAP*Valdf*(-GAP/KaldIGAP - 1)/(Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) + FBP*Valdf*(-GAP/KaldIGAP - 1)/pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2) + Valdf/(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP));
    dwdx[30] = sc*(-Glc*Vglctr*(-Glc/(KglctrGlc*KglctrIIG6P) - 1/KglctrIG6P)/(KglctrGlc*Vratio*pow(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1, 2)) - Glc0*Vglctr*(Glc/(KglctrGlc*KglctrIIG6P) + 1.0/KglctrIG6P)*(Glc0*Pglctr/KglctrGlc + 1)/(KglctrGlc*Vratio*(Glc*Pglctr/KglctrGlc + 1)*pow(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1), 2)));
    dwdx[31] = 0.69999999999999996*Pyr*V1PDC*sc*(-KaPDC/(KxPDC*(G6P*bet/(KxPDC*alp) + 1)) + KaPDC*bet*(G6P/KxPDC + 1)/(KxPDC*alp*pow(G6P*bet/(KxPDC*alp) + 1, 2)) + Pyr*bet*(G6P/(KxPDC*alp) + 1)/(KxPDC*alp*pow(G6P*bet/(KxPDC*alp) + 1, 2)) - Pyr/(KxPDC*alp*(G6P*bet/(KxPDC*alp) + 1)))/pow(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1), 2);
    dwdx[32] = sc*(F6P*Vpgi/(Kpgieq*pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2)) - G6P*Vpgi/pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2) + Vpgi/(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P));
    dwdx[33] = ATP*kg6pst*sc;
    dwdx[34] = sc*(-DHAP*GAP*Valdf*(-DHAP*Valdf/(Kaldeq*Valdr) - FBP/KaldIGAP - KaldDHAP*Valdf/(Kaldeq*Valdr))/(Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) - DHAP*Valdf/(Kaldeq*(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP)) + FBP*Valdf*(-DHAP*Valdf/(Kaldeq*Valdr) - FBP/KaldIGAP - KaldDHAP*Valdf/(Kaldeq*Valdr))/pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2));
    dwdx[35] = sc*(DPG*NADH*Vgapdh/(pow(KgapdhGAP, 2)*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)) - GAP*NAD*Vgapdh/(pow(KgapdhGAP, 2)*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)) + NAD*Vgapdh/(KgapdhGAP*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
    dwdx[36] = sc*(-DHAP*KtimDHAP*Vtimm/(KtimGAP*pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2)) + GAP*KtimDHAP*Vtimm/(KtimGAP*Ktimeq*pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2)) - Vtimm/(Ktimeq*(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP)));
    dwdx[37] = sc*(-Glc*Vglctr*(-G6P/(KglctrGlc*KglctrIIG6P) - Pglctr*(Glc0/KglctrGlc + 1)/(KglctrGlc*(Glc0*Pglctr/KglctrGlc + 1)) - 1/KglctrGlc)/(KglctrGlc*Vratio*pow(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1, 2)) + Glc0*Vglctr*(-(G6P/(KglctrGlc*KglctrIIG6P) + 1.0/KglctrGlc)*(Glc0*Pglctr/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1) + Pglctr*(Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(KglctrGlc*pow(Glc*Pglctr/KglctrGlc + 1, 2)))/(KglctrGlc*Vratio*pow(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1), 2)) - Vglctr/(KglctrGlc*Vratio*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1)));
    dwdx[38] = ATP*Glc*Vhkm*sc*(-ATP - KhkATP)/pow(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc, 2) + ATP*Vhkm*sc/(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc);
    dwdx[39] = sc*(-Glc*Vglctr*(Pglctr*(Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(KglctrGlc*pow(Glc0*Pglctr/KglctrGlc + 1, 2)) - (Glc*Pglctr/KglctrGlc + 1)/(KglctrGlc*(Glc0*Pglctr/KglctrGlc + 1)))/(KglctrGlc*Vratio*pow(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1, 2)) + Glc0*Vglctr*(-Pglctr*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(KglctrGlc*(Glc*Pglctr/KglctrGlc + 1)) - 1/KglctrGlc)/(KglctrGlc*Vratio*pow(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1), 2)) + Vglctr/(KglctrGlc*Vratio*(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1))));
    dwdx[40] = dGlyc*sc/Vratio;
    dwdx[41] = -dGlyc*sc/Vratio;
    dwdx[42] = DHAP*kdhapf;
    dwdx[43] = Pyr*kpyrf;
    dwdx[44] = dHCN*sc/Vratio;
    dwdx[45] = kVap;
    dwdx[46] = ACA0*kacaf;
    dwdx[47] = -dHCN*sc/Vratio;
    dwdx[48] = 2*EtOH*V1adh1*V2adh1*sc/(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1) + 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(-ACA*EtOH*V2adh1/Kiadh1ACA - ACA*Kadh1NADH*V1adh1/(Keq*Kiadh1NAD) - EtOH*V2adh1 - Kadh1EtOH*V2adh1)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
    dwdx[49] = ACA*kOAc*sc;
    dwdx[50] = sc*(DPG*NADH*Vgapdh/(KgapdhGAP*pow(KgapdhNAD, 2)*Kgapdheq*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + GAP*Vgapdh/(KgapdhGAP*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) - GAP*NAD*Vgapdh/(KgapdhGAP*pow(KgapdhNAD, 2)*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
    dwdx[51] = DHAP*Vlpglycm*sc*(-DHAP*KlpglycNADH/(KlpglycINAD*NADH) - KlpglycDHAP*KlpglycINADH/(KlpglycINAD*NADH))/pow(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1), 2);
    dwdx[52] = -2*ACA*V1adh1*V2adh1*sc/(Keq*(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1)) + 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(-ACA*EtOH*V1adh1/(Keq*Kiadh1EtOH) - ACA*V1adh1/Keq - EtOH*Kadh1NAD*V2adh1/Kiadh1NADH - Kadh1ACA*V1adh1/Keq)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
    dwdx[53] = sc*(-DPG*Vgapdh/(KgapdhGAP*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + DPG*NADH*Vgapdh/(KgapdhGAP*KgapdhNAD*KgapdhNADH*Kgapdheq*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) - GAP*NAD*Vgapdh/(KgapdhGAP*KgapdhNAD*KgapdhNADH*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
    dwdx[54] = DHAP*Vlpglycm*sc*(DHAP*KlpglycNADH*(1 + NAD/KlpglycINAD)/pow(NADH, 2) + KlpglycDHAP*KlpglycINADH*(1 + NAD/KlpglycINAD)/pow(NADH, 2))/pow(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1), 2);
    dwdx[55] = dOAc*sc/Vratio;
    dwdx[56] = -dOAc*sc/Vratio;
    dwdx[57] = -ADP*PEP*Vpkm*sc/((ADP + KpkADP)*pow(KpkPEP + PEP, 2)) + ADP*Vpkm*sc/((ADP + KpkADP)*(KpkPEP + PEP));
    dwdx[58] = -ATP*klppepr*sc;
    dwdx[59] = -0.69999999999999996*Pyr*V1PDC*sc*(G6P/(KxPDC*alp) + 1)/(pow(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1), 2)*(G6P*bet/(KxPDC*alp) + 1)) + 0.69999999999999996*V1PDC*sc/(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1));
    dwdx[60] = HCN*kpyrf;
    dwdx[61] = -1.1499999999999999*Pyr*VPYRd*sc/pow(KmPYRd + Pyr, 2) + 1.1499999999999999*VPYRd*sc/(KmPYRd + Pyr);
    dwdx[62] = -kpyrb;
    dwdx[63] = -kacab;
}