#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>

#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"

#include "sundials/sundials_types.h"

namespace amici {
class Solver;
}

/**
 * @brief Wrapper function to instantiate the linked Amici model without knowing
 * the name at compile time.
 * @return
 */
extern void J_Nakakuki2010(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Nakakuki2010(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Nakakuki2010(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Nakakuki2010(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Nakakuki2010(sunindextype *colptrs);
extern void JSparse_rowvals_Nakakuki2010(sunindextype *rowvals);
extern void JSparseB_Nakakuki2010(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Nakakuki2010(sunindextype *colptrs);
extern void JSparseB_rowvals_Nakakuki2010(sunindextype *rowvals);
extern void Jy_Nakakuki2010(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Nakakuki2010(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Nakakuki2010(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Nakakuki2010(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Nakakuki2010(sunindextype *rowvals, int index);
extern void dwdp_Nakakuki2010(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Nakakuki2010(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Nakakuki2010(sunindextype *colptrs);
extern void dwdx_rowvals_Nakakuki2010(sunindextype *rowvals);
extern void dxdotdw_Nakakuki2010(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Nakakuki2010(sunindextype *colptrs);
extern void dxdotdw_rowvals_Nakakuki2010(sunindextype *rowvals);
extern void dxdotdp_Nakakuki2010(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Nakakuki2010(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Nakakuki2010(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Nakakuki2010(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Nakakuki2010(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Nakakuki2010(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Nakakuki2010(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Nakakuki2010(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Nakakuki2010(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Nakakuki2010(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Nakakuki2010(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Nakakuki2010(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Nakakuki2010(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Nakakuki2010(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Nakakuki2010 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Nakakuki2010()
        : amici::Model_ODE(
              49,                                // nx_rdata
              49,                            // nxtrue_rdata
              49,                               // nx_solver
              49,                           // nxtrue_solver
              49,                                      // ny
              49,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              78,                                      // nw
              129,                                   // ndwdx
              10998,                                   // ndwdp
              147,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              193,                                     // nnz
              49,                                     // ubw
              49,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{0.22, 0.72, 160.0, 0.648, 60.0, 19.49872346, 29.94073716, 19.4987234631759, 29.9407371620698, 0.012, 0.018, 0.012, 0.018, 0.011, 0.013, 29.24109258, 169.0473748, 3.970849295, 0.000126129, 0.007875765, 0.001245747, 5.636949216, 34180.48, 2.992346912, 0.001172165, 0.000257, 4.81e-05, 0.024269764, 0.070467899, 0.024269764, 0.070467899, 0.157678678, 735598.6967, 0.005648117, 387.8377182, 0.000257, 4.81e-05, 0.550346114, 29516.06587, 0.342848369838443, 307.041525298866, 10.09063736, 0.913939859, 0.025925065, 0.129803956, 19.23118154, 441.5834425, 6.574759504, 14.99180922, 0.518529841, 21312.69109, 13.79479021, 15.04396629, 0.655214248, 185.9760682, 1.988003164, 0.003284434, 0.000601234209304622, 7.64816282169636e-05, 8.907637012, 8562.744184, 0.000597315, 528.552427, 1.745848179, 0.070379236, 0.000257, 4.81e-05, 0.54528521, 0.133249762, 0.54528521, 0.133249762, 0.909968714, 3992.061328, 0.076717457, 1157.116021, 0.078344305, 0.051168202, 0.000257, 4.81e-05, 0.001670815, 15.80783969, 0.686020478, 0.314470502, 2.335459127, 26.59483436, 0.01646825, 9.544308421, 0.7485, 0.001670815, 15.80783969, 0.686020478, 0.314470502, 2.335459127, 26.59483436, 0.01646825, 9.544308421, 1.026834758, 0.637490056, 3.584464176, 0.000270488, 0.001443889, 0.002448164, 3.49860901414122e-05, 0.019898797, 0.396950616, 4.13466150826031e-05, 0.013844393, 2.800340453, 350.0, 220.0, 940.0, 0.01807448, 3475.168, 0.09858154, 237.2001, 0.3573399, 1334.132, 4.635749, 4046.71, 0.05393704, 1.027895, 0.109304, 606.871, 5.291093, 424.6884, 0.03436149, 11.5048, 0.1374307, 7424.816, 0.08258693, 425.5268, 0.02487469, 858.3423, 0.8850982, 4665.217, 0.05377297, 20.50809, 0.03957055, 7.774197, 13.74244, 2122.045},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(49, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Nakakuki2010(*this);
    }

    /** model specific implementation for fJ
     * @param J Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJ(realtype *J, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w, const realtype *dwdx) override {
        J_Nakakuki2010(J, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation for fJB
     * @param JB Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJB(realtype *JB, const realtype t, const realtype *x,
                     const realtype *p, const realtype *k, const realtype *h,
                     const realtype *xB, const realtype *w,
                     const realtype *dwdx) override {
        JB_Nakakuki2010(JB, t, x, p, k, h, xB, w, dwdx);
    }

    /** model specific implementation for fJDiag
     * @param JDiag Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJDiag(realtype *JDiag, const realtype t, const realtype *x,
                        const realtype *p, const realtype *k, const realtype *h,
                        const realtype *w, const realtype *dwdx) override {
        JDiag_Nakakuki2010(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Nakakuki2010( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Nakakuki2010(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Nakakuki2010(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Nakakuki2010( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Nakakuki2010(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Nakakuki2010(rowvals);
    }


    /** model specific implementation of fJrz
     * @param nllh regularization for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fJrz(realtype *nllh, const int iz, const realtype *p,
                      const realtype *k, const realtype *rz,
                      const realtype *sigmaz) override {}

    /** model specific implementation of fJy
     * @param nllh negative log-likelihood for measurements y
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurements at timepoint
     **/
    virtual void fJy(realtype *nllh, const int iy, const realtype *p,
                     const realtype *k, const realtype *y,
                     const realtype *sigmay, const realtype *my) override {
        Jy_Nakakuki2010(nllh, iy, p, k, y, sigmay, my);
    }

    /** model specific implementation of fJz
     * @param nllh negative log-likelihood for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurements at timepoint
     **/
    virtual void fJz(realtype *nllh, const int iz, const realtype *p,
                     const realtype *k, const realtype *z,
                     const realtype *sigmaz, const realtype *mz) override {}

    /** model specific implementation of fdJrzdsigma
     * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
     * standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdsigma(realtype *dJrzdsigma, const int iz,
                             const realtype *p, const realtype *k,
                             const realtype *rz,
                             const realtype *sigmaz) override {}

    /** model specific implementation of fdJrzdz
     * @param dJrzdz partial derivative of event penalization Jrz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p,
                         const realtype *k, const realtype *rz,
                         const realtype *sigmaz) override {}

    /** model specific implementation of fdJydsigma
     * @param dJydsigma Sensitivity of time-resolved measurement
     * negative log-likelihood Jy w.r.t. standard deviation sigmay
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurement at timepoint
     **/
    virtual void fdJydsigma(realtype *dJydsigma, const int iy,
                            const realtype *p, const realtype *k,
                            const realtype *y, const realtype *sigmay,
                            const realtype *my) override {
        dJydsigmay_Nakakuki2010(dJydsigma, iy, p, k, y, sigmay, my);
    }


    /** model specific implementation of fdJzdsigma
     * @param dJzdsigma Sensitivity of event measurement
     * negative log-likelihood Jz w.r.t. standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdsigma(realtype *dJzdsigma, const int iz,
                            const realtype *p, const realtype *k,
                            const realtype *z, const realtype *sigmaz,
                            const realtype *mz) override {}

    /** model specific implementation of fdJzdz
     * @param dJzdz partial derivative of event measurement negative
     *log-likelihood Jz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdz(realtype *dJzdz, const int iz, const realtype *p,
                        const realtype *k, const realtype *z,
                        const realtype *sigmaz, const realtype *mz) override {}

    /** model specific implementation of fdeltasx
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     **/
    virtual void fdeltaqB(realtype *deltaqB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ip,
                          const int ie, const realtype *xdot,
                          const realtype *xdot_old,
                          const realtype *xB) override {}

    /** model specific implementation of fdeltasx
     * @param deltasx sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param w repeating elements vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param sx state sensitivity
     * @param stau event-time sensitivity
     **/
    virtual void fdeltasx(realtype *deltasx, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h,
                          const realtype *w, const int ip, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *sx, const realtype *stau) override {}

    /** model specific implementation of fdeltax
     * @param deltax state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     **/
    virtual void fdeltax(realtype *deltax, const realtype t, const realtype *x,
                         const realtype *p, const realtype *k,
                         const realtype *h, const int ie, const realtype *xdot,
                         const realtype *xdot_old) override {}

    /** model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     **/
    virtual void fdeltaxB(realtype *deltaxB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *xB) override {}

    /** model specific implementation of fdrzdp
     * @param drzdp partial derivative of root output rz w.r.t. model parameters
     *p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdrzdp(realtype *drzdp, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h, const int ip) override {}

    /** model specific implementation of fdrzdx
     * @param drzdx partial derivative of root output rz w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdrzdx(realtype *drzdx, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h) override {}

    /** model specific implementation of fsigmay
     * @param dsigmaydp partial derivative of standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmaydp(realtype *dsigmaydp, const realtype t,
                            const realtype *p, const realtype *k,
                            const int ip) override {
        dsigmaydp_Nakakuki2010(dsigmaydp, t, p, k, ip);
    }

    /** model specific implementation of fsigmaz
     * @param dsigmazdp partial derivative of standard deviation of event
     *measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmazdp(realtype *dsigmazdp, const realtype t,
                            const realtype *p, const realtype *k,
                            const int ip) override {}

    virtual void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_Nakakuki2010( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Nakakuki2010(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Nakakuki2010(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Nakakuki2010( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Nakakuki2010( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Nakakuki2010( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Nakakuki2010(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Nakakuki2010(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Nakakuki2010( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
    }


    /** model specific implementation of fdydx
     * @param dydx partial derivative of observables y w.r.t. model states x
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdydx(realtype *dydx, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w, const realtype *dwdx) override {
        dydx_Nakakuki2010(dydx, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation of fdydp
     * @param dydp partial derivative of observables y w.r.t. model parameters p
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdydp(realtype *dydp, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const int ip, const realtype *w,
                       const realtype *dwdp) override {
        dydp_Nakakuki2010(dydp, t, x, p, k, h, ip, w, dwdp);
    }

    /** model specific implementation of fdzdp
     * @param dzdp partial derivative of event-resolved output z w.r.t. model
     *parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdzdp(realtype *dzdp, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h, const int ip) override {}

    /** model specific implementation of fdzdx
     * @param dzdx partial derivative of event-resolved output z w.r.t. model
     *states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdzdx(realtype *dzdx, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h) override {}

    /** model specific implementation for froot
     * @param root values of the trigger function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     **/
    virtual void froot(realtype *root, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k,
                       const realtype *h) override {}

    /** model specific implementation of frz
     * @param rz value of root function at current timepoint (non-output events
     *not included)
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void frz(realtype *rz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h) override {}

    /** model specific implementation of fsigmay
     * @param sigmay standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmay(realtype *sigmay, const realtype t, const realtype *p,
                         const realtype *k) override {
        sigmay_Nakakuki2010(sigmay, t, p, k);
    }

    /** model specific implementation of fsigmaz
     * @param sigmaz standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p,
                         const realtype *k) override {}

    /** model specific implementation of fsrz
     * @param srz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param sx current state sensitivity
     * @param h heavyside vector
     * @param ip sensitivity index
     **/
    virtual void fsrz(realtype *srz, const int ie, const realtype t,
                      const realtype *x, const realtype *p, const realtype *k,
                      const realtype *h, const realtype *sx,
                      const int ip) override {}

    /** model specific implementation of fstau
     * @param stau total derivative of event timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     * @param ie event index
     **/
    virtual void fstau(realtype *stau, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *sx, const int ip,
                       const int ie) override {}

    /** model specific implementation of fsx0
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0(realtype *sx0, const realtype t, const realtype *x0,
                      const realtype *p, const realtype *k,
                      const int ip) override {
        sx0_Nakakuki2010(sx0, t, x0, p, k, ip);
    }

    /** model specific implementation of fsx0_fixedParameters
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0_fixedParameters(realtype *sx0, const realtype t,
                                      const realtype *x0, const realtype *p,
                                      const realtype *k,
                                      const int ip) override {
        sx0_fixedParameters_Nakakuki2010(sx0, t, x0, p, k, ip);
    }

    /** model specific implementation of fsz
     * @param sz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     **/
    virtual void fsz(realtype *sz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h, const realtype *sx,
                     const int ip) override {}

    virtual void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        w_Nakakuki2010( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Nakakuki2010(x0, t, p, k);
    }

    /** model specific implementation of fx0_fixedParameters
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0_fixedParameters(realtype *x0, const realtype t,
                                     const realtype *p,
                                     const realtype *k) override {
        x0_fixedParameters_Nakakuki2010(x0, t, p, k);
    }

    /** model specific implementation for fxdot
     * @param xdot residual function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     **/
    virtual void fxdot(realtype *xdot, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w) override {
        xdot_Nakakuki2010(xdot, t, x, p, k, h, w);
    }

    /** model specific implementation of fy
     * @param y model output at current timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fy(realtype *y, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w) override {
        y_Nakakuki2010(y, t, x, p, k, h, w);
    }

    /** model specific implementation of fz
     * @param z value of event output
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fz(realtype *z, const int ie, const realtype t,
                    const realtype *x, const realtype *p, const realtype *k,
                    const realtype *h) override {}

    

    virtual void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_Nakakuki2010( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Nakakuki2010( total_cl,  x_rdata);
    }


    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>{"V2",
"V3",
"K3",
"V4",
"K4",
"V5",
"K5",
"V6",
"K6",
"KimERK",
"KexERK",
"KimERKP",
"KexERKP",
"KimERKPP",
"KexERKPP",
"V10",
"K10",
"n10",
"p11",
"p12",
"p13",
"V14",
"K14",
"V15",
"K15",
"p16",
"p17",
"KimDUSP",
"KexDUSP",
"KimDUSPP",
"KexDUSPP",
"V20",
"K20",
"V21",
"K21",
"p22",
"p23",
"V24",
"K24",
"V1",
"K1",
"V25",
"K25",
"KimRSKP",
"KexRSKP",
"V27",
"K27",
"V28",
"K28",
"V29",
"K29",
"V30",
"K30",
"V31",
"K31",
"n31",
"p32",
"p33",
"p34",
"V35",
"K35",
"V36",
"K36",
"V37",
"K37",
"p38",
"p39",
"KimFOS",
"KexFOS",
"KimFOSP",
"KexFOSP",
"V42",
"K42",
"V43",
"K43",
"V44",
"K44",
"p45",
"p46",
"p47",
"m47",
"p48",
"p49",
"m49",
"p50",
"p51",
"m51",
"Fct",
"p52",
"m52",
"p53",
"p54",
"m54",
"p55",
"p56",
"m56",
"V57",
"K57",
"n57",
"p58",
"p59",
"p60",
"p61",
"KimF",
"KexF",
"p63",
"KF31",
"nF31",
"K2",
"Vn",
"Vc",
"V101",
"K101",
"V102",
"K102",
"V103",
"K103",
"V104",
"K104",
"V105",
"K105",
"V106",
"K106",
"V107",
"K107",
"V108",
"K108",
"V109",
"K109",
"V110",
"K110",
"V111",
"K111",
"V112",
"K112",
"V113",
"K113",
"V114",
"K114",
"V115",
"K115",};
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>{"EGF",
"HRG",
"A1",
"A1_2",
"A2",
"A2_2",
"A3",
"A3_2",
"DUSPmRNA",
"ERK_c",
"pERK_c",
"ppERK_c",
"F",
"c_FOS_c",
"pc_FOS_c",
"c_FOSmRNA",
"FmRNA",
"Kin",
"Kin_2",
"pMEK",
"MEK",
"DUSP_c",
"pDUSP_c",
"RSK_c",
"pRSK_c",
"RsD",
"RsT",
"CREB_n",
"pCREB_n",
"ERK_n",
"pERK_n",
"ppERK_n",
"Elk1_n",
"pElk1_n",
"FOSn",
"FOSn_2",
"Fn",
"DUSP_n",
"pDUSP_n",
"pDUSP_n_ERK_n",
"pDUSP_n_pERK_n",
"pDUSP_n_ppERK_n",
"DUSP_n_ERK_n",
"DUSP_n_pERK_n",
"DUSP_n_ppERK_n",
"PreDUSPmRNA",
"PreFOSmRNA",
"PreFmRNA",
"pRSK_n",};
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    virtual std::vector<std::string> getFixedParameterNames() const override {
        return std::vector<std::string>{
            };
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    virtual std::vector<std::string> getObservableNames() const override {
        return std::vector<std::string>{"x0",
"x1",
"x2",
"x3",
"x4",
"x5",
"x6",
"x7",
"x8",
"x9",
"x10",
"x11",
"x12",
"x13",
"x14",
"x15",
"x16",
"x17",
"x18",
"x19",
"x20",
"x21",
"x22",
"x23",
"x24",
"x25",
"x26",
"x27",
"x28",
"x29",
"x30",
"x31",
"x32",
"x33",
"x34",
"x35",
"x36",
"x37",
"x38",
"x39",
"x40",
"x41",
"x42",
"x43",
"x44",
"x45",
"x46",
"x47",
"x48",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"V2",
"V3",
"K3",
"V4",
"K4",
"V5",
"K5",
"V6",
"K6",
"KimERK",
"KexERK",
"KimERKP",
"KexERKP",
"KimERKPP",
"KexERKPP",
"V10",
"K10",
"n10",
"p11",
"p12",
"p13",
"V14",
"K14",
"V15",
"K15",
"p16",
"p17",
"KimDUSP",
"KexDUSP",
"KimDUSPP",
"KexDUSPP",
"V20",
"K20",
"V21",
"K21",
"p22",
"p23",
"V24",
"K24",
"V1",
"K1",
"V25",
"K25",
"KimRSKP",
"KexRSKP",
"V27",
"K27",
"V28",
"K28",
"V29",
"K29",
"V30",
"K30",
"V31",
"K31",
"n31",
"p32",
"p33",
"p34",
"V35",
"K35",
"V36",
"K36",
"V37",
"K37",
"p38",
"p39",
"KimFOS",
"KexFOS",
"KimFOSP",
"KexFOSP",
"V42",
"K42",
"V43",
"K43",
"V44",
"K44",
"p45",
"p46",
"p47",
"m47",
"p48",
"p49",
"m49",
"p50",
"p51",
"m51",
"Fct",
"p52",
"m52",
"p53",
"p54",
"m54",
"p55",
"p56",
"m56",
"V57",
"K57",
"n57",
"p58",
"p59",
"p60",
"p61",
"KimF",
"KexF",
"p63",
"KF31",
"nF31",
"K2",
"Vn",
"Vc",
"V101",
"K101",
"V102",
"K102",
"V103",
"K103",
"V104",
"K104",
"V105",
"K105",
"V106",
"K106",
"V107",
"K107",
"V108",
"K108",
"V109",
"K109",
"V110",
"K110",
"V111",
"K111",
"V112",
"K112",
"V113",
"K113",
"V114",
"K114",
"V115",
"K115",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"EGF",
"HRG",
"A1",
"A1_2",
"A2",
"A2_2",
"A3",
"A3_2",
"DUSPmRNA",
"ERK_c",
"pERK_c",
"ppERK_c",
"F",
"c_FOS_c",
"pc_FOS_c",
"c_FOSmRNA",
"FmRNA",
"Kin",
"Kin_2",
"pMEK",
"MEK",
"DUSP_c",
"pDUSP_c",
"RSK_c",
"pRSK_c",
"RsD",
"RsT",
"CREB_n",
"pCREB_n",
"ERK_n",
"pERK_n",
"ppERK_n",
"Elk1_n",
"pElk1_n",
"FOSn",
"FOSn_2",
"Fn",
"DUSP_n",
"pDUSP_n",
"pDUSP_n_ERK_n",
"pDUSP_n_pERK_n",
"pDUSP_n_ppERK_n",
"DUSP_n_ERK_n",
"DUSP_n_pERK_n",
"DUSP_n_ppERK_n",
"PreDUSPmRNA",
"PreFOSmRNA",
"PreFmRNA",
"pRSK_n",};
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getFixedParameterIds() const override {
        return std::vector<std::string>{
            };
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    virtual std::vector<std::string> getObservableIds() const override {
        return std::vector<std::string>{"y0",
"y1",
"y2",
"y3",
"y4",
"y5",
"y6",
"y7",
"y8",
"y9",
"y10",
"y11",
"y12",
"y13",
"y14",
"y15",
"y16",
"y17",
"y18",
"y19",
"y20",
"y21",
"y22",
"y23",
"y24",
"y25",
"y26",
"y27",
"y28",
"y29",
"y30",
"y31",
"y32",
"y33",
"y34",
"y35",
"y36",
"y37",
"y38",
"y39",
"y40",
"y41",
"y42",
"y43",
"y44",
"y45",
"y46",
"y47",
"y48",};
    }

    /**
     * @brief function indicating whether reinitialization of states depending on
     fixed parameters is permissible
     * @return flag inidication whether reinitialization of states depending on
     fixed parameters is permissible
     */
    virtual bool isFixedParameterStateReinitializationAllowed() const override {
        return true;
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return ver amici version string
     */
    virtual const std::string getAmiciVersion() const override {
        return "0.10.7";
    }

    /**
     & @brief returns the amici version that was used to generate the model
     * @return commit amici git commit hash
     */
    virtual const std::string getAmiciCommit() const override {
        return "unknown";
    }

    virtual bool wasPythonGenerated() const override {
        return true;
    }
};

#endif /* _amici_TPL_MODELNAME_h */
