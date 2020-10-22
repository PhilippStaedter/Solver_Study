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
extern void J_Holzhutter2004(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Holzhutter2004(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Holzhutter2004(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Holzhutter2004(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Holzhutter2004(sunindextype *colptrs);
extern void JSparse_rowvals_Holzhutter2004(sunindextype *rowvals);
extern void JSparseB_Holzhutter2004(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Holzhutter2004(sunindextype *colptrs);
extern void JSparseB_rowvals_Holzhutter2004(sunindextype *rowvals);
extern void Jy_Holzhutter2004(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Holzhutter2004(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Holzhutter2004(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Holzhutter2004(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Holzhutter2004(sunindextype *rowvals, int index);
extern void dwdp_Holzhutter2004(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Holzhutter2004(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Holzhutter2004(sunindextype *colptrs);
extern void dwdx_rowvals_Holzhutter2004(sunindextype *rowvals);
extern void dxdotdw_Holzhutter2004(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Holzhutter2004(sunindextype *colptrs);
extern void dxdotdw_rowvals_Holzhutter2004(sunindextype *rowvals);
extern void dxdotdp_Holzhutter2004(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Holzhutter2004(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Holzhutter2004(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Holzhutter2004(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Holzhutter2004(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Holzhutter2004(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Holzhutter2004(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Holzhutter2004(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Holzhutter2004(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Holzhutter2004(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Holzhutter2004(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Holzhutter2004(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Holzhutter2004(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Holzhutter2004(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Holzhutter2004 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Holzhutter2004()
        : amici::Model_ODE(
              45,                                // nx_rdata
              45,                            // nxtrue_rdata
              45,                               // nx_solver
              45,                           // nxtrue_solver
              45,                                      // ny
              45,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              38,                                      // nw
              133,                                   // ndwdx
              6308,                                   // ndwdp
              113,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              310,                                     // nnz
              45,                                     // ubw
              45,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{0.54, 6.9, 1.0, 1.7, 33.6, 3.44, 2.7, 0.0045, 1.03, 3900.0, 1.14, 33.2, 1.44, 15.8, 0.1, 1.0, 0.071, 0.182, 0.3925, 935.0, 0.033, 0.44, 0.01, 0.001072, 0.068, 0.1, 100000.0, 239.0, 0.176, 0.0364, 0.1906, 0.0572, 0.114, 0.0071, 98.91000366, 0.428, 0.838, 0.0407, 5456.600098, 0.0035, 0.0083, 0.000192, 3.9, 0.005, 0.05, 4300.0, 1.2, 0.48, 1455.0, 0.002, 0.35, 5000.0, 0.04, 100000.0, 76000.0, 0.2, 100000.0, 0.53, 1.0, 5.0, 0.145, 2000.0, 1.0, 1.0, 1.7, 1500.0, 0.005, 3.39, 19.0, 0.474, 0.225, 13790.0, 570.0, 9090.0, 2800000.0, 14181.8, 243.4, 1.68, 0.11, 0.25, 0.08, 0.09, 1380.0, 2.289, 0.00312, 0.749, 2000.0, 0.00367, 0.0667, 162.0, 0.0045, 0.058, 0.154, 0.12, 141.7, 0.018, 0.01, 1575.0, 1.04, 0.07, 20.0, 0.00852, 0.0652, 90.0, 0.03, 0.5, 0.19, 2.7, 4634.0, 2.2, 0.78, 3.0, 730.0, 48.8, 0.00496, 0.41139, 12.432, 0.00774, 0.3055, 0.4177, 1.05, 23.5, 2.524, 0.006095, 0.8683, 0.1733, 0.4653, 0.04765, 0.00823, 1.05, 27.2, 0.57, 0.03, 100000.0, 1.1, 0.215, 0.0003, 0.0287, 0.0548, 0.122, 0.3055, 0.00184, 1.2, 23.5, 1.0, 100.0, 1.0, 10000.0, 1.0, 10000.0, 0.072, 10000000.0, 0.76, 10000000.0, 16.64, 10000000.0, 1.667, 10000000.0, 0.0002, 10000000.0, 1e-05, 10000000.0, 1e-05, 10000000.0, 0.0002, 10000000.0},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(45, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Holzhutter2004(*this);
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
        J_Holzhutter2004(J, t, x, p, k, h, w, dwdx);
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
        JB_Holzhutter2004(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_Holzhutter2004(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Holzhutter2004( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Holzhutter2004(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Holzhutter2004(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Holzhutter2004( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Holzhutter2004(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Holzhutter2004(rowvals);
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
        Jy_Holzhutter2004(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_Holzhutter2004(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_Holzhutter2004(dsigmaydp, t, p, k, ip);
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
        dJydy_Holzhutter2004( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Holzhutter2004(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Holzhutter2004(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Holzhutter2004( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Holzhutter2004( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Holzhutter2004( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Holzhutter2004(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Holzhutter2004(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Holzhutter2004( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_Holzhutter2004(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_Holzhutter2004(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_Holzhutter2004(sigmay, t, p, k);
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
        sx0_Holzhutter2004(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_Holzhutter2004(sx0, t, x0, p, k, ip);
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
        w_Holzhutter2004( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Holzhutter2004(x0, t, p, k);
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
        x0_fixedParameters_Holzhutter2004(x0, t, p, k);
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
        xdot_Holzhutter2004(xdot, t, x, p, k, h, w);
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
        y_Holzhutter2004(y, t, x, p, k, h, w);
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
        x_solver_Holzhutter2004( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Holzhutter2004( total_cl,  x_rdata);
    }


    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>{"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",};
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>{"Glucose in",
"MgATP",
"Glucose 6-phosphate",
"MgADP",
"Fructose 6-phosphate",
"Fructose 1,6-phosphate",
"Glyceraldehyde 3-phosphate",
"Dihydroxyacetone phosphate",
"Phosphate",
"NAD",
"1,3-Bisphospho-D-glycerate",
"NADH",
"3-Phospho-D-glycerate",
"2,3-Bisphospho-D-glycerate",
"2-Phospho-D-glycerate",
"Phosphoenolpyruvate",
"Pyruvate",
"Lactate",
"NADPH",
"NADP",
"AMP",
"ADP",
"Phospho-D-glucono-1,5-lactone",
"Ribulose 5-phosphate",
"Oxidized Glutathione",
"Reduced Glutathione",
"Xylulose 5-phosphate",
"Ribose 5-phosphate",
"Sedoheptulose 7-phosphate",
"Erythrose 4-phosphate",
"MgAMP",
"ATP",
"Mg",
"MgGri23P2",
"Protein1 bound NADP",
"Protein1",
"Protein1 bound NADPH",
"Protein2 bound NADP",
"Protein2",
"Protein2 bound NADPH",
"PRPP",
"External Lactate",
"External Pyruvate",
"Glucose outside",
"Phosphate external",};
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
"x44",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"vGLT_alfav0",
"vGLT_KMinv0",
"vGLT_Keqv0",
"vGLT_KMoutv0",
"vGLT_Vmaxv0",
"vHEX_KMg23P2Gv1",
"vHEX_K23P2Gv1",
"vHEX_KGlc6Pv1",
"vHEX_KMgv1",
"vHEX_Keqv1",
"vHEX_KMgATPMgv1",
"vHEX_Vmax2v1",
"vHEX_KMgATPv1",
"vHEX_Vmax1v1",
"vHEX_KMGlcv1",
"vHEX_Inhibv1",
"vGPI_KFru6Pv2",
"vGPI_KGlc6Pv2",
"vGPI_Keqv2",
"vGPI_Vmaxv2",
"vPFK_KAMPv3",
"vPFK_KMgv3",
"vPFK_KATPv3",
"vPFK_L0v3",
"vPFK_KMgATPv3",
"vPFK_KFru6Pv3",
"vPFK_Keqv3",
"vPFK_Vmaxv3",
"vALD_KiiGraPv4",
"vALD_KDHAPv4",
"vALD_KGraPv4",
"vALD_KiGraPv4",
"vALD_Keqv4",
"vALD_KFru16P2v4",
"vALD_Vmaxv4",
"vTPI_KGraPv5",
"vTPI_KDHAPv5",
"vTPI_Keqv5",
"vTPI_Vmaxv5",
"vGAPDH_K13P2Gv6",
"vGAPDH_KNADHv6",
"vGAPDH_Keqv6",
"vGAPDH_KPv6",
"vGAPDH_KGraPv6",
"vGAPDH_KNADv6",
"vGAPDH_Vmaxv6",
"vPGK_K3PGv7",
"vPGK_KMgATPv7",
"vPGK_Keqv7",
"vPGK_K13P2Gv7",
"vPGK_KMgADPv7",
"vPGK_Vmaxv7",
"vBPGM_K23P2Gv8",
"vBPGM_Keqv8",
"vBPGM_kDPGMv8",
"vBPGP_K23P2Gv9",
"vBPGP_Keqv9",
"vBPGP_Vmaxv9",
"vPGM_K2PGv10",
"vPGM_K3PGv10",
"vPGM_Keqv10",
"vPGM_Vmaxv10",
"vENO_KPEPv11",
"vENO_K2PGv11",
"vENO_Keqv11",
"vENO_Vmaxv11",
"vPK_KFru16P2v12",
"vPK_KATPv12",
"vPK_L0v12",
"vPK_KMgADPv12",
"vPK_KPEPv12",
"vPK_Keqv12",
"vPK_Vmaxv12",
"vLDHNADH_Keqv13",
"vLDHNADH_Vmaxv13",
"vLDHNADPH_Keqv14",
"vLDHNADPH_kLDHv14",
"vATPase_kATPasev15",
"vAK_KADPv16",
"vAK_Keqv16",
"vAK_KAMPv16",
"vAK_KATPv16",
"vAK_Vmaxv16",
"vG6PDH_KPGA23v17",
"vG6PDH_KNADPHv17",
"vG6PDH_KATPv17",
"vG6PDH_Keqv17",
"vG6PDH_KNADPv17",
"vG6PDH_KG6Pv17",
"vG6PDH_Vmaxv17",
"vPGLDH_KNADPHv18",
"vPGLDH_K6PG2v18",
"vPGLDH_KATPv18",
"vPGLDH_KPGA23v18",
"vPGLDH_Keqv18",
"vPGLDH_KNADPv18",
"vPGLDH_K6PG1v18",
"vPGLDH_Vmaxv18",
"vGSSGRD_Keqv19",
"vGSSGRD_KNADPv19",
"vGSSGRD_KGSHv19",
"vGSSGRD_KNADPHv19",
"vGSSGRD_KGSSGv19",
"vGSSGRD_Vmaxv19",
"vGSHox_Kv20",
"vRibPepi_KX5Pv21",
"vRibPepi_KRu5Pv21",
"vRibPepi_Keqv21",
"vRibPepi_Vmaxv21",
"vRibPiso_KR5Pv22",
"vRibPiso_KRu5Pv22",
"vRibPiso_Keqv22",
"vRibPiso_Vmaxv22",
"vTrKet1_K7v23",
"vTrKet1_K4v23",
"vTrKet1_K5v23",
"vTrKet1_K3v23",
"vTrKet1_K6v23",
"vTrKet1_K2v23",
"vTrKet1_K1v23",
"vTrKet1_Keqv23",
"vTrKet1_Vmaxv23",
"vTrAld_K7v24",
"vTrAld_K4v24",
"vTrAld_K5v24",
"vTrAld_K3v24",
"vTrAld_K6v24",
"vTrAld_K2v24",
"vTrAld_K1v24",
"vTrAld_Keqv24",
"vTrAld_Vmaxv24",
"vPPRPPS_KR5Pv25",
"vPPRPPS_KATPv25",
"vPPRPPS_Keqv25",
"vPPRPPS_Vmaxv25",
"vTrKet2_K7v26",
"vTrKet2_K4v26",
"vTrKet2_K5v26",
"vTrKet2_K3v26",
"vTrKet2_K6v26",
"vTrKet2_K2v26",
"vTrKet2_K1v26",
"vTrKet2_Keqv26",
"vTrKet2_Vmaxv26",
"vPhiexch_Keqv27",
"vPhiexch_Vmaxv27",
"vLacexch_Keqv28",
"vLacexch_Vmaxv28",
"vPyrexch_Keqv29",
"vPyrexch_Vmaxv29",
"vMgATP_KdATP",
"vMgATP_EqMult",
"vMgADP_KdADP",
"vMgADP_EqMult",
"vMgAMP_KdAMP",
"vMgAMP_EqMult",
"vMgGri23P2_Kd23P2G",
"vMgGri23P2_EqMult",
"vP1NADP_Kd1",
"vP1NADP_EqMult",
"vP1NADPH_Kd3",
"vP1NADPH_EqMult",
"vP2NADP_Kd2",
"vP2NADP_EqMult",
"vP2NADPH_Kd4",
"vP2NADPH_EqMult",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"Glcin",
"MgATP",
"Glc6P",
"MgADP",
"Fru6P",
"Fru16P2",
"GraP",
"DHAP",
"Phi",
"NAD",
"Gri13P2",
"NADH",
"Gri3P",
"Gri23P2f",
"Gri2P",
"PEP",
"Pyr",
"Lac",
"NADPHf",
"NADPf",
"AMPf",
"ADPf",
"GlcA6P",
"Rul5P",
"GSSG",
"GSH",
"Xul5P",
"Rib5P",
"Sed7P",
"E4P",
"MgAMP",
"ATPf",
"Mgf",
"MgGri23P2",
"P1NADP",
"P1f",
"P1NADPH",
"P2NADP",
"P2f",
"P2NADPH",
"PRPP",
"Lacex",
"Pyrex",
"Glcout",
"Phiex",};
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
"y44",};
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
