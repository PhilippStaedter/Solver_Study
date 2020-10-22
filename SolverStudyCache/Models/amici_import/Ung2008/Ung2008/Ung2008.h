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
extern void J_Ung2008(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Ung2008(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Ung2008(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Ung2008(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Ung2008(sunindextype *colptrs);
extern void JSparse_rowvals_Ung2008(sunindextype *rowvals);
extern void JSparseB_Ung2008(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Ung2008(sunindextype *colptrs);
extern void JSparseB_rowvals_Ung2008(sunindextype *rowvals);
extern void Jy_Ung2008(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Ung2008(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Ung2008(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Ung2008(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Ung2008(sunindextype *rowvals, int index);
extern void dwdp_Ung2008(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Ung2008(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Ung2008(sunindextype *colptrs);
extern void dwdx_rowvals_Ung2008(sunindextype *rowvals);
extern void dxdotdw_Ung2008(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Ung2008(sunindextype *colptrs);
extern void dxdotdw_rowvals_Ung2008(sunindextype *rowvals);
extern void dxdotdp_Ung2008(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Ung2008(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Ung2008(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Ung2008(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Ung2008(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Ung2008(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Ung2008(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Ung2008(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Ung2008(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Ung2008(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Ung2008(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Ung2008(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Ung2008(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Ung2008(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Ung2008 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Ung2008()
        : amici::Model_ODE(
              194,                                // nx_rdata
              194,                            // nxtrue_rdata
              194,                               // nx_solver
              194,                           // nxtrue_solver
              194,                                      // ny
              194,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              205,                                      // nw
              420,                                   // ndwdx
              64165,                                   // ndwdp
              590,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              956,                                     // nnz
              194,                                     // ubw
              194,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{0.0038, 100.0, 0.02, 10.0, 2.014, 0.2, 3.14, 0.2661, 0.6, 90.0, 0.5838, 0.3, 4.481, 0.2, 3.114, 0.2661, 0.1, 3.0, 0.0214, 10.0, 0.0015, 0.1, 0.045, 10.0, 0.18, 202.9, 0.1434, 0.05, 3.0, 0.06, 10.0, 0.025, 2.734, 0.18, 202.9, 0.1434, 0.05, 1.754, 0.7624, 0.01833, 4.0, 3.5, 0.01833, 4.0, 2.9, 0.033, 3.0, 16.0, 0.033, 3.0, 5.7, 0.2, 71.7, 1.0, 0.8, 14.3, 0.058, 0.5, 0.25, 0.058, 0.6, 7.0, 0.27, 0.5, 5.0, 0.3, 1.667e-05, 0.96, 2.854, 7.6, 0.1, 8.898, 0.426, 0.1, 8.898, 0.426, 0.002, 0.1743, 14.0, 33.72, 0.038, 1.0, 0.595, 3.5, 25.0, 25.0, 1.0, 3.0, 1.0, 3.0, 3.0, 10.0, 0.001, 0.9, 1.0, 0.001, 0.5, 3.0, 3.0, 0.001, 0.033, 1.1, 16.0, 0.5, 5.0, 5.0, 0.1298, 17.0, 0.0214, 10.0, 0.18, 2.029, 0.1434, 0.96, 2.845, 0.262, 0.96, 2.845, 1.205, 0.18, 20.29, 0.18, 20.29, 4.98, 0.262, 0.01, 0.1, 0.96, 2.854, 7.76, 1.0, 10.0, 1.0, 10.0, 2.661, 2.661, 0.2, 3.114, 2.661, 0.2, 3.114, 2.661, 0.2, 3.114, 2.661, 0.2, 3.114, 2.661, 0.2, 3.114, 2.661, 0.9356, 40.0, 10.0, 0.1298, 0.9356, 40.0, 10.0, 0.96, 2.845, 1.205, 0.5, 1.754, 7.624, 0.9356, 40.0, 40.0, 0.0003302, 0.001, 1.0, 1.0, 0.01, 1.2987, 0.1, 0.005, 0.5, 0.005, 0.5, 0.01, 5.0, 0.01, 5.0, 0.001, 0.001, 0.005, 0.033, 1.1, 16.0, 0.005, 10.0, 129.8, 0.05, 1.754, 0.07624, 1.0, 8.898, 0.426, 0.05, 3.0, 0.1, 3.0, 0.5, 3.0, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 0.5, 1.667, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 0.5, 1.667, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 0.5, 1.667, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 0.5, 1.667, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 1.67, 5.0, 0.1002, 0.1002, 1.67, 1.67, 5.0, 0.1002, 0.1002, 1.67, 0.5, 1.667, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 0.5, 1.667, 1.67, 5.0, 1.693, 1.693, 5.0, 1.67, 0.05, 16.67, 0.05, 16.67, 0.05, 16.67, 0.05, 16.67, 1.67, 5.0, 0.1002, 0.1002, 1.67, 1.67, 5.0, 0.1002, 0.1002, 1.67, 1.67, 5.0, 0.1002, 0.1002, 1.67, 1.67, 5.0, 0.1002, 0.1002, 1.67, 0.5, 1.667, 1.205, 0.5, 1.667, 1.205, 0.5, 1.667, 1.205, 0.5, 1.667, 1.205},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(194, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Ung2008(*this);
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
        J_Ung2008(J, t, x, p, k, h, w, dwdx);
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
        JB_Ung2008(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_Ung2008(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Ung2008( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Ung2008(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Ung2008(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Ung2008( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Ung2008(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Ung2008(rowvals);
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
        Jy_Ung2008(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_Ung2008(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_Ung2008(dsigmaydp, t, p, k, ip);
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
        dJydy_Ung2008( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Ung2008(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Ung2008(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Ung2008( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Ung2008( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Ung2008( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Ung2008(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Ung2008(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Ung2008( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_Ung2008(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_Ung2008(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_Ung2008(sigmay, t, p, k);
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
        sx0_Ung2008(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_Ung2008(sx0, t, x0, p, k, ip);
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
        w_Ung2008( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Ung2008(x0, t, p, k);
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
        x0_fixedParameters_Ung2008(x0, t, p, k);
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
        xdot_Ung2008(xdot, t, x, p, k, h, w);
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
        y_Ung2008(y, t, x, p, k, h, w);
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
        x_solver_Ung2008( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Ung2008( total_cl,  x_rdata);
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
        return std::vector<std::string>{"EGF",
"EGFR",
"EGF-EGFR",
"EGF-EGFR-2",
"EGF-pEGFR-2",
"SHP",
"EGF-pEGFR-2-SHP",
"Shc",
"EGF-pEGFR-2-Shc",
"EGF-pEGFR-2-pShc",
"pShc",
"pShc-SHP",
"Grb2",
"EGF-pEGFR-2-pShc-Grb2",
"SOS",
"EGF-pEGFR-2-pShc-Grb2-SOS",
"Grb2-SOS",
"RasGDP",
"EGF-pEGFR-2-pShc-Grb2-SOS-RasGDP",
"RasGTP",
"EGF-pEGFR-2-Grb2",
"EGF-pEGFR-2-Grb2-SOS",
"EGF-pEGFR-2-Grb2-SOS-RasGDP",
"Raf",
"Raf-RasGTP",
"pRaf",
"MEK",
"pRaf-MEK",
"pMEK",
"pRaf-pMEK",
"ppMEK",
"ERK",
"ppMEK-ERK",
"pERK",
"ppMEK-pERK",
"ppERK",
"Pase",
"pRaf-Pase",
"PP2A",
"ppMEK-PP2A",
"pMEK-PP2A",
"MKP3",
"ppERK-MKP3",
"pERK-MKP3",
"RasGAP",
"RasGTP-RasGAP",
"ppERK-EGF-pEGFR-2-pShc-Grb2-SOS",
"pSOS",
"ppERK-EGF-pEGFR-2-Grb2-SOS",
"PI3K",
"EGF-pEGFR-2-PI3K",
"EGF-pEGFF-2",
"pPI3K",
"TP4",
"pPI3K-TP4",
"PIP2",
"pPI3K-PIP2",
"PIP3",
"Akt",
"Akt-PIP3",
"PDK1",
"Akt-PIP3-PDK1",
"pAkt-PIP3",
"pAkt",
"Takt",
"pAkt-PIP3-Takt",
"pRaf-pAkt-PIP3",
"ppRaf",
"pROK",
"PTEN",
"pROK-PTEN",
"pPTEN",
"pPTEN-PIP3",
"RacGEF",
"PIP3-RacGEF",
"RacGDP",
"PIP3-RacGEF-RacGDP",
"RacGTP",
"RhoGDI",
"RhoGDI-RacGDP",
"RacGAP",
"RacGTP-RacGAP",
"RhoGDP",
"RhoGDP-RhoGDI",
"pRhoGEF",
"RhoGDP-pRhoGEF",
"RhoGTP",
"EGF-pEGFR-2-RasGAP",
"EGF-pEGFR-2-RasGAP-RasGTP",
"EGF-pEGFR2-RasGAP",
"SHP2",
"EGF-pEGFR-2-pShc-Grb2-SHP2",
"EGF-pEGFR-2-Grb2-SHP2",
"EGF-pEGFR-2-pShc-Grb2-SHP2-pRhoGEF",
"RhoGEF",
"pRhoGAP",
"EGF-pEGFR-2-pShc-Grb2-SHP2-pRhoGAP",
"RhoGAP",
"EGF-pEGFR-2-Grb2-SHP2-pRhoGEF",
"EGF-pEGFR-2-Grb2-SHP2-pRhoGAP",
"EGF-pEGFR-2-RasGAP-SHP2",
"pSrc",
"pSrc-RhoGEF",
"pSrc-RhoGAP",
"pRhoGAP-RhoGTP",
"ROK",
"RhoGTP-ROK",
"Src",
"EGF-pEGFR-2-Src",
"EGF-pEGFR-2-pSrc",
"EGF-pEGRF-2",
"TP7",
"pSrc-TP7",
"Src-TP7",
"Cbl-CIN85",
"EGF-pEGFR-2-pShc-Grb2-SOS-Cbl-CIN85",
"EGF-pEGFR-2-Grb2-SOS-Cbl-CIN85",
"EPn",
"EGF-pEGFR-2-pShc-Grb2-SOS-Cbl-CIN85-EPn",
"EGF-pEGFR-2-Grb2-SOS-Cbl-CIN85-EPn",
"EGF-pEGFR-2-degrade",
"pShc-Grb2-SOS",
"Pro-EGFR",
"pROK-EPn",
"pEPn",
"MPase",
"pEPn-MPase",
"pEPn-Mpase",
"Ras-GTP-RhoGEF",
"ppERK-pROK",
"MEKK1abcdef",
"Grb2-MEKK1abcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abcdef",
"EGF-pEGFR-2-Grb2-MEKK1abcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abpMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1abMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1abpMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbpMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbppMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbcdef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbpMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbppMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abcdefRasGTP",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abMEKcdefRasGTP",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abpMEKcdefRasGTP",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcdefRasGTP",
"EGF-pEGFR-2-Grb2-MEKK1abcdefRasGTP",
"EGF-pEGFR-2-Grb2-MEKK1abMEKcdefRasGTP",
"EGF-pEGFR-2-Grb2-MEKK1abpMEKcdefRasGTP",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcdefRasGTP",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcERKdefRasGTP",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcpERKdefRasGTP",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcppERKdefRasGTP",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcERKdefRasGTP",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcpERKdefRasGTP",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcppERKdefRasGTP",
"EGF-pEGFR-2-pShc-Grb2-MEKK1apRafbcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1apRafbMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1apRafbpMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1apRafbppMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1apRafbcdef",
"EGF-pEGFR-2-Grb2-MEKK1apRafbMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1apRafbpMEKcdef",
"EGF-pEGFR-2-Grb2-MEKK1apRafbppMEKcdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abMEKcdRhoGTPef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcdRhoGTPef",
"EGF-pEGFR-2-Grb2-MEKK1abMEKcdRhoGTPef",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcdRhoGTPef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbMEKcdRhoGTPef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbcdRhoGTPef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbMEKcdRhoGTPef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbcdRhoGTPef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcERKdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcpERKdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1abppMEKcppERKdef",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcERKdef",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcpERKdef",
"EGF-pEGFR-2-Grb2-MEKK1abppMEKcppERKdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbppMEKcERKdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbppMEKcpERKdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbppMEKcppERKdef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbppMEKcERKdef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbppMEKcpERKdef",
"EGF-pEGFR-2-Grb2-MEKK1aRafbppMEKcppERKdef",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbppMEKcdRhoGTPepRhoGAPf",
"EGF-pEGFR-2-Grb2-MEKK1aRafbppMEKcdRhoGTPepRhoGAPf",
"EGF-pEGFR-2-pShc-Grb2-MEKK1aRafbcdRhoGTPepRhoGAPf",
"EGF-pEGFR-2-Grb2-MEKK1aRafbcdRhoGTPepRhoGAPf",};
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
"x48",
"x49",
"x50",
"x51",
"x52",
"x53",
"x54",
"x55",
"x56",
"x57",
"x58",
"x59",
"x60",
"x61",
"x62",
"x63",
"x64",
"x65",
"x66",
"x67",
"x68",
"x69",
"x70",
"x71",
"x72",
"x73",
"x74",
"x75",
"x76",
"x77",
"x78",
"x79",
"x80",
"x81",
"x82",
"x83",
"x84",
"x85",
"x86",
"x87",
"x88",
"x89",
"x90",
"x91",
"x92",
"x93",
"x94",
"x95",
"x96",
"x97",
"x98",
"x99",
"x100",
"x101",
"x102",
"x103",
"x104",
"x105",
"x106",
"x107",
"x108",
"x109",
"x110",
"x111",
"x112",
"x113",
"x114",
"x115",
"x116",
"x117",
"x118",
"x119",
"x120",
"x121",
"x122",
"x123",
"x124",
"x125",
"x126",
"x127",
"x128",
"x129",
"x130",
"x131",
"x132",
"x133",
"x134",
"x135",
"x136",
"x137",
"x138",
"x139",
"x140",
"x141",
"x142",
"x143",
"x144",
"x145",
"x146",
"x147",
"x148",
"x149",
"x150",
"x151",
"x152",
"x153",
"x154",
"x155",
"x156",
"x157",
"x158",
"x159",
"x160",
"x161",
"x162",
"x163",
"x164",
"x165",
"x166",
"x167",
"x168",
"x169",
"x170",
"x171",
"x172",
"x173",
"x174",
"x175",
"x176",
"x177",
"x178",
"x179",
"x180",
"x181",
"x182",
"x183",
"x184",
"x185",
"x186",
"x187",
"x188",
"x189",
"x190",
"x191",
"x192",
"x193",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"reaction_0_k2",
"reaction_0_k1",
"reaction_1_k2",
"reaction_1_k1",
"reaction_2_k1",
"reaction_3_k2",
"reaction_3_k1",
"reaction_4_k1",
"reaction_5_k2",
"reaction_5_k1",
"reaction_6_k1",
"reaction_7_k2",
"reaction_7_k1",
"reaction_8_k2",
"reaction_8_k1",
"reaction_9_k1",
"reaction_10_k2",
"reaction_10_k1",
"reaction_11_k2",
"reaction_11_k1",
"reaction_12_k2",
"reaction_12_k1",
"reaction_13_k2",
"reaction_13_k1",
"reaction_14_k2",
"reaction_14_k1",
"reaction_15_k1",
"reaction_16_k2",
"reaction_16_k1",
"reaction_17_k2",
"reaction_17_k1",
"reaction_18_k2",
"reaction_18_k1",
"reaction_19_k2",
"reaction_19_k1",
"reaction_20_k1",
"reaction_21_k2",
"reaction_21_k1",
"reaction_22_k1",
"reaction_23_k2",
"reaction_23_k1",
"reaction_24_k1",
"reaction_25_k2",
"reaction_25_k1",
"reaction_26_k1",
"reaction_27_k2",
"reaction_27_k1",
"reaction_28_k1",
"reaction_29_k2",
"reaction_29_k1",
"reaction_30_k1",
"reaction_31_k2",
"reaction_31_k1",
"reaction_32_k1",
"reaction_33_k2",
"reaction_33_k1",
"reaction_34_k1",
"reaction_35_k2",
"reaction_35_k1",
"reaction_36_k1",
"reaction_37_k2",
"reaction_37_k1",
"reaction_38_k1",
"reaction_39_k2",
"reaction_39_k1",
"reaction_40_k1",
"reaction_41_k1",
"reaction_42_k2",
"reaction_42_k1",
"reaction_43_k1",
"reaction_44_k2",
"reaction_44_k1",
"reaction_45_k1",
"reaction_46_k2",
"reaction_46_k1",
"reaction_47_k1",
"reaction_48_k1",
"reaction_49_k2",
"reaction_49_k1",
"reaction_50_k1",
"reaction_51_k2",
"reaction_51_k1",
"reaction_52_k1",
"reaction_53_k2",
"reaction_53_k1",
"reaction_54_k1",
"reaction_55_k2",
"reaction_55_k1",
"reaction_56_k2",
"reaction_56_k1",
"reaction_57_k1",
"reaction_58_k2",
"reaction_58_k1",
"reaction_59_k2",
"reaction_59_k1",
"reaction_60_k1",
"reaction_61_k2",
"reaction_61_k1",
"reaction_62_k1",
"reaction_63_k1",
"reaction_64_k2",
"reaction_64_k1",
"reaction_65_k1",
"reaction_66_k2",
"reaction_66_k1",
"reaction_67_k1",
"reaction_68_k1",
"reaction_69_k1",
"reaction_70_k2",
"reaction_70_k1",
"reaction_71_k2",
"reaction_71_k1",
"reaction_72_k1",
"reaction_73_k2",
"reaction_73_k1",
"reaction_74_k1",
"reaction_75_k2",
"reaction_75_k1",
"reaction_76_k1",
"reaction_77_k2",
"reaction_77_k1",
"reaction_78_k2",
"reaction_78_k1",
"reaction_79_k1",
"reaction_80_k1",
"reaction_81_k2",
"reaction_81_k1",
"reaction_82_k2",
"reaction_82_k1",
"reaction_83_k1",
"reaction_84_k2",
"reaction_84_k1",
"reaction_85_k2",
"reaction_85_k1",
"reaction_86_k1",
"reaction_87_k1",
"reaction_88_k2",
"reaction_88_k1",
"reaction_89_k1",
"reaction_90_k2",
"reaction_90_k1",
"reaction_91_k1",
"reaction_92_k2",
"reaction_92_k1",
"reaction_93_k1",
"reaction_94_k2",
"reaction_94_k1",
"reaction_95_k1",
"reaction_96_k2",
"reaction_96_k1",
"reaction_97_k1",
"reaction_98_k2",
"reaction_98_k1",
"reaction_99_k1",
"reaction_100_k1",
"reaction_101_k2",
"reaction_101_k1",
"reaction_102_k1",
"reaction_103_k2",
"reaction_103_k1",
"reaction_104_k1",
"reaction_105_k2",
"reaction_105_k1",
"reaction_106_k1",
"reaction_107_k2",
"reaction_107_k1",
"reaction_108_k1",
"reaction_109_k2",
"reaction_109_k1",
"reaction_110_k2",
"reaction_110_k1",
"reaction_111_k1",
"reaction_112_k2",
"reaction_112_k1",
"reaction_113_k2",
"reaction_113_k1",
"reaction_114_k2",
"reaction_114_k1",
"reaction_115_k2",
"reaction_115_k1",
"reaction_116_k2",
"reaction_116_k1",
"reaction_117_k1",
"reaction_118_k1",
"reaction_119_k1",
"reaction_120_k2",
"reaction_120_k1",
"reaction_121_k1",
"reaction_122_k2",
"reaction_122_k1",
"reaction_123_k1",
"reaction_124_k2",
"reaction_124_k1",
"reaction_125_k1",
"reaction_126_k2",
"reaction_126_k1",
"reaction_127_k1",
"reaction_128_k2",
"reaction_128_k1",
"reaction_129_k2",
"reaction_129_k1",
"reaction_130_k2",
"reaction_130_k1",
"reaction_131_k2",
"reaction_131_k1",
"reaction_132_k1",
"reaction_133_k1",
"reaction_134_k2",
"reaction_134_k1",
"reaction_135_k2",
"reaction_135_k1",
"reaction_136_k1",
"reaction_137_k1",
"reaction_138_k2",
"reaction_138_k1",
"reaction_139_k2",
"reaction_139_k1",
"reaction_140_k2",
"reaction_140_k1",
"reaction_141_k1",
"reaction_142_k1",
"reaction_143_k2",
"reaction_143_k1",
"reaction_144_k2",
"reaction_144_k1",
"reaction_145_k2",
"reaction_145_k1",
"reaction_146_k1",
"reaction_147_k1",
"reaction_148_k2",
"reaction_148_k1",
"reaction_149_k2",
"reaction_149_k1",
"reaction_150_k2",
"reaction_150_k1",
"reaction_151_k1",
"reaction_152_k1",
"reaction_153_k2",
"reaction_153_k1",
"reaction_154_k2",
"reaction_154_k1",
"reaction_155_k2",
"reaction_155_k1",
"reaction_156_k1",
"reaction_157_k1",
"reaction_158_k2",
"reaction_158_k1",
"reaction_159_k2",
"reaction_159_k1",
"reaction_160_k1",
"reaction_161_k1",
"reaction_162_k1",
"reaction_163_k2",
"reaction_163_k1",
"reaction_164_k1",
"reaction_165_k1",
"reaction_166_k1",
"reaction_167_k2",
"reaction_167_k1",
"reaction_168_k2",
"reaction_168_k1",
"reaction_169_k1",
"reaction_170_k1",
"reaction_171_k2",
"reaction_171_k1",
"reaction_172_k2",
"reaction_172_k1",
"reaction_173_k2",
"reaction_173_k1",
"reaction_174_k1",
"reaction_175_k1",
"reaction_176_k2",
"reaction_176_k1",
"reaction_177_k2",
"reaction_177_k1",
"reaction_178_k2",
"reaction_178_k1",
"reaction_179_k2",
"reaction_179_k1",
"reaction_180_k2",
"reaction_180_k1",
"reaction_181_k2",
"reaction_181_k1",
"reaction_182_k1",
"reaction_183_k1",
"reaction_184_k1",
"reaction_185_k2",
"reaction_185_k1",
"reaction_186_k1",
"reaction_187_k1",
"reaction_188_k1",
"reaction_189_k2",
"reaction_189_k1",
"reaction_190_k1",
"reaction_191_k1",
"reaction_192_k1",
"reaction_193_k2",
"reaction_193_k1",
"reaction_194_k1",
"reaction_195_k1",
"reaction_196_k1",
"reaction_197_k2",
"reaction_197_k1",
"reaction_198_k1",
"reaction_199_k2",
"reaction_199_k1",
"reaction_200_k1",
"reaction_201_k2",
"reaction_201_k1",
"reaction_202_k1",
"reaction_203_k2",
"reaction_203_k1",
"reaction_204_k1",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"species_0",
"species_1",
"species_2",
"species_3",
"species_4",
"species_5",
"species_6",
"species_7",
"species_8",
"species_9",
"species_10",
"species_11",
"species_12",
"species_13",
"species_14",
"species_15",
"species_16",
"species_17",
"species_18",
"species_19",
"species_20",
"species_21",
"species_22",
"species_23",
"species_24",
"species_25",
"species_26",
"species_27",
"species_28",
"species_29",
"species_30",
"species_31",
"species_32",
"species_33",
"species_34",
"species_35",
"species_36",
"species_37",
"species_38",
"species_39",
"species_40",
"species_41",
"species_42",
"species_43",
"species_44",
"species_45",
"species_46",
"species_47",
"species_48",
"species_49",
"species_50",
"species_51",
"species_52",
"species_53",
"species_54",
"species_55",
"species_56",
"species_57",
"species_58",
"species_59",
"species_60",
"species_61",
"species_62",
"species_63",
"species_64",
"species_65",
"species_66",
"species_67",
"species_68",
"species_69",
"species_70",
"species_71",
"species_72",
"species_73",
"species_74",
"species_75",
"species_76",
"species_77",
"species_78",
"species_79",
"species_80",
"species_81",
"species_82",
"species_83",
"species_84",
"species_85",
"species_86",
"species_87",
"species_88",
"species_89",
"species_90",
"species_91",
"species_92",
"species_93",
"species_94",
"species_95",
"species_96",
"species_97",
"species_98",
"species_99",
"species_100",
"species_101",
"species_102",
"species_103",
"species_104",
"species_105",
"species_106",
"species_107",
"species_108",
"species_109",
"species_110",
"species_111",
"species_112",
"species_113",
"species_114",
"species_115",
"species_116",
"species_117",
"species_118",
"species_119",
"species_120",
"species_121",
"species_122",
"species_123",
"species_124",
"species_125",
"species_126",
"species_127",
"species_128",
"species_129",
"species_130",
"species_131",
"species_132",
"species_133",
"species_134",
"species_135",
"species_136",
"species_137",
"species_138",
"species_139",
"species_140",
"species_141",
"species_142",
"species_143",
"species_144",
"species_145",
"species_146",
"species_147",
"species_148",
"species_149",
"species_150",
"species_151",
"species_152",
"species_153",
"species_154",
"species_155",
"species_156",
"species_157",
"species_158",
"species_159",
"species_160",
"species_161",
"species_162",
"species_163",
"species_164",
"species_165",
"species_166",
"species_167",
"species_168",
"species_169",
"species_170",
"species_171",
"species_172",
"species_173",
"species_174",
"species_175",
"species_176",
"species_177",
"species_178",
"species_179",
"species_180",
"species_181",
"species_182",
"species_183",
"species_184",
"species_185",
"species_186",
"species_187",
"species_188",
"species_189",
"species_190",
"species_191",
"species_192",
"species_193",};
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
"y48",
"y49",
"y50",
"y51",
"y52",
"y53",
"y54",
"y55",
"y56",
"y57",
"y58",
"y59",
"y60",
"y61",
"y62",
"y63",
"y64",
"y65",
"y66",
"y67",
"y68",
"y69",
"y70",
"y71",
"y72",
"y73",
"y74",
"y75",
"y76",
"y77",
"y78",
"y79",
"y80",
"y81",
"y82",
"y83",
"y84",
"y85",
"y86",
"y87",
"y88",
"y89",
"y90",
"y91",
"y92",
"y93",
"y94",
"y95",
"y96",
"y97",
"y98",
"y99",
"y100",
"y101",
"y102",
"y103",
"y104",
"y105",
"y106",
"y107",
"y108",
"y109",
"y110",
"y111",
"y112",
"y113",
"y114",
"y115",
"y116",
"y117",
"y118",
"y119",
"y120",
"y121",
"y122",
"y123",
"y124",
"y125",
"y126",
"y127",
"y128",
"y129",
"y130",
"y131",
"y132",
"y133",
"y134",
"y135",
"y136",
"y137",
"y138",
"y139",
"y140",
"y141",
"y142",
"y143",
"y144",
"y145",
"y146",
"y147",
"y148",
"y149",
"y150",
"y151",
"y152",
"y153",
"y154",
"y155",
"y156",
"y157",
"y158",
"y159",
"y160",
"y161",
"y162",
"y163",
"y164",
"y165",
"y166",
"y167",
"y168",
"y169",
"y170",
"y171",
"y172",
"y173",
"y174",
"y175",
"y176",
"y177",
"y178",
"y179",
"y180",
"y181",
"y182",
"y183",
"y184",
"y185",
"y186",
"y187",
"y188",
"y189",
"y190",
"y191",
"y192",
"y193",};
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
