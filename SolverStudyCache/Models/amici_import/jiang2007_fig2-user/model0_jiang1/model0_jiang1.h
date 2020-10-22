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
extern void J_model0_jiang1(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_model0_jiang1(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_model0_jiang1(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_model0_jiang1(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_model0_jiang1(sunindextype *colptrs);
extern void JSparse_rowvals_model0_jiang1(sunindextype *rowvals);
extern void JSparseB_model0_jiang1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_model0_jiang1(sunindextype *colptrs);
extern void JSparseB_rowvals_model0_jiang1(sunindextype *rowvals);
extern void Jy_model0_jiang1(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_model0_jiang1(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_model0_jiang1(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_model0_jiang1(sunindextype *colptrs, int index);
extern void dJydy_rowvals_model0_jiang1(sunindextype *rowvals, int index);
extern void dwdp_model0_jiang1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_model0_jiang1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_model0_jiang1(sunindextype *colptrs);
extern void dwdx_rowvals_model0_jiang1(sunindextype *rowvals);
extern void dxdotdw_model0_jiang1(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_model0_jiang1(sunindextype *colptrs);
extern void dxdotdw_rowvals_model0_jiang1(sunindextype *rowvals);
extern void dxdotdp_model0_jiang1(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_model0_jiang1(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_model0_jiang1(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_model0_jiang1(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_model0_jiang1(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_model0_jiang1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_model0_jiang1(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_model0_jiang1(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_model0_jiang1(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_model0_jiang1(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_model0_jiang1(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_model0_jiang1(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_model0_jiang1(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_model0_jiang1(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_model0_jiang1 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_model0_jiang1()
        : amici::Model_ODE(
              59,                                // nx_rdata
              59,                            // nxtrue_rdata
              59,                               // nx_solver
              59,                           // nxtrue_solver
              59,                                      // ny
              59,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              45,                                      // nw
              127,                                   // ndwdx
              13770,                                   // ndwdp
              162,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              388,                                     // nnz
              59,                                     // ubw
              59,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{0.011, 64.941, 10000.0, 10000.0, 6.3e-05, 0.0001, 0.0005, 3.8617e-07, 4e-06, 1e-05, 4.8e-06, 1.18e-05, 0.005267, 3.8617e-07, 31.44, 20.47, 0.00011, 0.0005, 3.8617e-07, 0.00036, 0.00064, 7.8e-05, 0.00023, 29.6, 105.0, 3.8617e-07, 177.0, 2.5e-05, 8.1e-05, 1.1e-06, 0.0001, 0.00074, 0.00072, 0.0006, 0.0003, 5e-05, 2.5e-05, 0.00022, 3.8617e-07, 100.0, 100.0, 6.7e-06, 5e-06, 0.07, 3e-05, 2e-05, 0.0004, 8.375, 0.0006, 0.00045, 7.5e-06, 0.0006, 0.00045, 3.5e-05, 5e-06, 9.9211e-05, 1.73, 69.3, 0.037, 5.6e-06, 4.1e-06, 1.5e-06, 3e-07, 6.9e-05, 3e-05, 3.8617e-07, 900.0, 800.0, 5e-06, 2.5e-05, 3.8617e-07, 0.04, 0.39, 0.0019, 0.0071, 0.0001, 1.1e-05, 0.00017, 0.0016, 0.00011, 7.2e-05, 1e-05, 0.017, 1.6e-09, 0.0015, 3.8617e-07, 0.15, 337.0, 0.69, 0.012, 0.0087, 0.0004, 0.032, 0.0004, 0.002, 3.8617e-07, 1000.0, 300.0, 6.2, 0.0083, 0.002, 0.004, 4e-05, 0.0001, 0.0009, 0.00033211, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 0.0028, 0.00018, 0.0032, 8e-05, 0.00033211, 229.0, 498.0, 407.9, 9.8e-08, 2.1e-08, 5.9e-05, 9.9e-06, 0.00026, 9.2e-06, 9.963e-09, 426.8, 622.1, 1.9e-06, 2.8e-06, 5.7e-06, 5.4e-06, 3e-06, 2.8e-05, 0.002325, 0.000935, 0.00011, 3.8617e-07, 18000000.0, 0.0043, 4.5e-06, 5e-06, 3.9e-05, 4.5e-06, 5e-06, 0.004833, 0.0033211, 0.047, 0.0045, 0.075, 3.8617e-07, 31.44, 20.47, 0.00011, 0.0005, 0.05, 1.0, 0.00033211, 1.0, 1.0, 1.0, 1.0, 4.83, 3.675, 0.00017, 0.0014, 0.0007, 0.0003, 3.8617e-07, 260000.0, 570000.0, 1400.0, 26.0, 214.0, 4650.0, 35000000.0, 34000000.0, 3.8617e-07, 1000.0, 300.0, 6.2, 0.0083, 0.002, 0.004, 4e-05, 0.0001, 0.0009, 0.00033211, 1.0, 1.0, 1.0, 1.0, 3.5, 5.6, 4.18e-05, 0.00033, 0.00044, 0.00013, 3.3211e-05, 101.0, 78.0, 0.66, 3e-07, 3.1e-07, 4.2e-09, 3.2e-07, 3.9e-07, 3.1e-07, 3.3211e-05, 0.3, 2.18, 8.99, 1.19e-05, 7.53e-05, 2.4e-07, 7.6e-05, 2.42e-05, 1.08e-06, 1.2e-07, 3.9e-05, 3.8617e-07, 20.0, 200.0, 0.00024, 0.00019, 0.0079, 0.00013, 0.0016, 0.00015, 9.0, 5.1e-05, 0.00024, 0.016, 0.00037, 0.00163, 0.00011, 0.001, 3.99e-08, 3.4e-05, 0.001, 0.0399, 34.0, 3.8617e-07, 0.011, 0.000125, 0.333, 0.001, 0.001, 0.01, 0.00033211, 0.012, 0.1667, 3.8617e-07, 4.6e-14, 9.4e-08, 1.3e-10, 6e-12, 0.00025, 2.9e-05, 3.7e-07, 0.066, 9e-08, 9.6e-07, 9.5e-08, 0.051, 0.00033211, 1.0, 1.0, 1.0, 1.0, 3.5, 5.6, 4.18e-05, 0.00033, 0.00044, 0.00013, 0.00033211, 0.005, 1.11667, 3.8617e-07, 0.01295, 130.5, 500.0, 1000.0, 0.0003, 0.0002, 0.01, 0.143, 1000.0, 0.00033211, 1.5e-07, 1e-08, 3.8617e-07, 856.0, 3.6e-05, 3.5e-05, 6e-05, 0.00018, 0.0003, 0.00055, 6.9e-07, 5.9e-07, 5e-05, 1.3e-05, 2.5e-05},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(59, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_model0_jiang1(*this);
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
        J_model0_jiang1(J, t, x, p, k, h, w, dwdx);
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
        JB_model0_jiang1(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_model0_jiang1(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_model0_jiang1( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_model0_jiang1(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_model0_jiang1(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_model0_jiang1( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_model0_jiang1(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_model0_jiang1(rowvals);
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
        Jy_model0_jiang1(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_model0_jiang1(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_model0_jiang1(dsigmaydp, t, p, k, ip);
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
        dJydy_model0_jiang1( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_model0_jiang1(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_model0_jiang1(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_model0_jiang1( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_model0_jiang1( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_model0_jiang1( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_model0_jiang1(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_model0_jiang1(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_model0_jiang1( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_model0_jiang1(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_model0_jiang1(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_model0_jiang1(sigmay, t, p, k);
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
        sx0_model0_jiang1(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_model0_jiang1(sx0, t, x0, p, k, ip);
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
        w_model0_jiang1( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_model0_jiang1(x0, t, p, k);
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
        x0_fixedParameters_model0_jiang1(x0, t, p, k);
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
        xdot_model0_jiang1(xdot, t, x, p, k, h, w);
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
        y_model0_jiang1(y, t, x, p, k, h, w);
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
        x_solver_model0_jiang1( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_model0_jiang1( total_cl,  x_rdata);
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
"CS",
"",
"",
"",
"",
"",
"ACO",
"",
"",
"",
"",
"IDHa",
"",
"",
"",
"",
"",
"",
"OGDC",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"SCS",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"SDH",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"FM",
"",
"",
"",
"",
"MDH",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"AlaTA",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"AspTA",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"AGC",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"Complex_I",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"Complex_III",
"",
"",
"",
"",
"",
"",
"",
"",
"Complex_IV",
"",
"",
"CS",
"",
"",
"",
"",
"",
"",
"",
"",
"Complex_V",
"",
"",
"",
"ACO",
"",
"",
"",
"",
"",
"",
"OGC",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"MDH",
"",
"",
"",
"",
"",
"",
"",
"",
"AspTA",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"CIC",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"PC",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"Glycerol-3-phosphate dehydrogenase",
"",
"",
"Glycerol-3-phosphate dehydrogenase",
"",
"",
"MDH",
"",
"",
"",
"",
"",
"",
"AAC",
"",
"",
"IDHc",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"CIC",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"AAC",
"",
"",
"MDH",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"PYC",
"",
"",
"PDC",
"",
"",
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
        return std::vector<std::string>{"adenine diphosphate",
"adenine diphosphate",
"adenine monophosphate",
"adenine triphosphate",
"adenine triphosphate",
"acetyl CoA",
"acetyl CoA",
"slanine",
"aspartate",
"aspartate",
"carbon dioxide",
"citrate",
"citrate",
"coenzyme A",
"coenzyme A",
"ferricytochrome c",
"ferrocytochrome c",
"dihydrohxyacetone-phosphate",
"1,2-bisphospho-D-glycerate",
"electron transfer flavoprotein (oxidised form)",
"electron transfer flavoprotein (reduced form)",
"fructose-6-phosphate",
"FAD",
"FADH2",
"fructose-1,6-bisphosphate",
"fumarate",
"glycerol-3-phosphate",
"glyceraldehyde 3-phosphate",
"guanosine diphosphate",
"glucose",
"guanosine triphosphate",
"glutamate",
"glutamate",
"water",
"isocitrate",
"isocitrate",
"lactate",
"malate",
"malate",
"NAD",
"NADH",
"NADH",
"NADPH",
"NADPH",
"NADP",
"NADP+",
"NAD+",
"oxoglutarate",
"oxoglutarate",
"oxaloacetate",
"oxaloacetate",
"phosphoenolpyruvate",
"pyruvate",
"phosphate",
"pyruvate",
"ubiquinone",
"ubiquinol",
"succinyl-CoA",
"succinate",};
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
"x58",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"flow",
"GLCflow_Glc_F",
"hidden_1_k9b",
"hidden_1_k9f",
"v1_K1ATP",
"v1_K1GLC",
"v1_V1",
"v10_v10_CS",
"v10_Kib",
"v10_Kia",
"v10_Kb",
"v10_Ka",
"v10_V",
"v11_v11_ACO",
"v11_KcR",
"v11_KcF",
"v11_Kp",
"v11_Ks",
"v12_v12_IDHa",
"v12_f",
"v12_e",
"v12_d",
"v12_c",
"v12_b",
"v12_KcF",
"v14_v14_OGDC",
"v14_KcF",
"v14_Kir",
"v14_Kiq",
"v14_Kip",
"v14_Kic",
"v14_Kib",
"v14_Kia",
"v14_KmR",
"v14_KmP",
"v14_KmC",
"v14_KmB",
"v14_KmA",
"v15_v15_SCS",
"v15_Kc2",
"v15_Kc1",
"v15_Kir",
"v15_Kiq",
"v15_Kip",
"v15_Kic",
"v15_Kib",
"v15_Kia",
"v15_Keq",
"v15_KmP2",
"v15_KmC2",
"v15_KmQ",
"v15_KmP",
"v15_KmC",
"v15_KmB",
"v15_KmA",
"v16_v16_SDH",
"v16_KcR",
"v16_KcF",
"v16_Keq",
"v16_KiP2",
"v16_KiS1",
"v16_KmP2",
"v16_KmP1",
"v16_KmS2",
"v16_KmS1",
"v17_v17_FM",
"v17_KcR",
"v17_KcF",
"v17_Ks",
"v17_Kp",
"v18_v18_MDH",
"v18_KcR",
"v18_KcF",
"v18_KiP2",
"v18_KiP1",
"v18_KiS2",
"v18_KiS1",
"v18_KmP2",
"v18_KmP1",
"v18_KmS2",
"v18_KmS1",
"v2_K2ATP",
"v2_k2",
"v2_K2",
"v2_V2",
"v20_v20_AlaTA",
"v20_KcR",
"v20_KcF",
"v20_Keq",
"v20_KiP2",
"v20_KiS1",
"v20_KmP2",
"v20_KmP1",
"v20_KmS2",
"v20_KmS1",
"v21_v21_AspTA",
"v21_KcR",
"v21_KcF",
"v21_Keq",
"v21_KiP2",
"v21_KiS1",
"v21_KmP2",
"v21_KmP1",
"v21_KmS2",
"v21_KmS1",
"v22_v22_AGC",
"v22_delta",
"v22_gamma",
"v22_beta",
"v22_alpha",
"v22_KcR",
"v22_KcF",
"v22_KiP2",
"v22_KiP1",
"v22_KiS2",
"v22_KiS1",
"v24_v24_Complex_I",
"v24_KcR",
"v24_KcF",
"v24_Keq",
"v24_KiP2",
"v24_KiS1",
"v24_KmP2",
"v24_KmP1",
"v24_KmS2",
"v24_KmS1",
"v25_v25_Complex_III",
"v25_KcF",
"v25_k8",
"v25_Kq2",
"v25_Kq1",
"v25_Kb2",
"v25_Kb1",
"v25_KmB",
"v25_KmA",
"v26_v26_Complex_IV",
"v26_KcF",
"v26_Ks",
"v27_v10_CS",
"v27_Keq",
"v27_Kid",
"v27_Kib",
"v27_Kia",
"v27_Kc",
"v27_Kb",
"v27_Ka",
"v27_V",
"v28_v28_Complex_V",
"v28_Ki",
"v28_Km",
"v28_V",
"v29_v29_ACO",
"v29_KcR",
"v29_KcF",
"v29_Kp",
"v29_Ks",
"v3_k3b",
"v3_k3f",
"v30_v30_OGC",
"v30_delta",
"v30_gamma",
"v30_beta",
"v30_alpha",
"v30_KcR",
"v30_KcF",
"v30_KiP2",
"v30_KiP1",
"v30_KiS2",
"v30_KiS1",
"v31_v31_MDH",
"v31_kminus4",
"v31_kminus3",
"v31_kminus2",
"v31_kminus1",
"v31_k4",
"v31_k3",
"v31_k2",
"v31_k1",
"v32_v32_AspTA",
"v32_KcR",
"v32_KcF",
"v32_Keq",
"v32_KiP2",
"v32_KiS1",
"v32_KmP2",
"v32_KmP1",
"v32_KmS2",
"v32_KmS1",
"v33_v33_CIC",
"v33_delta",
"v33_gamma",
"v33_beta",
"v33_alpha",
"v33_KcR",
"v33_KcF",
"v33_KiP2",
"v33_KiP1",
"v33_KiS2",
"v33_KiS1",
"v34_v34_ETF_QO",
"v34_KcR",
"v34_KcF",
"v34_Keq",
"v34_KiP2",
"v34_KiS1",
"v34_KmP2",
"v34_KmP1",
"v34_KmS2",
"v34_KmS1",
"v35_v35_ACD",
"v35_KcR",
"v35_KcF",
"v35_Keq",
"v35_KiP2",
"v35_KiP1",
"v35_KiS2",
"v35_KiS1",
"v35_KmP2",
"v35_KmP1",
"v35_KmS2",
"v35_KmS1",
"v36_v36_PC",
"v36_KcR",
"v36_KcF",
"v36_Kir",
"v36_Kiq",
"v36_Kip",
"v36_Kic",
"v36_Kib",
"v36_Kia",
"v36_Keq",
"v36_KmR",
"v36_KmQ",
"v36_KmP",
"v36_KmC",
"v36_KmB",
"v36_KmA",
"v37_v37_GUT2P",
"v37_V",
"v37_K",
"v38_v38_GUT2P",
"v38_V",
"v38_K",
"v39_v39_MDH",
"v39_Knadp",
"v39_Kmal",
"v39_Kcat",
"v4_K4NAD",
"v4_K4GAP",
"v4_V4",
"v40_v40_AAC",
"v40_K",
"v40_V",
"v41_v41_IDHc",
"v41_phir123",
"v41_phir23",
"v41_phir13",
"v41_phir12",
"v41_phir3",
"v41_phir2",
"v41_phir1",
"v41_phir0",
"v41_phi12",
"v41_phi2",
"v41_phi1",
"v41_phi0",
"v42_v42_CIC",
"v42_delta",
"v42_gamma",
"v42_beta",
"v42_alpha",
"v42_KcR",
"v42_KcF",
"v42_KiP2",
"v42_KiP1",
"v42_KiS2",
"v42_KiS1",
"v43_v43_AAC",
"v43_K",
"v43_V",
"v44_v44_MDH",
"v44_Km",
"v44_Kcat",
"v5_k5b",
"v5_k5f",
"v6_K6ADP",
"v6_K6PEP",
"v6_V6",
"v7_k8b",
"v7_k8f",
"v8_v8_PYC",
"v8_K",
"v8_V",
"v9_v9_PDC",
"v9_KcF",
"v9_Kir",
"v9_Kiq",
"v9_Kip",
"v9_Kic",
"v9_Kib",
"v9_Kia",
"v9_KmR",
"v9_KmP",
"v9_KmC",
"v9_KmB",
"v9_KmA",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"ADP",
"ADP_cyt",
"AMP",
"ATP",
"ATP_cyt",
"Acetyl_CoA",
"Acetyl_CoA_cyt",
"Ala",
"Asp",
"Asp_cyt",
"CO2",
"Cit",
"Cit_cyt",
"CoA",
"CoA_cyt",
"Cytc2p",
"Cytc3p",
"DHAP",
"DPG",
"ETFox",
"ETFred",
"F6P",
"FAD",
"FADH2",
"FBP",
"Fum",
"G3P",
"GAP",
"GDP",
"GLC",
"GTP",
"Glu",
"Glu_cyt",
"H2O",
"IsoCit",
"IsoCitcyt",
"LAC",
"Mal",
"Mal_cyt",
"NAD",
"NADH",
"NADH_cyt",
"NADPH",
"NADPH_cyt",
"NADP_cyt",
"NADP_p",
"NAD_p",
"OG",
"OG_cyt",
"OXA",
"OXA_cyt",
"PEP",
"PYR_cyt",
"Pi",
"Pyr",
"Q",
"QH2",
"SCoA",
"Suc",};
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
"y58",};
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
