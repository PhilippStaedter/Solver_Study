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
extern void J_Sasagawa2005(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Sasagawa2005(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Sasagawa2005(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Sasagawa2005(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Sasagawa2005(sunindextype *colptrs);
extern void JSparse_rowvals_Sasagawa2005(sunindextype *rowvals);
extern void JSparseB_Sasagawa2005(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Sasagawa2005(sunindextype *colptrs);
extern void JSparseB_rowvals_Sasagawa2005(sunindextype *rowvals);
extern void Jy_Sasagawa2005(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Sasagawa2005(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Sasagawa2005(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Sasagawa2005(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Sasagawa2005(sunindextype *rowvals, int index);
extern void dwdp_Sasagawa2005(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Sasagawa2005(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Sasagawa2005(sunindextype *colptrs);
extern void dwdx_rowvals_Sasagawa2005(sunindextype *rowvals);
extern void dxdotdw_Sasagawa2005(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Sasagawa2005(sunindextype *colptrs);
extern void dxdotdw_rowvals_Sasagawa2005(sunindextype *rowvals);
extern void dxdotdp_Sasagawa2005(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Sasagawa2005(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Sasagawa2005(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Sasagawa2005(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Sasagawa2005(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Sasagawa2005(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Sasagawa2005(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Sasagawa2005(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Sasagawa2005(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Sasagawa2005(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Sasagawa2005(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Sasagawa2005(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Sasagawa2005(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Sasagawa2005(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Sasagawa2005 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Sasagawa2005()
        : amici::Model_ODE(
              99,                                // nx_rdata
              99,                            // nxtrue_rdata
              99,                               // nx_solver
              99,                           // nxtrue_solver
              99,                                      // ny
              99,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              150,                                      // nw
              318,                                   // ndwdx
              35100,                                   // ndwdp
              392,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              623,                                     // nnz
              99,                                     // ubw
              99,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{0.0001, 0.0001, 0.0029666, 2.2833, 0.02, 10.0, 0.0168, 0.03, 0.0168, 0.03, 0.001, 4.0, 0.2, 0.5, 0.2, 10.0, 1e-05, 0.002, 0.2, 0.5, 0.002, 0.002, 0.2, 10.0, 1.0, 0.05, 0.001, 0.2, 0.5, 0.05, 0.001, 0.05, 0.001, 1.0, 0.2, 10.0, 0.2, 10.0, 0.2, 10.0, 0.2, 0.5, 0.2, 10.0, 0.05, 0.001, 0.005, 0.1, 0.2, 0.005, 0.005, 0.002, 1.0, 0.2, 1.0, 0.2, 1.0, 1.0, 0.2, 1.0, 0.2, 0.5, 0.2, 0.5, 0.05, 0.05, 1.0, 0.2, 1.0, 0.05, 0.001, 0.001, 0.01, 0.12, 25.641, 1.0, 25.641, 1.0, 0.5, 60.0, 0.5, 60.0, 0.5, 60.0, 15.657, 3.0, 15.657, 3.0, 15.657, 3.0, 15.657, 3.0, 0.075, 10.0, 0.0001667, 0.0001166, 0.01, 0.024, 0.02, 2.0, 6.4e-05, 6.2, 1.0, 0.00063, 0.00042, 0.0022, 0.2, 10.0, 0.2, 10.0, 0.1, 5.0, 0.1, 5.0, 0.2, 10.0, 0.2, 10.0, 0.1, 0.1, 2.0, 0.1, 5.0, 0.1, 5.0, 2.0, 0.0022, 0.0022, 0.0022, 0.0022, 0.00042, 0.00042, 0.00042, 0.2, 10.0, 0.2, 10.0, 0.00063, 0.00063, 0.00063, 0.00063, 0.00063, 0.00063, 0.2, 1.0, 0.2, 1.0, 0.2, 10.0, 0.2, 10.0, 0.0022, 0.00042, 0.0022, 0.00042, 0.00042, 0.00027778, 0.0008333, 0.2, 10.0, 0.2, 10.0, 0.2, 10.0, 0.2, 0.5, 0.2, 1.0, 0.2, 1.0, 1.0, 10.0, 1.0, 2.0, 0.1, 0.02, 25.641, 1.0, 0.6, 16.304, 0.6, 16.304, 0.6, 16.304, 0.15, 1.0, 10.0, 1.0, 10.0, 1.0, 2.0, 2.0, 15.625, 2.0, 15.625, 2.0, 15.625, 2.0, 15.625, 0.8, 6.25, 0.8, 6.25, 0.8, 6.25, 0.8, 6.25, 1.2, 9.375, 1.2, 9.375, 1.2, 9.375, 1.2, 9.375, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.001, 0.24, 15.0, 0.24, 15.0, 0.06, 0.06},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(99, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Sasagawa2005(*this);
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
        J_Sasagawa2005(J, t, x, p, k, h, w, dwdx);
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
        JB_Sasagawa2005(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_Sasagawa2005(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Sasagawa2005( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Sasagawa2005(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Sasagawa2005(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Sasagawa2005( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Sasagawa2005(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Sasagawa2005(rowvals);
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
        Jy_Sasagawa2005(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_Sasagawa2005(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_Sasagawa2005(dsigmaydp, t, p, k, ip);
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
        dJydy_Sasagawa2005( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Sasagawa2005(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Sasagawa2005(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Sasagawa2005( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Sasagawa2005( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Sasagawa2005( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Sasagawa2005(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Sasagawa2005(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Sasagawa2005( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_Sasagawa2005(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_Sasagawa2005(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_Sasagawa2005(sigmay, t, p, k);
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
        sx0_Sasagawa2005(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_Sasagawa2005(sx0, t, x0, p, k, ip);
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
        w_Sasagawa2005( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Sasagawa2005(x0, t, p, k);
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
        x0_fixedParameters_Sasagawa2005(x0, t, p, k);
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
        xdot_Sasagawa2005(xdot, t, x, p, k, h, w);
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
        y_Sasagawa2005(y, t, x, p, k, h, w);
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
        x_solver_Sasagawa2005( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Sasagawa2005( total_cl,  x_rdata);
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
"",};
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
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
"proteasome",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
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
"x98",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"re1_re1_k2",
"re1_re1_k1",
"re2_re2_k2",
"re2_re2_k1",
"re8_re8_k2",
"re8_re8_k1",
"J3_J3_k2",
"J3_J3_k1",
"J4_J4_k2",
"J4_J4_k1",
"J5_J5_k2",
"J5_J5_k1",
"J6_J6_k2",
"J6_J6_k1",
"J7_J7_k2",
"J7_J7_k1",
"J8_J8_k2",
"J8_J8_k1",
"J9_J9_k2",
"J9_J9_k1",
"J10_J10_k",
"J11_J11_k",
"J12_J12_k2",
"J12_J12_k1",
"J13_J13_k",
"J14_J14_k",
"J15_J15_k",
"J16_J16_k2",
"J16_J16_k1",
"J17_J17_k",
"J18_J18_k",
"J19_J19_k",
"J20_J20_k",
"J21_J21_k",
"J22_J22_k2",
"J22_J22_k1",
"J23_J23_k2",
"J23_J23_k1",
"J24_J24_k2",
"J24_J24_k1",
"J25_J25_k2",
"J25_J25_k1",
"J27_J27_k2",
"J27_J27_k1",
"J28_J28_k",
"J29_J29_k",
"J30_J30_k",
"J31_J31_Km1",
"J31_J31_Vmax",
"J32_J32_k",
"J33_J33_k",
"J34_J34_k2",
"J34_J34_k1",
"J35_J35_k2",
"J35_J35_k1",
"J36_J36_k2",
"J36_J36_k1",
"J37_J37_k",
"J38_J38_k2",
"J38_J38_k1",
"J39_J39_k2",
"J39_J39_k1",
"J40_J40_k2",
"J40_J40_k1",
"J41_J41_k",
"J42_J42_k",
"J43_J43_k",
"J44_J44_k2",
"J44_J44_k1",
"J45_J45_k",
"J46_J46_k",
"J47_J47_k",
"J49_J49_k2",
"J49_J49_k1",
"J50_J50_Km1",
"J50_J50_Vmax",
"J51_J51_Km1",
"J51_J51_Vmax",
"J52_J52_k2",
"J52_J52_k1",
"J53_J53_k2",
"J53_J53_k1",
"J54_J54_k2",
"J54_J54_k1",
"J57_J57_Km1",
"J57_J57_Vmax",
"J58_J58_Km1",
"J58_J58_Vmax",
"J61_J61_Km1",
"J61_J61_Vmax",
"J62_J62_Km1",
"J62_J62_Vmax",
"J63_J63_k2",
"J63_J63_k1",
"J66_J66_k",
"J67_J67_k",
"J68_J68_Km1",
"J68_J68_Vmax",
"J69_J69_Km1",
"J69_J69_Vmax",
"J70_J70_k2",
"J70_J70_k1",
"J71_J71_k",
"J72_J72_k",
"J73_J73_k",
"J74_J74_k",
"J75_J75_k2",
"J75_J75_k1",
"J76_J76_k2",
"J76_J76_k1",
"J77_J77_k2",
"J77_J77_k1",
"J78_J78_k2",
"J78_J78_k1",
"J79_J79_k2",
"J79_J79_k1",
"J80_J80_k2",
"J80_J80_k1",
"J81_J81_k",
"J82_J82_k",
"J83_J83_k",
"J84_J84_k2",
"J84_J84_k1",
"J85_J85_k2",
"J85_J85_k1",
"J86_J86_k",
"J87_J87_k",
"J88_J88_k",
"J89_J89_k",
"J90_J90_k",
"J92_J92_k",
"J93_J93_k",
"J94_J94_k",
"J95_J95_k2",
"J95_J95_k1",
"J96_J96_k2",
"J96_J96_k1",
"J97_J97_k",
"J98_J98_k",
"J99_J99_k",
"J100_J100_k",
"J101_J101_k",
"J102_J102_k",
"J103_J103_k2",
"J103_J103_k1",
"J104_J104_k2",
"J104_J104_k1",
"J105_J105_k2",
"J105_J105_k1",
"J106_J106_k2",
"J106_J106_k1",
"J107_J107_k",
"J108_J108_k",
"J109_J109_k",
"J110_J110_k",
"J112_J112_k",
"J113_J113_k2",
"J113_J113_k1",
"J115_J115_k2",
"J115_J115_k1",
"J116_J116_k2",
"J116_J116_k1",
"J117_J117_k2",
"J117_J117_k1",
"J119_J119_k2",
"J119_J119_k1",
"J118_J118_k2",
"J118_J118_k1",
"J120_J120_k2",
"J120_J120_k1",
"J121_J121_Km1",
"J121_J121_Vmax",
"J122_J122_Km1",
"J122_J122_Vmax",
"J123_J123_Km1",
"J123_J123_Vmax",
"J124_J124_Km1",
"J124_J124_Vmax",
"J133_J133_k2",
"J133_J133_k1",
"J134_J134_k2",
"J134_J134_k1",
"J135_J135_k2",
"J135_J135_k1",
"J136_J136_k",
"J137_J137_Km1",
"J137_J137_Vmax",
"J138_J138_Km1",
"J138_J138_Vmax",
"J139_J139_Km1",
"J139_J139_Vmax",
"J140_J140_k2",
"J140_J140_k1",
"J141_J141_k2",
"J141_J141_k1",
"J142_J142_k2",
"J142_J142_k1",
"J143_J143_k2",
"J143_J143_k1",
"J144_J144_k2",
"J144_J144_k1",
"J145_J145_k2",
"J145_J145_k1",
"J146_J146_k2",
"J146_J146_k1",
"J147_J147_k2",
"J147_J147_k1",
"J148_J148_k2",
"J148_J148_k1",
"J149_J149_k2",
"J149_J149_k1",
"J150_J150_k2",
"J150_J150_k1",
"J151_J151_k2",
"J151_J151_k1",
"J152_J152_k",
"J153_J153_k",
"J154_J154_k",
"J155_J155_k",
"J156_J156_k",
"J157_J157_k",
"J158_J158_k",
"J159_J159_k",
"J160_J160_k",
"J161_J161_k",
"J162_J162_k",
"J163_J163_k",
"J164_J164_k",
"J165_J165_k2",
"J165_J165_k1",
"J166_J166_k2",
"J166_J166_k1",
"J167_J167_k",
"J168_J168_k",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"EGFR",
"L_EGFR",
"L_EGFR_dimer",
"SOS",
"L_dpEGFR",
"pSOS",
"SOS_Grb2",
"Grb2",
"Dok",
"pDok",
"Crk",
"FRS2",
"Shc",
"pSOS_Grb2",
"Rap1_GDP",
"MEK",
"MKP3",
"pShc_dpEGFR",
"dpEGFR_c_Cbl",
"B_Raf_Rap1_GTP",
"pShc_dpEGFR_c_Cbl",
"pFRS2_dpEGFR_c_Cbl",
"Shc_dpEGFR",
"c_Cbl",
"RasGAP",
"c_Raf",
"B_Raf",
"ERK",
"PP2A",
"Ras_GDP",
"Rap1GAP",
"C3G",
"NGFR",
"pShc",
"pFRS2_dpEGFR",
"pTrkA_endo",
"MEK_ERK",
"pMEK_ERK",
"FRS2_dpEGFR_c_Cbl_ubiq",
"Crk_C3G_pFRS2_dpEGFR_c_Cbl",
"pShc_dpEGFR_c_Cbl_ubiq",
"Crk_C3G_pFRS2_dpEGFR",
"Grb2_SOS_pShc_dpEGFR_c_Cbl_ubiq",
"Grb2_SOS_pShc_dpEGFR_c_Cbl",
"Shc_dpEGFR_c_Cbl_ubiq",
"dpEGFR_c_Cbl_ubiq",
"proteosome",
"Grb2_SOS_pShc",
"Shc_dpEGFR_c_Cbl",
"Grb2_SOS_pShc_dpEGFR",
"pFRS2",
"FRS2_dpEGFR",
"pDok_RasGAP",
"pMEK",
"FRS2_dpEGFR_c_Cbl",
"pFRS2_dpEGFR_c_Cbl_ubiq",
"Ras_GTP",
"Crk_C3G_pFRS2_dpEGFR_c_Cbl_ubiq",
"c_Raf_Ras_GTP",
"B_Raf_Ras_GTP",
"ppMEK",
"ppERK",
"pTrkA",
"Crk_C3G",
"Rap1_GTP",
"L_NGFR",
"ppMEK_ERK",
"dppERK",
"Shc_pTrkA",
"Shc_pTrkA_endo",
"pShc_pTrkA",
"pFRS2_pTrkA",
"FRS2_pTrkA",
"pShc_pTrkA_endo",
"FRS2_pTrkA_endo",
"pFRS2_pTrkA_endo",
"Crk_C3G_pFRS2_pTrkA_endo",
"Grb2_SOS_pShc_pTrkA",
"Crk_C3G_pFRS2_pTrkA",
"Grb2_SOS_pShc_pTrkA_endo",
"c_Raf_Ras_GTP_MEK",
"c_Raf_Ras_GTP_pMEK",
"c_Raf_Ras_GTP_MEK_ERK",
"c_Raf_Ras_GTP_pMEK_ERK",
"B_Raf_Ras_GTP_MEK",
"B_Raf_Ras_GTP_pMEK",
"B_Raf_Ras_GTP_MEK_ERK",
"B_Raf_Ras_GTP_pMEK_ERK",
"B_Raf_Rap1_GTP_MEK",
"B_Raf_Rap1_GTP_pMEK",
"B_Raf_Rap1_GTP_MEK_ERK",
"B_Raf_Rap1_GTP_pMEK_ERK",
"ppERK_MKP3",
"dppERK_MKP3",
"pro_TrkA",
"NGF",
"EGF",
"pro_EGFR",
"degradation",};
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
"y98",};
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
