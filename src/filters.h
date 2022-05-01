#ifndef _FILTERS_H
#define _FILTERS_H
#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <RcppDist.h>
Rcpp::List vets_cpp_llh_equicor(Rcpp::NumericVector , arma::mat , arma::sp_mat , arma::sp_mat , arma::sp_mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat, arma::vec);
Rcpp::List vets_cpp_llh_diagonal(Rcpp::NumericVector , arma::mat , arma::sp_mat , arma::sp_mat , arma::sp_mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat, arma::vec);
Rcpp::List vets_cpp_llh_full(Rcpp::NumericVector , arma::mat , arma::sp_mat , arma::sp_mat , arma::sp_mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat, arma::vec);
Rcpp::List vets_cpp_llh_shrink(Rcpp::NumericVector , arma::mat , arma::sp_mat , arma::sp_mat , arma::sp_mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat, arma::vec);
Rcpp::List vets_cpp_filter(Rcpp::NumericVector , arma::mat , arma::sp_mat , arma::sp_mat , arma::sp_mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat);
Rcpp::List vets_cpp_predict(Rcpp::NumericVector , arma::mat , arma::mat , arma::sp_mat , arma::sp_mat , arma::sp_mat , arma::vec , arma::mat , arma::mat );
Rcpp::List vets_cpp_simulate(Rcpp::NumericVector , arma::mat , arma::mat , arma::sp_mat , arma::sp_mat , arma::sp_mat , arma::vec , arma::mat , arma::mat );

// constants

#ifndef VETS_LN_2PI
#define VETS_LN_2PI 1.837877066409345483560659472811235279722794947275566825634
#endif

#ifndef VETS_LARGE_POS_NUM
#define VETS_LARGE_POS_NUM 1.0E8
#endif

#ifndef VETS_MAX_EIGVAL_TOL
#define VETS_MAX_EIGVAL_TOL 1.01
#endif

// inline/template code

inline
arma::mat
shrinkcov(const arma::mat& E, const int n, const double rho)
{
    arma::mat V = arma::cov(E);
    arma::mat S = (rho * arma::trace(V)/static_cast<double>(n)) + (1.0 - rho) * V;

    return S;
}

inline
arma::mat
equicorrelation(const arma::rowvec& V, const int n, const double rho)
{
    arma::mat R(n, n);
    R.fill(rho);
    R.diag().ones();

    arma::mat S = arma::diagmat(V);

    return(S * R * S.t());
}

inline
double
power_iterations(const arma::mat& inp_mat, const double err_tol = 1.0e-04, const int max_iter = 1000)
{
    // arma::vec b = arma::randu(inp_mat.n_cols) + 0.2;
    arma::vec b = arma::ones(inp_mat.n_cols) - 0.02;

    int iter = 0;
    double err = 2*err_tol;

    double spec_rad = arma::dot(b,inp_mat*b) / arma::dot(b,b);
    double spec_rad_old = spec_rad;

    while (err > err_tol && iter < max_iter) {
        iter++;

        b = inp_mat * b;
        b /= arma::norm(b, 2);

        if ((iter > 100) && (iter % 20 == 0)) {
            spec_rad = arma::dot(b,inp_mat*b) / arma::dot(b,b);

            err = std::abs(spec_rad - spec_rad_old);
            spec_rad_old = spec_rad;
        }
    }

    if (iter >= max_iter) {
        spec_rad = arma::dot(b,inp_mat*b) / arma::dot(b,b);
    }

    return spec_rad;
}

inline
double
power_iterations_fast(const arma::mat& inp_mat)
{
    // arma::vec b = arma::randu(inp_mat.n_cols) + 0.2;
    arma::vec b = arma::ones(inp_mat.n_cols) - 0.02;

    for (int i = 0; i < 200; ++i) {
        b = inp_mat * b;
        b /= arma::norm(b, 2);
    }

    double spec_rad = arma::dot(b,inp_mat*b) / arma::dot(b,b);

    return spec_rad;
}

#endif
