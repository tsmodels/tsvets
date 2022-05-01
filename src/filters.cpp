#include "filters.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::export]]
Rcpp::List vets_cpp_llh_equicor(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, 
                                arma::mat Y, arma::mat X, arma::mat beta, arma::mat good, arma::vec select)
{
  // model: time[t] series[n] xreg[0,1] rho (correlation)
    try {
        int t = static_cast<int>(model[0]);
        int n = static_cast<int>(model[1]);
        int use_x = static_cast<int>(model[2]);
        double rho = static_cast<double>(model[3]);

        arma::mat GA_mat = Gmat * Amat;

        arma::mat Error = arma::zeros(Y.n_rows, Y.n_cols-1);
        arma::mat Aux = arma::zeros(Y.n_rows, select.n_elem);
        // check stability condition

        arma::mat Cond = Fmat - GA_mat * Hmat;

        // double spec_radius = power_iterations(Cond);
        double spec_radius = power_iterations_fast(Cond);

        bool stability_test = (spec_radius < VETS_MAX_EIGVAL_TOL);

        if (!stability_test) {
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = VETS_LARGE_POS_NUM,
                                                   Rcpp::Named("condition") = 1);
            return(output);
        } else {
            int kselect = 0;
            for (int i = 1; i < t; i++) {
                arma::vec Yhat = Hmat * States.col(i-1);

                if (use_x == 1) {
                    Yhat += beta * X.col(i);
                }
                Error.col(i-1) = Y.col(i) - Yhat;
                // zero out the errors of missing data
                Error.col(i-1) = Error.col(i-1) % good.col(i-1);
                if (select(i - 1) == 1) {
                  // select non missing data errors
                  Aux.col(kselect) =  Error.col(i-1);
                  kselect+=1;
                }
                States.col(i) = Fmat * States.col(i-1) + GA_mat * Error.col(i-1);
            }
            Aux = Aux.t();
            Error = Error.t();
            // select complete non missingness matrix
            arma::rowvec Sig = arma::stddev(Aux, 0, 0);
            arma::mat S = equicorrelation(Sig, n, rho);

            arma::vec eigvalS;
            arma::mat eigvecS;
            arma::eig_sym(eigvalS, eigvecS, S);

            arma::mat E = arma::pow(Aux * eigvecS, 2);
            double sL = arma::accu( E * (1.0/eigvalS) );

            double ldet = arma::accu(arma::log(eigvalS));

            double lconstant = - 0.5 * t * (n*VETS_LN_2PI + ldet);

            // Negative log likelihood
            double loglik = - (lconstant - 0.5 * sL);

            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                                                   Rcpp::Named("condition") = 0);
            return(output);
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsvets--> vets filter exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List vets_cpp_llh_diagonal(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, 
                                 arma::mat X, arma::mat beta, arma::mat good, arma::vec select)
{
    // model: time[t] series[n] xreg[0,1]
    try {
        int t =  static_cast<int>(model[0]);
        int n =  static_cast<int>(model[1]);
        int use_x = static_cast<int>(model[2]);

        arma::mat GA_mat = Gmat * Amat;

        arma::mat Error = arma::zeros(Y.n_rows, Y.n_cols-1);
        arma::mat Aux = arma::zeros(Y.n_rows, select.n_elem);
        // check stability condition

        arma::mat Cond = Fmat - GA_mat * Hmat;

        // double spec_radius = power_iterations(Cond);
        double spec_radius = power_iterations_fast(Cond);
        // std::cout << "spec_radius = " << spec_radius << "\n";

        bool stability_test = (spec_radius < VETS_MAX_EIGVAL_TOL);

        if (!stability_test) {
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = VETS_LARGE_POS_NUM,
                                                   Rcpp::Named("condition") = 1);
            return(output);
        } else {
            int kselect = 0;
            for (int i = 1; i < t; i++) {
                arma::vec Yhat = Hmat * States.col(i-1);
                
                if (use_x == 1) {
                    Yhat += beta * X.col(i);
                }

                Error.col(i-1) = Y.col(i) - Yhat;
                // zero out the errors of missing data
                Error.col(i-1) = Error.col(i-1) % good.col(i-1);
                if (select(i - 1) == 1) {
                  // select non missing data errors
                  Aux.col(kselect) =  Error.col(i-1);
                  kselect+=1;
                }
                States.col(i) = Fmat * States.col(i-1) + GA_mat * Error.col(i-1);
            }
            Aux = Aux.t();
            Error = Error.t();

            arma::rowvec V = arma::var(Aux, 0, 0);
            arma::mat S = arma::diagmat(V);
            double ldet = arma::accu(arma::log(V));
            
            double lconstant = -0.5 * t * (n*VETS_LN_2PI + ldet);

            arma::vec eigvalS;
            arma::mat eigvecS;
            arma::eig_sym(eigvalS, eigvecS, S);
            
            arma::mat E = arma::pow(Aux * eigvecS, 2);
            double sL = arma::accu( E * (1.0/eigvalS) );
            
            // Negative log likelihood
            double loglik = -1.0 * (lconstant - 0.5 * sL);
            
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                                                Rcpp::Named("condition") = 0);
            return(output);
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsvets--> vets filter exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List vets_cpp_llh_full(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, arma::mat X, 
                             arma::mat beta, arma::mat good, arma::vec select)
{
  // model: time[t] series[n]
    try {
        int t = static_cast<int>(model[0]);
        int n = static_cast<int>(model[1]);
        int use_x = static_cast<int>(model[2]);

        arma::mat GA_mat = Gmat * Amat;

        arma::mat Error = arma::zeros(Y.n_rows, Y.n_cols-1);
        arma::mat Aux = arma::zeros(Y.n_rows, select.n_elem);
        // check stability condition

        arma::mat Cond = Fmat - GA_mat * Hmat;

        // double spec_radius = power_iterations(Cond);
        double spec_radius = power_iterations_fast(Cond);

        bool stability_test = (spec_radius < VETS_MAX_EIGVAL_TOL);

        if (!stability_test) {
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = VETS_LARGE_POS_NUM,
                                                   Rcpp::Named("condition") = 1);
            return(output);
        } else {
            int kselect = 0;
            for (int i = 1; i < t; i++) {
                arma::vec Yhat = Hmat * States.col(i-1);

                if (use_x == 1) {
                    Yhat += beta * X.col(i);
                }
                    
                Error.col(i-1) = Y.col(i) - Yhat;
                // zero out the errors of missing data
                Error.col(i-1) = Error.col(i-1) % good.col(i-1);
                if (select(i - 1) == 1) {
                  // select non missing data errors
                  Aux.col(kselect) =  Error.col(i-1);
                  kselect+=1;
                }
                States.col(i) = Fmat * States.col(i-1) + GA_mat * Error.col(i-1);
            }

            Error = Error.t();
            Aux = Aux.t();
            // needs fixing
            arma::mat S = arma::cov(Aux, 0);
            
            arma::vec eigvalS;
            arma::mat eigvecS;
            arma::eig_sym(eigvalS, eigvecS, S);

            arma::mat E = arma::pow(Aux * eigvecS, 2);
            double sL = arma::accu( E * (1.0/eigvalS) );

            double ldet = arma::accu(arma::log(eigvalS));

            double lconstant = - 0.5 * t * (n*VETS_LN_2PI + ldet);

            // Negative log likelihood
            double loglik = - (lconstant - 0.5 * sL);

            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                                                Rcpp::Named("condition") = 0);
            return(output);
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsvets--> vets filter exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List vets_cpp_llh_shrink(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, 
                               arma::mat X, arma::mat beta, arma::mat good, arma::vec select)
{
  // model: time[t] series[n]
    try {
        int t = static_cast<int>(model[0]);
        int n = static_cast<int>(model[1]);
        int use_x = static_cast<int>(model[2]);
        double rho = static_cast<double>(model[3]);

        arma::mat GA_mat = Gmat * Amat;

        arma::mat Error = arma::zeros(Y.n_rows, Y.n_cols-1);
        arma::mat Aux = arma::zeros(Y.n_rows, select.n_elem);
        // check stability condition

        arma::mat Cond = Fmat - GA_mat * Hmat;

        // double spec_radius = power_iterations(Cond);
        double spec_radius = power_iterations_fast(Cond);

        bool stability_test = (spec_radius < VETS_MAX_EIGVAL_TOL);

        if (!stability_test) {
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = VETS_LARGE_POS_NUM,
                                                   Rcpp::Named("condition") = 1);
            return(output);
        } else {
            int kselect = 0;
            for (int i = 1; i < t; i++) {
                arma::vec Yhat = Hmat * States.col(i-1);

                if (use_x == 1) {
                    Yhat += beta * X.col(i);
                }
                Error.col(i-1) = Y.col(i) - Yhat;
                // zero out the errors of missing data
                Error.col(i-1) = Error.col(i-1) % good.col(i-1);
                if (select(i - 1) == 1) {
                  // select non missing data errors
                  Aux.col(kselect) =  Error.col(i-1);
                  kselect+=1;
                }
                States.col(i) = Fmat * States.col(i-1) + GA_mat * Error.col(i-1);
            }

            Error = Error.t();
            Aux = Aux.t();
            // needs fixing
            arma::mat S = shrinkcov(Aux, n, rho);
            
            arma::vec eigvalS;
            arma::mat eigvecS;
            arma::eig_sym(eigvalS, eigvecS, S);

            arma::mat E = arma::pow(Aux * eigvecS, 2);
            double sL = arma::accu( E * (1.0/eigvalS) );

            double ldet = arma::accu(arma::log(eigvalS));

            double lconstant = - 0.5 * t * (n*VETS_LN_2PI + ldet);

            // Negative log likelihood
            double loglik = - (lconstant - 0.5 * sL);

            Rcpp::List output = Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                                                Rcpp::Named("condition") = 0);
            return(output);
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsvets--> vets filter exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List vets_cpp_filter(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, arma::mat X, arma::mat beta, 
                           arma::mat good)
{
  // model: time[t] series[n] xreg[0,1] rho (correlation)
  try {
    int t = static_cast<int>(model[0]);
    int use_x = static_cast<int>(model[2]);
    arma::mat Yhat(Y.n_rows, Y.n_cols);
    Yhat.fill(0);
    arma::mat Error(Y.n_rows, Y.n_cols);
    Error.fill(0);
    arma::mat Cond = Fmat - Gmat * Amat * Hmat;
    int i;
    for(i=1;i<t;i++) {
      Yhat.col(i) = Hmat * States.col(i-1);
      if (use_x == 1) {
        Yhat.col(i)+= beta * X.col(i);
      }
      Error.col(i) = Y.col(i) - Yhat.col(i);
      // zero out the errors of missing data
      Error.col(i) = Error.col(i) % good.col(i-1);
      States.col(i) = Fmat * States.col(i-1) +  Gmat * Amat * Error.col(i);
    }
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("fitted") = Yhat.t(),
                                           Rcpp::Named("States") = States.t(),
                                           Rcpp::Named("Error") = Error.t());
    return(output);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "tsvets--> vets filter exception (unknown reason)" );
  }
  return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List vets_cpp_predict(Rcpp::NumericVector model, arma::mat S, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::vec istate, arma::mat X, arma::mat beta)
{
  try {
    int h = static_cast<int>(model[0]);
    int nsim = static_cast<int>(model[1]);
    int n = static_cast<int>(model[2]);
    int use_x = static_cast<int>(model[3]);
    int j, i;
    arma::vec mu(n);
    mu.fill(0);
    int n_states = istate.n_rows;
    arma::cube simY(n, h, nsim);
    arma::cube simStates(n_states, h, nsim);
    arma::mat Y = arma::mat(n, h + 1);
    arma::cube E(n, h, nsim);
    arma::mat States = arma::mat(n_states, h + 1);
    for(j=0;j<nsim;j++){
      arma::mat R = rmvnorm(h + 1, mu, S);
      R = R.t();
      E.slice(j) = R.cols(1,h);
      Y.fill(0);
      States.fill(0);
      States.col(0) = istate;
      for(i=1;i<=h;i++) {
        Y.col(i) = Hmat * States.col(i-1) + R.col(i);
        if (use_x == 1) {
          Y.col(i)+= beta * X.col(i);
        }
        States.col(i) = Fmat * States.col(i-1) +  Gmat * Amat * R.col(i);
      }
      simY.slice(j) = Y.cols(1, h);
      simStates.slice(j) = States.cols(1, h);
      
    }
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("Y") = simY,
                                           Rcpp::Named("States") = simStates,
                                           Rcpp::Named("Error") = E);
    return(output);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "tsvets--> vets predict exception (unknown reason)" );
  }
  return R_NilValue;
}


// [[Rcpp::export]]
Rcpp::List vets_cpp_simulate(Rcpp::NumericVector model, arma::cube E, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::vec istate, arma::mat X, arma::mat beta)
{
  try {
    int h = static_cast<int>(model[0]);
    int nsim = static_cast<int>(model[1]);
    int n = static_cast<int>(model[2]);
    int use_x = static_cast<int>(model[3]);
    int j, i;
    arma::vec mu(n);
    mu.fill(0);
    int n_states = istate.n_rows;
    arma::cube simY(n, h, nsim);
    arma::cube simStates(n_states, h, nsim);
    arma::mat Y = arma::mat(n, h + 1);
    arma::mat States = arma::mat(n_states, h + 1);
    for(j=0;j<nsim;j++){
      Y.fill(0);
      States.fill(0);
      States.col(0) = istate;
      arma::mat R = E.slice(j);
      for(i=1;i<=h;i++) {
        Y.col(i) = Hmat * States.col(i-1) + R.col(i);
        if (use_x == 1) {
          Y.col(i)+= beta * X.col(i);
        }
        States.col(i) = Fmat * States.col(i-1) +  Gmat * Amat * R.col(i);
      }
      simY.slice(j) = Y.cols(1, h);
      simStates.slice(j) = States.cols(1, h);
    }
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("Y") = simY,
                                           Rcpp::Named("States") = simStates);
    return(output);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "tsvets--> vets simulation exception (unknown reason)" );
  }
  return R_NilValue;
}
