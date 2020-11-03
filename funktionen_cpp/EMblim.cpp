#include <RcppEigen.h>

#include <iostream>

#include <stdio.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map; 
using Eigen::MatrixXd; 
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using namespace Eigen;

// [[Rcpp::export]]
Rcpp::List emBLIMcpp(
  const Map < MatrixXd > R,
    const Map < MatrixXd > W,
      const Map < MatrixXd > K,
        const Map < VectorXd > NR,
          const Map < VectorXd > PKr,
            const Map < VectorXd > betar,
              const Map < VectorXd > etar,
                const int maxiter = 10000,
                  const double tol = 1e-07,
                    const bool fdb = true

) {

  //copy mapped input

  VectorXd PK = PKr;
  VectorXd eta = etar;
  VectorXd beta = betar;

  // declare objects including type      
  MatrixXd PRK;
  MatrixXd PR;
  MatrixXd PKR;

  MatrixXd diff(3, 1);

  VectorXd PKold;
  VectorXd etaold;
  VectorXd betaold;

  int iter = 0;
  double eps = 1e-9;
  bool converged = false;

  do {

    // E-Step
    //prevent log zero // maybe .array() <= eps is faster....
    for (int i = 0; i < eta.rows(); i++) {
      if (eta(i) < eps) {
        eta(i) += eps;
      } else if (eta(i) > (1 - eps)) {
        eta(i) -= eps;
      }

      if (beta(i) < eps) {
        beta(i) += eps;
      } else if (beta(i) > (1 - eps)) {
        beta(i) -= eps;
      }
    }

    PRK =
      (R.cwiseProduct(
          ((1 - beta.array()).log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        K.transpose() +

        R.cwiseProduct(
          (eta.array().log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        (
          ((1 - K.array()).matrix()).transpose()
        ) +

        (1 - R.array()).matrix().cwiseProduct(
          (beta.array().log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        K.transpose() +

        (1 - R.array()).matrix().cwiseProduct(
          ((1 - eta.array()).log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        (((1 - K.array()).matrix()).transpose())
      ).array().exp().matrix();

    PR = PRK * PK;

    PKR = (PK * PR.cwiseInverse().transpose()).cwiseProduct(PRK.transpose());

    // M - Step

    // save old estimates
    PKold = PK;
    etaold = eta;
    betaold = beta;

    // PK
    PK = (PKR * NR);
    PK /= NR.sum();

    // beta, eta
    beta = (
      (((PKR.transpose() * K).cwiseProduct(W)).transpose()) * NR
    ).cwiseQuotient(
      ((PKR.transpose() * K).transpose()) * NR);


    eta = (
      (((PKR.transpose() * ((1 - K.array()).matrix())).cwiseProduct(R)).transpose()) * NR
    ).cwiseQuotient(
      ((PKR.transpose() * (1 - K.array()).matrix()).transpose()) * NR);

    diff(0, 0) = ((PKold - PK).cwiseAbs().maxCoeff()) < tol;
    diff(1, 0) = ((etaold - eta).cwiseAbs().maxCoeff()) < tol;
    diff(2, 0) = ((betaold - beta).cwiseAbs().maxCoeff()) < tol;

    iter++;

    if (fdb) {
      if (iter % 50 == 0) {
        Rcpp::Rcout << ".  " << "Iteration #: " << iter << std::endl;
        Rcpp::checkUserInterrupt();
      } else {
        Rcpp::Rcout << ".";
      }
    }
  }
  while (diff.sum() < 3 && iter < maxiter);

  if (diff.sum() >= 3) {
    converged = true;
  }

  if (fdb) {
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "     **** DONE ****" << std::endl;

    if (iter >= maxiter) {
      Rcpp::Rcout << "Maximum # of " << maxiter << " Iterations reached." << std::endl;
    }

    if (converged == false) {
      Rcpp::Rcout << "EM-Algorithm did NOT converged!" << std::endl;
    }

  }



 // compute PKR with final estimates
   for (int i = 0; i < eta.rows(); i++) {
      if (eta(i) < eps) {
        eta(i) += eps;
      } else if (eta(i) > (1 - eps)) {
        eta(i) -= eps;
      }

      if (beta(i) < eps) {
        beta(i) += eps;
      } else if (beta(i) > (1 - eps)) {
        beta(i) -= eps;
      }
    } 
   PRK =
      (R.cwiseProduct(
          ((1 - beta.array()).log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        K.transpose() +

        R.cwiseProduct(
          (eta.array().log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        (
          ((1 - K.array()).matrix()).transpose()
        ) +

        (1 - R.array()).matrix().cwiseProduct(
          (beta.array().log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        K.transpose() +

        (1 - R.array()).matrix().cwiseProduct(
          ((1 - eta.array()).log().matrix()).replicate(1, R.rows()).transpose()
        ) *
        (((1 - K.array()).matrix()).transpose())
      ).array().exp().matrix();

    PR = PRK * PK;

    PKR = (PK * PR.cwiseInverse().transpose()).cwiseProduct(PRK.transpose());

  //Rcpp::Rcout << "1" << std::endl;
  //Rcpp::Rcout << typeid(beta).name() << std::endl;
  //Rcpp::Rcout << typeid(eta).name() << std::endl;
  //Rcpp::Rcout << typeid(mu0).name() << std::endl;
  //Rcpp::Rcout << typeid(mu1).name() << std::endl;
  //Rcpp::Rcout << typeid(iter).name() << std::endl;
  //Rcpp::Rcout << typeid(diff).name() << std::endl;

  return Rcpp::List::create(Rcpp::Named("P.K") = PK,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("maxiter") = iter >= maxiter,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("PK.R") = PKR
  );

}
