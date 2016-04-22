\docType{package}
\name{fusedanova-package}
\alias{fusedanova-package}
\title{Fused-ANOVA package: general presentation}
\description{
This package is designed to fit accurately the Fused-ANOVA
model, a penalized method to solve the one-way ANOVA
problem by collapsing the coefficients of \eqn{K}{K}
conditions.  For a large class of weights implemented here,
our homotopy algorithm is in
\eqn{\mathcal{O}(K\log(K))}{O(klog(K))}.  These weights
induce a balanced tree structure and simplify the
interpretation of the results. The package contains an
illustrating phenotypic data set: given a trait, we
reconstruct a balanced tree structure and assess its
agreement with the known phylogeny.
}
\section{Problem solved}{
  The optimization problem solved by fused-ANOVA is
  \if{latex}{\deqn{% \hat{\beta}_{\lambda} = \arg
  \min_{\beta} \left\{\sum_{k=1}^K \sum_{i=1}^{n_k}
  \left(Y_{ik}-\beta_k \right)^2 + \lambda \sum_{k,\ell}
  w_{kl} \left|\beta_k - \beta_\ell \right|\right\}}}
  \if{html}{\out{ <center> &beta;<sup>hat</sup>
  <sub>&lambda;<sub>1</sub></sub> = argmin<sub>&beta;</sub>
  sum<sub>k</sub> sum_i (Y<sub>ik</sub> -
  &beta<sub>k</sub>)<sup>2</sup> + &lambda;
  sum<sub>k,l</sub> w<sub>k,l</sub> &#124;
  &beta;<sub>k</sub> - &beta;<sub>l</sub> &#124;, </center>
  }} \if{text}{\deqn{beta.hat(lambda) = argmin_beta sum_k
  sum_i (Y_ik - beta_k)^2 + lambda sum_k sum_l w_kl |
  beta_k - beta_l|,}}

  where \eqn{Y_{ik}}{Y_ik} is the intensity of a continuous
  random variable for sample \eqn{i}{i} in condition
  \eqn{k}{k} and \eqn{\beta_k}{beta_k} is the mean
  parameter of condition \eqn{k}{k}. We denote by
  \eqn{K}{K} the total number of conditions and
  \eqn{n_k}{n_k} the number of sample in each condition.
}

\section{Choice of weights and performance of the algorithm}{
  For various weights in the fused-penalty (entailing
  "laplace", "gaussian", "default", "adaptive" - see the
  corresponding documentation), the homotopy algorithm
  produces a path that contains no split, which is highly
  desirable since in this case \enumerate{%

  \item the order of the \eqn{\beta_k}{beta_k} always
  matches the order of the empirical mean of each
  condition;

  \item the recovered structure is a tree which simplifies
  the interpretation;

  \item the total number of iterations is guaranteed to be
  small and equal to \eqn{K}{K};

  \item we avoid maximum flow problems whose resolution is
  computationally demanding.

  }

  The associated algorithm is in
  \eqn{\mathcal{O}(K\log(K))}{O(klog(K))}. In this
  perspective, we extend the work of Hocking et al. to a
  larger class of weights.

  For other weights, split can occur along the path of
  solution. We adapted the algorithm developed by Hoefling
  (reference below) to the fused-ANOVA problem.
}

\section{Efficient cross-validation procedure}{
  We provide a fast cross validation (CV) procedure to
  select \eqn{\lambda}{lambda} for both the general and the
  no split algorithms.  The idea behind this procedure is
  to take advantage of the DAG structure of the path of
  solutions along \eqn{\lambda}{lambda}.  Rather than
  computing the CV error for each condition separately, we
  traverse each edge of the DAG once and only once and
  compute simultaneously the error of all conditions going
  through this edge.  If we consider a perfectly balanced
  tree and a grid of \eqn{P}{P} values of
  \eqn{\lambda}{lambda} we achieve \eqn{\mathcal{O}(P \log
  (P))}{O(P log (P))} rather than a
  \eqn{\mathcal{O}(P^2)}{O(P^2)} complexity.
}

\section{Technical remarks}{
  Most of the numerical work is done in C++, relying on the
  \pkg{Rcpp} package. We also use the multi-core capability
  of the computer through the \pkg{parallel} package when
  multiple variables are to be classified. This feature is
  not available for Windows user, though.
}
\author{
Pierre Gutierrez, Julien Chiquet, Guillem Rigaill.
}
\references{
Fused-ANOVA: shortly coming

H. Hoefling. A path algorithm for the fused lasso signal
approximator, technical report, arXiv, 2010.

T. Hocking, J.-P. Vert, F. Bach, and A. Joulin.
Clusterpath: an Algorithm for Clustering using Convex
Fusion Penalties, ICML, 2011.
}

