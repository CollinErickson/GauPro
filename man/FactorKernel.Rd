% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_Factor.R
\docType{class}
\name{FactorKernel}
\alias{FactorKernel}
\alias{k_FactorKernel}
\title{Factor Kernel R6 class}
\format{
\code{\link[R6]{R6Class}} object.
}
\usage{
k_FactorKernel(
  s2 = 1,
  D,
  nlevels,
  xindex,
  p_lower = 0,
  p_upper = 0.9,
  p_est = TRUE,
  s2_lower = 1e-08,
  s2_upper = 1e+08,
  s2_est = TRUE,
  p,
  useC = TRUE,
  offdiagequal = 1 - 1e-06
)
}
\arguments{
\item{s2}{Initial variance}

\item{D}{Number of input dimensions of data}

\item{nlevels}{Number of levels for the factor}

\item{xindex}{Index of the factor (which column of X)}

\item{p_lower}{Lower bound for p}

\item{p_upper}{Upper bound for p}

\item{p_est}{Should p be estimated?}

\item{s2_lower}{Lower bound for s2}

\item{s2_upper}{Upper bound for s2}

\item{s2_est}{Should s2 be estimated?}

\item{p}{Vector of correlations}

\item{useC}{Should C code used? Not implemented for FactorKernel yet.}

\item{offdiagequal}{What should offdiagonal values be set to when the
indices are the same? Use to avoid decomposition errors, similar to
adding a nugget.}
}
\value{
Object of \code{\link[R6]{R6Class}} with methods for fitting GP model.
}
\description{
Initialize kernel object
}
\details{
For a factor that has been converted to its indices.
Each factor will need a separate kernel.
}
\examples{
kk <- FactorKernel$new(D=1, nlevels=5, xindex=1)
kk$p <- (1:10)/100
kmat <- outer(1:5, 1:5, Vectorize(kk$k))
kmat
kk$plot()


# 2D, Gaussian on 1D, index on 2nd dim
if (requireNamespace("dplyr", quietly=TRUE)) {
library(dplyr)
n <- 20
X <- cbind(matrix(runif(n,2,6), ncol=1),
           matrix(sample(1:2, size=n, replace=TRUE), ncol=1))
X <- rbind(X, c(3.3,3))
n <- nrow(X)
Z <- X[,1] - (X[,2]-1.8)^2 + rnorm(n,0,.1)
tibble(X=X, Z) \%>\% arrange(X,Z)
k2a <- IgnoreIndsKernel$new(k=Gaussian$new(D=1), ignoreinds = 2)
k2b <- FactorKernel$new(D=2, nlevels=3, xind=2)
k2 <- k2a * k2b
k2b$p_upper <- .65*k2b$p_upper
gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5,
                              nug.min=1e-2, restarts=0)
gp$kernel$k1$kernel$beta
gp$kernel$k2$p
gp$kernel$k(x = gp$X)
tibble(X=X, Z=Z, pred=gp$predict(X)) \%>\% arrange(X, Z)
tibble(X=X[,2], Z) \%>\% group_by(X) \%>\% summarize(n=n(), mean(Z))
curve(gp$pred(cbind(matrix(x,ncol=1),1)),2,6, ylim=c(min(Z), max(Z)))
points(X[X[,2]==1,1], Z[X[,2]==1])
curve(gp$pred(cbind(matrix(x,ncol=1),2)), add=TRUE, col=2)
points(X[X[,2]==2,1], Z[X[,2]==2], col=2)
curve(gp$pred(cbind(matrix(x,ncol=1),3)), add=TRUE, col=3)
points(X[X[,2]==3,1], Z[X[,2]==3], col=3)
legend(legend=1:3, fill=1:3, x="topleft")
# See which points affect (5.5, 3 themost)
data.frame(X, cov=gp$kernel$k(X, c(5.5,3))) \%>\% arrange(-cov)
plot(k2b)
}

}
\section{Super class}{
\code{\link[GauPro:GauPro_kernel]{GauPro::GauPro_kernel}} -> \code{GauPro_kernel_FactorKernel}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{p}}{Parameter for correlation}

\item{\code{p_est}}{Should p be estimated?}

\item{\code{p_lower}}{Lower bound of p}

\item{\code{p_upper}}{Upper bound of p}

\item{\code{p_length}}{length of p}

\item{\code{s2}}{variance}

\item{\code{s2_est}}{Is s2 estimated?}

\item{\code{logs2}}{Log of s2}

\item{\code{logs2_lower}}{Lower bound of logs2}

\item{\code{logs2_upper}}{Upper bound of logs2}

\item{\code{xindex}}{Index of the factor (which column of X)}

\item{\code{nlevels}}{Number of levels for the factor}

\item{\code{offdiagequal}}{What should offdiagonal values be set to when the
indices are the same? Use to avoid decomposition errors, similar to
adding a nugget.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-GauPro_kernel_FactorKernel-new}{\code{FactorKernel$new()}}
\item \href{#method-GauPro_kernel_FactorKernel-k}{\code{FactorKernel$k()}}
\item \href{#method-GauPro_kernel_FactorKernel-kone}{\code{FactorKernel$kone()}}
\item \href{#method-GauPro_kernel_FactorKernel-dC_dparams}{\code{FactorKernel$dC_dparams()}}
\item \href{#method-GauPro_kernel_FactorKernel-C_dC_dparams}{\code{FactorKernel$C_dC_dparams()}}
\item \href{#method-GauPro_kernel_FactorKernel-dC_dx}{\code{FactorKernel$dC_dx()}}
\item \href{#method-GauPro_kernel_FactorKernel-param_optim_start}{\code{FactorKernel$param_optim_start()}}
\item \href{#method-GauPro_kernel_FactorKernel-param_optim_start0}{\code{FactorKernel$param_optim_start0()}}
\item \href{#method-GauPro_kernel_FactorKernel-param_optim_lower}{\code{FactorKernel$param_optim_lower()}}
\item \href{#method-GauPro_kernel_FactorKernel-param_optim_upper}{\code{FactorKernel$param_optim_upper()}}
\item \href{#method-GauPro_kernel_FactorKernel-set_params_from_optim}{\code{FactorKernel$set_params_from_optim()}}
\item \href{#method-GauPro_kernel_FactorKernel-s2_from_params}{\code{FactorKernel$s2_from_params()}}
\item \href{#method-GauPro_kernel_FactorKernel-print}{\code{FactorKernel$print()}}
\item \href{#method-GauPro_kernel_FactorKernel-clone}{\code{FactorKernel$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="GauPro" data-topic="GauPro_kernel" data-id="plot"><a href='../../GauPro/html/GauPro_kernel.html#method-GauPro_kernel-plot'><code>GauPro::GauPro_kernel$plot()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-new"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-new}{}}}
\subsection{Method \code{new()}}{
Initialize kernel object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$new(
  s2 = 1,
  D,
  nlevels,
  xindex,
  p_lower = 0,
  p_upper = 0.9,
  p_est = TRUE,
  s2_lower = 1e-08,
  s2_upper = 1e+08,
  s2_est = TRUE,
  p,
  useC = TRUE,
  offdiagequal = 1 - 1e-06
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{s2}}{Initial variance}

\item{\code{D}}{Number of input dimensions of data}

\item{\code{nlevels}}{Number of levels for the factor}

\item{\code{xindex}}{Index of the factor (which column of X)}

\item{\code{p_lower}}{Lower bound for p}

\item{\code{p_upper}}{Upper bound for p}

\item{\code{p_est}}{Should p be estimated?}

\item{\code{s2_lower}}{Lower bound for s2}

\item{\code{s2_upper}}{Upper bound for s2}

\item{\code{s2_est}}{Should s2 be estimated?}

\item{\code{p}}{Vector of correlations}

\item{\code{useC}}{Should C code used? Not implemented for FactorKernel yet.}

\item{\code{offdiagequal}}{What should offdiagonal values be set to when the
indices are the same? Use to avoid decomposition errors, similar to
adding a nugget.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-k"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-k}{}}}
\subsection{Method \code{k()}}{
Calculate covariance between two points
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$k(x, y = NULL, p = self$p, s2 = self$s2, params = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector.}

\item{\code{y}}{vector, optional. If excluded, find correlation
of x with itself.}

\item{\code{p}}{Correlation parameters.}

\item{\code{s2}}{Variance parameter.}

\item{\code{params}}{parameters to use instead of beta and s2.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-kone"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-kone}{}}}
\subsection{Method \code{kone()}}{
Find covariance of two points
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$kone(x, y, p, s2, isdiag = TRUE, offdiagequal = self$offdiagequal)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector}

\item{\code{y}}{vector}

\item{\code{p}}{correlation parameters on regular scale}

\item{\code{s2}}{Variance parameter}

\item{\code{isdiag}}{Is this on the diagonal of the covariance?}

\item{\code{offdiagequal}}{What should offdiagonal values be set to when the
indices are the same? Use to avoid decomposition errors, similar to
adding a nugget.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-dC_dparams"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-dC_dparams}{}}}
\subsection{Method \code{dC_dparams()}}{
Derivative of covariance with respect to parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$dC_dparams(params = NULL, X, C_nonug, C, nug)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{params}}{Kernel parameters}

\item{\code{X}}{matrix of points in rows}

\item{\code{C_nonug}}{Covariance without nugget added to diagonal}

\item{\code{C}}{Covariance with nugget}

\item{\code{nug}}{Value of nugget}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-C_dC_dparams"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-C_dC_dparams}{}}}
\subsection{Method \code{C_dC_dparams()}}{
Calculate covariance matrix and its derivative
 with respect to parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$C_dC_dparams(params = NULL, X, nug)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{params}}{Kernel parameters}

\item{\code{X}}{matrix of points in rows}

\item{\code{nug}}{Value of nugget}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-dC_dx"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-dC_dx}{}}}
\subsection{Method \code{dC_dx()}}{
Derivative of covariance with respect to X
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$dC_dx(XX, X, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{XX}}{matrix of points}

\item{\code{X}}{matrix of points to take derivative with respect to}

\item{\code{...}}{Additional args, not used}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-param_optim_start"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-param_optim_start}{}}}
\subsection{Method \code{param_optim_start()}}{
Starting point for parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$param_optim_start(
  jitter = F,
  y,
  p_est = self$p_est,
  s2_est = self$s2_est
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{jitter}}{Should there be a jitter?}

\item{\code{y}}{Output}

\item{\code{p_est}}{Is p being estimated?}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-param_optim_start0"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-param_optim_start0}{}}}
\subsection{Method \code{param_optim_start0()}}{
Starting point for parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$param_optim_start0(
  jitter = F,
  y,
  p_est = self$p_est,
  s2_est = self$s2_est
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{jitter}}{Should there be a jitter?}

\item{\code{y}}{Output}

\item{\code{p_est}}{Is p being estimated?}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-param_optim_lower"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-param_optim_lower}{}}}
\subsection{Method \code{param_optim_lower()}}{
Lower bounds of parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$param_optim_lower(p_est = self$p_est, s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{p_est}}{Is p being estimated?}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-param_optim_upper"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-param_optim_upper}{}}}
\subsection{Method \code{param_optim_upper()}}{
Upper bounds of parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$param_optim_upper(p_est = self$p_est, s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{p_est}}{Is p being estimated?}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-set_params_from_optim"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-set_params_from_optim}{}}}
\subsection{Method \code{set_params_from_optim()}}{
Set parameters from optimization output
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$set_params_from_optim(
  optim_out,
  p_est = self$p_est,
  s2_est = self$s2_est
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{optim_out}}{Output from optimization}

\item{\code{p_est}}{Is p being estimated?}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-s2_from_params"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-s2_from_params}{}}}
\subsection{Method \code{s2_from_params()}}{
Get s2 from params vector
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$s2_from_params(params, s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{params}}{parameter vector}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-print"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-print}{}}}
\subsection{Method \code{print()}}{
Print this object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GauPro_kernel_FactorKernel-clone"></a>}}
\if{latex}{\out{\hypertarget{method-GauPro_kernel_FactorKernel-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FactorKernel$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
