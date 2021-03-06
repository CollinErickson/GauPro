% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_White.R
\docType{class}
\name{White}
\alias{White}
\title{White noise Kernel R6 class}
\format{
\code{\link{R6Class}} object.
}
\value{
Object of \code{\link{R6Class}} with methods for fitting GP model.
}
\description{
White noise Kernel R6 class

White noise Kernel R6 class
}
\examples{
k1 <- White$new(s2=1e-8)
}
\section{Super class}{
\code{\link[GauPro:GauPro_kernel]{GauPro::GauPro_kernel}} -> \code{GauPro_kernel_White}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{s2}}{variance}

\item{\code{logs2}}{Log of s2}

\item{\code{logs2_lower}}{Lower bound of logs2}

\item{\code{logs2_upper}}{Upper bound of logs2}

\item{\code{s2_est}}{Should s2 be estimated?}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{White$new()}}
\item \href{#method-k}{\code{White$k()}}
\item \href{#method-kone}{\code{White$kone()}}
\item \href{#method-dC_dparams}{\code{White$dC_dparams()}}
\item \href{#method-C_dC_dparams}{\code{White$C_dC_dparams()}}
\item \href{#method-dC_dx}{\code{White$dC_dx()}}
\item \href{#method-param_optim_start}{\code{White$param_optim_start()}}
\item \href{#method-param_optim_start0}{\code{White$param_optim_start0()}}
\item \href{#method-param_optim_lower}{\code{White$param_optim_lower()}}
\item \href{#method-param_optim_upper}{\code{White$param_optim_upper()}}
\item \href{#method-set_params_from_optim}{\code{White$set_params_from_optim()}}
\item \href{#method-s2_from_params}{\code{White$s2_from_params()}}
\item \href{#method-clone}{\code{White$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="GauPro" data-topic="GauPro_kernel" data-id="plot">}\href{../../GauPro/html/GauPro_kernel.html#method-plot}{\code{GauPro::GauPro_kernel$plot()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize kernel object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$new(s2 = 1, D, s2_lower = 1e-08, s2_upper = 1e+08, s2_est = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{s2}}{Initial variance}

\item{\code{D}}{Number of input dimensions of data}

\item{\code{s2_lower}}{Lower bound for s2}

\item{\code{s2_upper}}{Upper bound for s2}

\item{\code{s2_est}}{Should s2 be estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-k"></a>}}
\if{latex}{\out{\hypertarget{method-k}{}}}
\subsection{Method \code{k()}}{
Calculate covariance between two points
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$k(x, y = NULL, s2 = self$s2, params = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector.}

\item{\code{y}}{vector, optional. If excluded, find correlation
of x with itself.}

\item{\code{s2}}{Variance parameter.}

\item{\code{params}}{parameters to use instead of beta and s2.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-kone"></a>}}
\if{latex}{\out{\hypertarget{method-kone}{}}}
\subsection{Method \code{kone()}}{
Find covariance of two points
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$kone(x, y, s2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector}

\item{\code{y}}{vector}

\item{\code{s2}}{Variance parameter}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dC_dparams"></a>}}
\if{latex}{\out{\hypertarget{method-dC_dparams}{}}}
\subsection{Method \code{dC_dparams()}}{
Derivative of covariance with respect to parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$dC_dparams(params = NULL, X, C_nonug, C, nug)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-C_dC_dparams"></a>}}
\if{latex}{\out{\hypertarget{method-C_dC_dparams}{}}}
\subsection{Method \code{C_dC_dparams()}}{
Calculate covariance matrix and its derivative
 with respect to parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$C_dC_dparams(params = NULL, X, nug)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-dC_dx"></a>}}
\if{latex}{\out{\hypertarget{method-dC_dx}{}}}
\subsection{Method \code{dC_dx()}}{
Derivative of covariance with respect to X
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$dC_dx(XX, X, s2 = self$s2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{XX}}{matrix of points}

\item{\code{X}}{matrix of points to take derivative with respect to}

\item{\code{s2}}{Variance parameter}

\item{\code{theta}}{Correlation parameters}

\item{\code{beta}}{log of theta}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-param_optim_start"></a>}}
\if{latex}{\out{\hypertarget{method-param_optim_start}{}}}
\subsection{Method \code{param_optim_start()}}{
Starting point for parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$param_optim_start(jitter = F, y, s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{jitter}}{Should there be a jitter?}

\item{\code{y}}{Output}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-param_optim_start0"></a>}}
\if{latex}{\out{\hypertarget{method-param_optim_start0}{}}}
\subsection{Method \code{param_optim_start0()}}{
Starting point for parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$param_optim_start0(jitter = F, y, s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{jitter}}{Should there be a jitter?}

\item{\code{y}}{Output}

\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-param_optim_lower"></a>}}
\if{latex}{\out{\hypertarget{method-param_optim_lower}{}}}
\subsection{Method \code{param_optim_lower()}}{
Lower bounds of parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$param_optim_lower(s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-param_optim_upper"></a>}}
\if{latex}{\out{\hypertarget{method-param_optim_upper}{}}}
\subsection{Method \code{param_optim_upper()}}{
Upper bounds of parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$param_optim_upper(s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{s2_est}}{Is s2 being estimated?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_params_from_optim"></a>}}
\if{latex}{\out{\hypertarget{method-set_params_from_optim}{}}}
\subsection{Method \code{set_params_from_optim()}}{
Set parameters from optimization output
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$set_params_from_optim(optim_out, s2_est = self$s2_est)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{optim_out}}{Output from optimization}

\item{\code{s2_est}}{s2 estimate}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-s2_from_params"></a>}}
\if{latex}{\out{\hypertarget{method-s2_from_params}{}}}
\subsection{Method \code{s2_from_params()}}{
Get s2 from params vector
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$s2_from_params(params, s2_est = self$s2_est)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{White$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
