% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_ignoreinds.R
\docType{class}
\name{IgnoreIndsKernel}
\alias{IgnoreIndsKernel}
\title{Kernel R6 class}
\format{
\code{\link{R6Class}} object.
}
\value{
Object of \code{\link{R6Class}} with methods for fitting GP model.
}
\description{
Kernel R6 class

Kernel R6 class
}
\examples{
kg <- Gaussian$new(D=3)
kig <- GauPro::IgnoreIndsKernel$new(k = Gaussian$new(D=3), ignoreinds = 2)
Xtmp <- as.matrix(expand.grid(1:2, 1:2, 1:2))
cbind(Xtmp, kig$k(Xtmp))
cbind(Xtmp, kg$k(Xtmp))
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{D}}{Number of input dimensions of data}

\item{\code{kernel}}{Kernel to use on indices that aren't ignored}

\item{\code{ignoreinds}}{Indices to ignore. For a matrix X, these are the columns
to ignore. For example, when those dimensions will be given a different
kernel, such as for factors.}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{s2_est}}{Is s2 being estimated?}

\item{\code{s2}}{Value of s2 (variance)}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{IgnoreIndsKernel$new()}}
\item \href{#method-k}{\code{IgnoreIndsKernel$k()}}
\item \href{#method-kone}{\code{IgnoreIndsKernel$kone()}}
\item \href{#method-dC_dparams}{\code{IgnoreIndsKernel$dC_dparams()}}
\item \href{#method-C_dC_dparams}{\code{IgnoreIndsKernel$C_dC_dparams()}}
\item \href{#method-param_optim_start}{\code{IgnoreIndsKernel$param_optim_start()}}
\item \href{#method-param_optim_start0}{\code{IgnoreIndsKernel$param_optim_start0()}}
\item \href{#method-param_optim_lower}{\code{IgnoreIndsKernel$param_optim_lower()}}
\item \href{#method-param_optim_upper}{\code{IgnoreIndsKernel$param_optim_upper()}}
\item \href{#method-set_params_from_optim}{\code{IgnoreIndsKernel$set_params_from_optim()}}
\item \href{#method-s2_from_params}{\code{IgnoreIndsKernel$s2_from_params()}}
\item \href{#method-clone}{\code{IgnoreIndsKernel$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize kernel object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$new(k, ignoreinds)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{k}}{Kernel to use on the non-ignored indices}

\item{\code{ignoreinds}}{Indices of columns of X to ignore.}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$k(x, y = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector.}

\item{\code{y}}{vector, optional. If excluded, find correlation
of x with itself.}

\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$kone(x, y, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector}

\item{\code{y}}{vector}

\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$dC_dparams(params = NULL, X, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{params}}{Kernel parameters}

\item{\code{X}}{matrix of points in rows}

\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$C_dC_dparams(params = NULL, X, nug)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-param_optim_start"></a>}}
\if{latex}{\out{\hypertarget{method-param_optim_start}{}}}
\subsection{Method \code{param_optim_start()}}{
Starting point for parameters for optimization
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$param_optim_start(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$param_optim_start0(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$param_optim_lower(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$param_optim_upper(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$set_params_from_optim(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$s2_from_params(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Passed to kernel}
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
\if{html}{\out{<div class="r">}}\preformatted{IgnoreIndsKernel$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
