% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stomp-par.R, R/stomp.R
\name{stomp_par}
\alias{stomp_par}
\alias{stomp}
\title{Univariate STOMP algorithm}
\usage{
stomp_par(
  ...,
  window_size,
  exclusion_zone = getOption("tsmp.exclusion_zone", 1/2),
  verbose = getOption("tsmp.verbose", 2),
  n_workers = 2
)

stomp(
  ...,
  window_size,
  exclusion_zone = getOption("tsmp.exclusion_zone", 1/2),
  verbose = getOption("tsmp.verbose", 2)
)
}
\arguments{
\item{\dots}{a \code{matrix} or a \code{vector}. If a second time series is supplied it will be a join matrix
profile.}

\item{window_size}{an \code{int}. Size of the sliding window.}

\item{exclusion_zone}{a \code{numeric}. Size of the exclusion zone, based on window size (default is
\code{1/2}). See details.}

\item{verbose}{an \code{int}. See details. (Default is \code{2}).}

\item{n_workers}{an \code{int}. Number of workers for parallel. (Default is \code{2}).}
}
\value{
Returns a \code{MatrixProfile} object, a \code{list} with the matrix profile \code{mp}, profile index \code{pi}
left and right matrix profile \code{lmp}, \code{rmp} and profile index \code{lpi}, \code{rpi}, window size \code{w} and
exclusion zone \code{ez}.
}
\description{
Computes the Matrix Profile and Profile Index for Univariate Time Series.
}
\details{
The Matrix Profile, has the potential to revolutionize time series data mining because of its
generality, versatility, simplicity and scalability. In particular it has implications for time
series motif discovery, time series joins, shapelet discovery (classification), density
estimation, semantic segmentation, visualization, rule discovery, clustering etc. \code{verbose}
changes how much information is printed by this function; \code{0} means nothing, \code{1} means text, \code{2}
adds the progress bar, \code{3} adds the finish sound. \code{exclusion_zone} is used to avoid  trivial
matches; if a query data is provided (join similarity), this parameter is ignored.
}
\section{Functions}{
\itemize{
\item \code{stomp_par()}: Parallel version.

\item \code{stomp()}: Single thread version.

}}
\examples{
mp <- stomp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
\donttest{
#' # using threads
mp <- stomp_par(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)

ref_data <- mp_toy_data$data[, 1]
query_data <- mp_toy_data$data[, 2]
# self similarity
mp <- stomp(ref_data, window_size = 30)
# join similarity
mp2 <- stomp(ref_data, query_data, window_size = 30)
}
}
\references{
\itemize{
\item Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting
a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and
Joins. Icdm. 2016 Jan 22;54(1):739-48.
}

Website: \url{http://www.cs.ucr.edu/~eamonn/MatrixProfile.html}
}
\seealso{
Other matrix profile computations: 
\code{\link{mstomp_par}()},
\code{\link{scrimp}()},
\code{\link{stamp_par}()},
\code{\link{tsmp}()},
\code{\link{valmod}()}
}
\concept{matrix profile computations}
