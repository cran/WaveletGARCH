\name{WaveletGARCHFit}
\alias{WaveletGARCHFit}
\alias{print.WaveletGARCHFit}
\title{Fitting of Wavelet-GARCH model}
\description{
Fitting of Wavelet-GARCH model based on ARCH LM test.
}
\usage{
WaveletGARCHFit(series,filtern,level)
\method{print}{WaveletGARCHFit}(x,\dots)
}
\arguments{
  \item{series}{univariate time series}
  \item{filtern}{The name of wavelet filter}
  \item{level}{The level of wavelet decomposition}
  \item{x}{An object of WaveletGARCHFit}
  \item{\dots}{Additional arguments if any}
     }
\value{
    \item{fittedobject}{The fitted value of the series by Waveetl-GARCH model}
  }

\references{
Percival D. B. and Walden A. T. 2000. Wavelet Methods for Time-Series Analysis. Cambridge Univ. Press, U.K.

Paul R. K., Prajneshu and Ghosh H. 2013. Wavelet Frequency Domain Approach for Modelling and Forecasting of Indian Monsoon Rainfall Time-Series Data. Journal of the Indian society of agricultural statistics, 67, 319 to 327.

Paul, R.K. and Birthal, P.S. 2015. Investigating rainfall trend over India using wavelet technique. Journal of Water and Climate Change, 7, 365 to 378.

Paul, R. K. 2015. ARIMAX-GARCH-WAVELET Model for forecasting volatile data. Model Assisted Statistics and Application, 10, 243 to252.
}
\examples{

data(mtcars)
ab<-mtcars$qsec

objfit<-WaveletGARCHFit(ab,"d4",4)
}

