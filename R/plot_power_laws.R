
#' Plot a power law's probability density function (PDF)
#' 
#' Given a vector of discrete or continuous values, it fits a power law
#' distribution to the data and plots the PDF in log-log scale.
#' 
#' @param x numeric; A vector of discrete or continuous values.
#' @param b numeric; Bins have size b^(n-1).
#' @param ticks_x integer; Number of ticks in the x axis.
#' @param ticks_y integer; Number of ticks in the y axis.
#' @param xlab character; Label for the x axis.
#' @param ylab character; Label for the y axis.
#' 
#' @return A ggplot object with the plot.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Use log-binning to plot PDFs for the generated power law distributions
#' p_disc <- plot_pl_pdf(rnd_disc_pl(10000, 2, 2.3), 2)
#' p_cont <- plot_pl_pdf(rnd_cont_pl(10000, 2, 2.3), 2)
#' 
#' @export
#' @import ggplot2
#' @importFrom igraph fit_power_law
#' @importFrom dplyr tibble
#' @importFrom scales trans_breaks trans_format math_format
#'
plot_pl_pdf<- function(x, b = 2, 
                          ticks_x = 5, ticks_y = 5, 
                          xlab = "x", ylab = "p(x)"){
  
  pl <- fit_power_law(x)
  gma <- pl$alpha
  xmin <- pl$xmin
  
  max_pwr <- ceiling(log10(max(x))/log10(b))
  
  h <- graphics::hist(x, breaks = b^seq(0, max_pwr, 1), plot = FALSE)
  
  distr <- tibble(x = h$mids, p = h$density)
  
  intercept <- log10(distr$p[which(distr$x >= xmin)[1]]) + gma * log10(xmin)

  pdf_plot <- ggplot(distr, aes_(~x, ~p)) + geom_point() +
    geom_abline(slope = -gma, intercept = intercept, colour = "red") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = ticks_x), 
                  labels = trans_format("log10", math_format())) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = ticks_y), 
                  labels = trans_format("log10", math_format())) + 
    annotation_logticks() + labs(x = xlab, y = ylab) +
    theme_bw() + theme(panel.grid.minor = element_blank())

  return(pdf_plot)
}

#' Plot a power law's cumulative density function (CDF)
#' 
#' Given a vector of discrete or continuous values, it fits a power law
#' distribution to the data and plots the CDF in log-log scale.
#' 
#' @param x numeric; A vector of discrete or continuous values.
#' @param ticks_x integer; Number of ticks in the x axis.
#' @param ticks_y integer; Number of ticks in the y axis.
#' @param xlab character; Label for the x axis.
#' @param ylab character; Label for the y axis.
#' 
#' @return A ggplot object with the plot.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Use log-binning to plot CDFs for the generated power law distributions
#' p_disc <- plot_pl_cdf(rnd_disc_pl(10000, 2, 2.3))
#' p_cont <- plot_pl_cdf(rnd_cont_pl(10000, 2, 2.3))
#' 
#' @export
#' @import ggplot2
#' @import poweRlaw
#' @importFrom igraph fit_power_law
#' @importFrom dplyr tibble
#' @importFrom scales trans_breaks trans_format
#'
plot_pl_cdf <- function(x, ticks_x = 5, ticks_y = 5, 
                          xlab = "x", ylab = "P(X >= x)"){
  
  if(is.integer(x)){
    m_pl <- displ$new(x)
  }else{
    m_pl <- conpl$new(x)
  }
  fit <- fit_power_law(x)
  m_pl$setXmin(fit$xmin)
  m_pl$setPars(fit$alpha)
  xmin <- m_pl$getXmin()
  gma <- m_pl$getPars()
  
  cdf <- plot(m_pl, draw = FALSE)
  
  # Construct the linear fit (continuous expressions are used for both the
  # continuous and discrete cases)
  lfit <- tibble(k = 1:max(x), p = (1:max(x) / xmin)^(1 - gma))
  lfit$p <- lfit$p * cdf$y[which(cdf$x >= xmin)[1]]
  
  cdf_plot <- ggplot(cdf, aes_(~x, ~y)) + geom_point() +
    geom_line(data = lfit, aes_(~k, ~p), colour = "red") +
    coord_cartesian(xlim = c(min(cdf$x), max(cdf$x)),
                    ylim = c(min(cdf$y), 1)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = ticks_x), 
                  labels = trans_format("log10", math_format())) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = ticks_y),
                  labels = trans_format("log10", math_format())) + 
    annotation_logticks() + labs(x = xlab, y = ylab) +
    theme_bw() + theme(panel.grid.minor = element_blank())

  return(cdf_plot)
}

#' Random numbers from continuous power law distribution
#' 
#' Generate n power-law distributed continuous values.
#' 
#' @param n integer; The number of values to generate.
#' @param xmin numeric; Start of the power law behaviour.
#' @param gma numeric; Exponent of the power law.
#' 
#' @return A vector of n power-law distributed values.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Generate 10000 power-law distributed values
#' pl_cont <- rnd_cont_pl(10000, 2, 2.3)
#' 
#' @export
#'
rnd_cont_pl <- function(n, xmin = 2, gma = 2.5){
  return(xmin * (1 - stats::runif(n))^(-1/(gma - 1)))
}

#' Random numbers from discrete power law distribution
#' 
#' Generate n power-law distributed discrete values.
#' 
#' @param n integer; The number of values to generate.
#' @param xmin numeric; Start of the power law behaviour.
#' @param gma numeric; Exponent of the power law.
#' 
#' @return A vector of n power-law distributed values.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Generate 10000 power-law distributed values
#' pl_disc <- rnd_disc_pl(10000, 2, 2.3)
#' 
#' @export
#'
rnd_disc_pl <- function(n, xmin = 2, gma = 2.5){
  return(round(xmin * (1 - stats::runif(n))^(-1/(gma - 1))))
}