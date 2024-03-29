% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotatePeaksNearRepeat.R
\name{annotatePeaksNearRepeat}
\alias{annotatePeaksNearRepeat}
\title{Obtain the position, type, name of the repeat elements that overlap a list of
peaks with variable widths, like ChIP-seq peaks}
\usage{
annotatePeaksNearRepeat(
  myPeakList,
  AnnotationData,
  minOverlap = 1L,
  output = c("hits", "all"),
  ...
)
}
\arguments{
\item{myPeakList}{A \link[GenomicRanges:GRanges-class]{GRanges} object}

\item{AnnotationData}{A \link[GenomicRanges:GRanges-class]{GRanges}}

\item{minOverlap}{The minimal degree of overlap between peak and Repeat element}

\item{output}{hits or all}
}
\description{
Obtain the position, type, name of the repeat elements that overlap a list of
peaks with variable widths, like ChIP-seq peaks
}
\examples{


    ## example 1: annotate myPeakList by Repeat.hg38, and output all peaks
    library(TECookbook)
    data(Reoeats.Human.hg38)
    annotatePeak = annotatePeaksNearRepeat(myPeakList, AnnotationData=Repeats.Human.hg38,
                   minOverlap=1, output = "all")

    ## example 2: annotate myPeakList (GRanges) and only output hits
    annotatePeak = annotatePeaksNearRepeat(myPeakList, AnnotationData=Repeats.Human.hg38,
                   minOverlap=1, output = "hits")

}
