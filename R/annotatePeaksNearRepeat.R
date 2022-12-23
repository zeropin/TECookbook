#' Obtain the position, type, name of the repeat elements that overlap a list of
#' peaks with variable widths, like ChIP-seq peaks
#'
#' @param myPeakList A \link[GenomicRanges:GRanges-class]{GRanges} object
#' @param AnnotationData A \link[GenomicRanges:GRanges-class]{GRanges}
#' @param minOverlap The minimal degree of overlap between peak and Repeat element
#' @param output hits or all
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom BiocGenerics start end width strand
#' @importFrom S4Vectors mcols
#' @importFrom dplyr %>% filter group_by top_n
#' @examples
#'
#'
#'     ## example 1: annotate myPeakList by Repeat.hg38, and output all peaks
#'     library(TECookbook)
#'     data(Reoeats.Human.hg38)
#'     annotatePeak = annotatePeaksNearRepeat(myPeakList, AnnotationData=Repeats.Human.hg38,
#'                    minOverlap=1, output = "all")
#'
#'     ## example 2: annotate myPeakList (GRanges) and only output hits
#'     annotatePeak = annotatePeaksNearRepeat(myPeakList, AnnotationData=Repeats.Human.hg38,
#'                    minOverlap=1, output = "hits")
#'
annotatePeaksNearRepeat <-
    function (myPeakList,
              AnnotationData,
              minOverlap = 1L,
              output = c("hits", "all"),
              ...)
    {
        if (missing(myPeakList))        stop("Missing required argument myPeakList!")

        if (!is(myPeakList, "GRanges")) stop("No valid myPeakList passed in. It needs to be GRanges object")

        if (missing(AnnotationData))    stop("No AnnotationData as GRanges is passed in")

        if (!inherits(AnnotationData, c("GRanges")))
            stop("AnnotationData needs to be GRanges object")

       if (!length(names(AnnotationData)) || any(is.na(names(AnnotationData))))
       {
             names(AnnotationData) <- paste0("ann", seq_along(AnnotationData))
       }

        nAnno <- length(AnnotationData)

        rm(nAnno)


      if (is.null(names(AnnotationData))){
            names(AnnotationData) <- formatC(seq_along(AnnotationData),
                                          width=nchar(length(AnnotationData)),
                                          flag="0")
        }
      if (is.null(names(myPeakList))) {
            names(myPeakList) <- formatC(seq_along(myPeakList),
                                         width = nchar(length(myPeakList)),
                                         flag = "0")
        }
      if(any(duplicated(names(myPeakList)))){
            warning("Found duplicated names in myPeakList.
                    Changing the peak names ...")
            names(myPeakList) <- formatC(seq_along(myPeakList),
                                         width = nchar(length(myPeakList)),
                                         flag = "0")
       }
        savedNames <- names(myPeakList)

        ##clear seqnames, the format should be chr+NUM
        ##TODO, how about the seqname not start with chr?
        ## fix by seqlevelsStyle
        AnnotationData <- .formatSeqnames(AnnotationData, myPeakList)
        #myPeakList <- formatSeqnames(myPeakList)
      if(!all(seqlevels(myPeakList) %in% seqlevels(AnnotationData))){
            warning("not all the seqnames of myPeakList is
                    in the AnnotationData.")
      }

      if(length(myPeakList)>10000){
            ##huge dataset
            myPeakList <- split(myPeakList,
                                cut(seq_along(myPeakList),
                                    ceiling(length(myPeakList)/5000)))
            myPeakList <- lapply(myPeakList, annotatePeaksNearRepeat,
                                 AnnotationData=AnnotationData,
                                 minOverlap = minOverlap,
                                 output = output)
            names(myPeakList) <- NULL
            myPeakList <- unlist(GRangesList(myPeakList))
            names(myPeakList) <- make.names(paste(myPeakList$peak,
                                                  myPeakList$feature))
            ##myPeakList
            return(myPeakList)
      }
        ###STEP1 get nearst annotation for each peak,
      ###use distanceToNearest(query, subject, ignore.strand=T/F, select)
      ## the distance got here is the shortestDistance
      ## select could only be arbitrary or all,
      ## if it is "first" or "last", use "all" instead.
      ## if output=="nearest", annotation should only consider the start point
      ##         ignore.strand <- all(strand(myPeakList)=="*") ||
      ##             all(strand(AnnotationData)=="*") ||
      ##             all(strand(myPeakList)=="+")

      featureGR <- AnnotationData
      end(featureGR) <-
        start(featureGR) <- round(rowMeans(cbind(start(featureGR),
                                             end(featureGR))))

      myPeaksGR <- myPeakList
      end(myPeaksGR) <-
        start(myPeaksGR) <- round(rowMeans(cbind(start(myPeaksGR),
                                             end(myPeaksGR))))

      dist <- as.data.frame(findOverlaps(myPeakList, AnnotationData,
                                               maxgap=-minOverlap,
                                               ignore.strand=TRUE,
                                               select="all"))

      if(nrow(dist)==0) dist[1,] <- NA
      if(ncol(dist)==1){
              dist <- cbind(queryHits=seq_along(myPeakList), subjectHits=dist)
              colnames(dist) <- c("queryHits", "subjectHits")
      }
            dist$output <- rep("Overlapping", nrow(dist))

        ##  nearest is NOT filtered by maxgap, is this should be changed?
        ##        dist <- dist[!is.na(dist$subjectHits), ]
        ##        distance <- distance(myPeakList[dist$queryHits],
        ##                             AnnotationData[dist$subjectHits],
        ##                             ignore.strand=ignore.strand)
        ##        dist <- dist[abs(distance) <= maxgap, ]

        myPeakList.Hit <-
            myPeakList[dist$queryHits[!is.na(dist$subjectHits)]]
        myPeakList.NA <-
            myPeakList[!names(myPeakList) %in% names(myPeakList.Hit)]
        subjectHits <-
            AnnotationData[dist$subjectHits[!is.na(dist$subjectHits)]]
        mcols(subjectHits)$output <-
            dist[!is.na(dist$subjectHits),"output"]
        #    myPeakList.Hit$distanceToNearest <-
        #    dist$distance[!is.na(dist$subjectHits)]

        ###STEP2 get distance for each peak and nearest annotation by
        ## distance(x, y). the distance is calculated by
        ##        PeakLocForDistance=c("start","middle","end"),
        ##        RepeatLocForDistance=c("TSS","middle","start",
        ##                              "end", "geneEnd")
        RepeatLoc <- as.integer(round(rowMeans(cbind(start(subjectHits),
                                                          end(subjectHits)))))
        PeakLoc   <- as.integer(round(rowMeans(cbind(start(myPeakList.Hit),
                                                     end(myPeakList.Hit)))))

        distancetoFeature <- abs(RepeatLoc - PeakLoc)

        ###STEP3 relationship between query and subject:
        ###   "inside", "overlapEnd", "overlapStart",
        ###    "includeFeature", "upstream", "downstream"
        ### insideFeature <- getRelationship(myPeakList.Hit, subjectHits)
        myPeakList.Hit$Peak           <- names(myPeakList.Hit)
        myPeakList.Hit$Repeat_Class   <- subjectHits$Class
        myPeakList.Hit$Repeat_Family  <- subjectHits$Family
        myPeakList.Hit$Repeat_Name    <- subjectHits$Name
        myPeakList.Hit$Repeat_ID      <- names(subjectHits)
        myPeakList.Hit$distancetoFeature <- distancetoFeature

        #myPeakList.Hit$fromOverlappingOrNearest <- subjectHits$output
        ##save oid for select == "first" or "last" filter
        #myPeakList.Hit$oid <- insideFeature[, "ss"]
        ###combind data of NA with Hits

        if(output=="all"){
          myPeakList.NA$Peak <- names(myPeakList.NA)
          for(ncol in c("Repeat_Class", "Repeat_Family", "Repeat_Name",
                        "Repeat_ID",
                        "distancetoFeature"
                  ))
        mcols(myPeakList.NA)[, ncol] <- NA

        ## in case duplicated colnames
        colnames(mcols(myPeakList.NA)) <- colnames(mcols(myPeakList.Hit))
        myPeakList <- c(myPeakList.Hit, myPeakList.NA)
        }else{
          myPeakList <- myPeakList.Hit
        }

        ##remove column oid
        myPeakList$oid <- NULL
        ##re-order myPeakList as the original order
        oid <- seq_along(savedNames)
        names(oid) <- savedNames
        oid <- oid[names(myPeakList)]
        if(!any(is.na(oid))){
            myPeakList <- myPeakList[order(oid)]
        }

        names(myPeakList) <-
            make.names(paste(myPeakList$peak, myPeakList$feature))

        ##myPeakList
        return(myPeakList)
    }

globalVariables(c("fromOverlappingOrNearest", "peak"))
