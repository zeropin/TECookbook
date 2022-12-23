#' Obtain the position, type, name of the repeat elements that contain a list of
#' sites
#'
#' @param mySiteList A \link[GenomicRanges:GRanges-class]{GRanges} object
#' @param AnnotationData A \link[GenomicRanges:GRanges-class]{GRanges} or
#' @param SiteLocForDistance middle, start, or end
#' @param RepeatLocForDistance middle, start, or end
#' @param output either "all" or "hits" only
#' @param ignore.strand whether or not strand gets ignored
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
#'     ## example 1: annotate mySiteList by Repeat.hg38, and output all sites
#'     library(TECookbook)
#'     data(Reoeats.Human.hg38)
#'     annotateSite = annotateSitesInRepeat(mySiteList, AnnotationData=Repeats.Human.hg38,
#'                    SiteLocForDistance = "middle", RepeatLocForDistance = "middle",
#'                    output = "all")
#'
#'     ## example 2: annotate mySiteList (GRanges) and only output hits
#'     annotateSite = annotateSitesInRepeat(mySiteList, AnnotationData=Repeats.Human.hg38,
#'                    SiteLocForDistance = "middle", RepeatLocForDistance = "middle",
#'                    output = "hits")
#'
annotateSitesInRepeat <-
  function (mySiteList,
            AnnotationData,
            SiteLocForDistance=c("middle","start","end"),
            RepeatLocForDistance=c("middle","start","end"),
            output = c("hits", "all"),
            ignore.strand=TRUE, ...)
  {
    SiteLocForDistance = match.arg(SiteLocForDistance)
    RepeatLocForDistance = match.arg(RepeatLocForDistance)
   # out = match.arg(out)
    select = "all"

    if (missing(mySiteList))        stop("Missing required argument mySiteList!")

    if (!is(mySiteList, "GRanges")) stop("No valid mySiteList passed in. It needs to be GRanges object")

    if (missing(AnnotationData))    stop("No AnnotationData as GRanges is passed in")

    if (!inherits(AnnotationData, c("GRanges")))
      stop("AnnotationData needs to be GRanges object")

    if (!length(names(AnnotationData)) || any(is.na(names(AnnotationData))))
      names(AnnotationData) <- paste0("ann", seq_along(AnnotationData))


    if (is.null(names(AnnotationData))){
      names(AnnotationData) <- formatC(seq_along(AnnotationData),
                                       width=nchar(length(AnnotationData)),
                                       flag="0")
    }
    if (is.null(names(mySiteList))) {
      names(mySiteList) <- formatC(seq_along(mySiteList),
                                   width = nchar(length(mySiteList)),
                                   flag = "0")
    }
    if(any(duplicated(names(mySiteList)))){
      warning("Found duplicated names in mySiteList.
                    Changing the peak names ...")
      names(mySiteList) <- formatC(seq_along(mySiteList),
                                   width = nchar(length(mySiteList)),
                                   flag = "0")
    }
    savedNames <- names(mySiteList)

    ##clear seqnames, the format should be chr+NUM
    ##TODO, how about the seqname not start with chr?
    ## fix by seqlevelsStyle
    AnnotationData <- .formatSeqnames(AnnotationData, mySiteList)
    #mySiteList <- formatSeqnames(mySiteList)
    if(!all(seqlevels(mySiteList) %in% seqlevels(AnnotationData))){
      warning("not all the seqnames of mySiteList is
                    in the AnnotationData.")
    }

    if(length(mySiteList)>10000){
      ##huge dataset
      mySiteList <- split(mySiteList,
                          cut(seq_along(mySiteList),
                              ceiling(length(mySiteList)/5000)))
      mySiteList <- lapply(mySiteList, annotateSitesInRepeat,
                           AnnotationData=AnnotationData,
                           SiteLocForDistance=SiteLocForDistance,
                           RepeatLocForDistance=RepeatLocForDistance,
                           output = output,
                           ignore.strand=ignore.strand)
      names(mySiteList) <- NULL
      mySiteList <- unlist(GRangesList(mySiteList))
      names(mySiteList) <- make.names(paste(mySiteList$peak,
                                            mySiteList$feature))
      ##mySiteList
      return(mySiteList)
    }
    ###STEP1 get nearst annotation for each peak,
    ###use distanceToNearest(query, subject, ignore.strand=T/F, select)
    ## the distance got here is the shortestDistance
    ## select could only be arbitrary or all,
    ## if it is "first" or "last", use "all" instead.
    ## if output=="nearest", annotation should only consider the start point
    ##         ignore.strand <- all(strand(mySiteList)=="*") ||
    ##             all(strand(AnnotationData)=="*") ||
    ##             all(strand(mySiteList)=="+")

    featureGR <- AnnotationData
    end(featureGR) <-
      start(featureGR) <-
      switch(RepeatLocForDistance,
             middle=round(rowMeans(cbind(start(featureGR),
                                         end(featureGR)))),
             start=start(featureGR),
             end=end(featureGR)
      )
    myPeaksGR <- mySiteList
    end(myPeaksGR) <-
      start(myPeaksGR) <-
      switch(SiteLocForDistance,
             middle=round(rowMeans(cbind(start(myPeaksGR),
                                         end(myPeaksGR)))),
             start=start(myPeaksGR),
             end=end(myPeaksGR)
      )

    dist <- as.data.frame(findOverlaps(query = mySiteList, subject = AnnotationData,
                                       ignore.strand=ignore.strand,
                                       select="all",
                                       type="within"))
    dist$output <- rep("inside", nrow(dist))


    mySiteList.Hit <-
      mySiteList[dist$queryHits[!is.na(dist$subjectHits)]]
    mySiteList.NA <-
      mySiteList[!names(mySiteList) %in% names(mySiteList.Hit)]
    subjectHits <-
      AnnotationData[dist$subjectHits[!is.na(dist$subjectHits)]]
    mcols(subjectHits)$output <-
      dist[!is.na(dist$subjectHits),"output"]
    #    mySiteList.Hit$distanceToNearest <-
    #    dist$distance[!is.na(dist$subjectHits)]

    ###STEP2 get distance for each peak and nearest annotation by
    ## distance(x, y). the distance is calculated by
    ##        SiteLocForDistance=c("start","middle","end"),
    ##        RepeatLocForDistance=c("TSS","middle","start",
    ##                              "end", "geneEnd")
    RepeatLoc <-
      switch(RepeatLocForDistance,
             middle=as.integer(round(rowMeans(cbind(start(subjectHits),
                                                    end(subjectHits))))),
             start=start(subjectHits),
             end=end(subjectHits),
             geneEnd=as.integer(ifelse(strand(subjectHits)=="-",
                                       start(subjectHits),
                                       end(subjectHits))),
             TSS=as.integer(ifelse(strand(subjectHits)=="-",
                                   end(subjectHits),
                                   start(subjectHits))))

    SiteLoc <-
      switch(SiteLocForDistance,
             start=start(mySiteList.Hit),
             end=end(mySiteList.Hit),
             middle=as.integer(round(rowMeans(
               cbind(start(mySiteList.Hit), end(mySiteList.Hit))))))

    distancetoFeature <- as.integer(ifelse(strand(subjectHits)=="-",
                                           RepeatLoc - SiteLoc,
                                           SiteLoc - RepeatLoc))

    ###STEP3 relationship between query and subject:
    ###   "inside", "overlapEnd", "overlapStart",
    ###    "includeFeature", "upstream", "downstream"
    ### insideFeature <- getRelationship(mySiteList.Hit, subjectHits)
    mySiteList.Hit$Site <- names(mySiteList.Hit)
    mySiteList.Hit$Repeat_Class   <- subjectHits$Class
    mySiteList.Hit$Repeat_Family  <- subjectHits$Family
    mySiteList.Hit$Repeat_Name    <- subjectHits$Name
    mySiteList.Hit$Repeat_ID      <- names(subjectHits)
    #mySiteList.Hit$Repeat_start <- start(subjectHits)
    #mySiteList.Hit$Repeat_end <- end(subjectHits)
    #mySiteList.Hit$Repeat_strand <- as.character(strand(subjectHits))
    #mySiteList.Hit$insideFeature <- insideFeature[, "insideFeature"]
    mySiteList.Hit$distancetoFeature <- distancetoFeature
    #mySiteList.Hit$shortestDistance <-
     # as.integer(insideFeature[,"shortestDistance"])

    #mySiteList.Hit$fromOverlappingOrNearest <- subjectHits$output
    ##save oid for select == "first" or "last" filter
    #mySiteList.Hit$oid <- insideFeature[, "ss"]
    ###combind data of NA with Hits
    if(output=="all"){
      mySiteList.NA$Site <- names(mySiteList.NA)
      for(ncol in c("Repeat_Class", "Repeat_Family", "Repeat_Name",
                    "Repeat_ID",
                   # "Repat_start",
                  #  "Repeat_end",
                  #  "Repeat_strand",
                  #  "insideFeature",
                    "distancetoFeature"
                  # "shortestDistance",
                  #  "fromOverlappingOrNearest",
                  #  "oid"
                  ))
      mcols(mySiteList.NA)[, ncol] <- NA

      ## in case duplicated colnames
      colnames(mcols(mySiteList.NA)) <- colnames(mcols(mySiteList.Hit))
      mySiteList <- c(mySiteList.Hit, mySiteList.NA)
    }else{
      mySiteList <- mySiteList.Hit
    }


    ##remove column oid
    mySiteList$oid <- NULL
    ##re-order mySiteList as the original order
    oid <- seq_along(savedNames)
    names(oid) <- savedNames
    oid <- oid[names(mySiteList)]
    if(!any(is.na(oid))){
      mySiteList <- mySiteList[order(oid)]
    }


    names(mySiteList) <- mySiteList$Site
      #make.names(paste(mySiteList$peak, mySiteList$feature))

    ##mySiteList
    return(mySiteList)
  }

globalVariables(c("fromOverlappingOrNearest", "peak"))
