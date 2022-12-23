###clear seqnames, the format should be chr+NUM
#' @importFrom GenomeInfoDb `seqlevels<-` `seqlevelsStyle` `seqlevelsStyle<-`
#' @NoRd
.formatSeqnames <- function(from, to) {
  forceFormatSeqnames <- function(from, to) {
    message("\n Try to keep the seqname style consistent.")
    seql <- seqlevels(to)
    getPrefix <- function(x, seql) {
      seql <- table(c(grepl(x, seql), "TRUE", "FALSE"))
      if (seql['TRUE'] > seql['FALSE']) {
        prefix <- x
      } else{
        prefix <- ""
      }
      prefix
    }

    prefix <- ""
    for (i in c("chr", "Chr")) {
      if (prefix == "") {
        prefix <- getPrefix(i, seql)
      }
    }
    if (prefix != "") {
      if (length(seqlevels(from)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$",
                                       seqlevels(from))]) > 0) {
        if (is(from, "GRanges")) {
          seqlevels(from)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$", seqlevels(from))] <-
            paste(prefix,
                  seqlevels(from)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$",
                                        seqlevels(from))], sep = "")
          #seqlevels(from)[seqlevels(from)=="chrMT"] <- "chrM"
        } else{
          seqnames(from)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$", seqnames(from))] <-
            paste(prefix,
                  seqnames(from)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$",
                                       seqnames(from))], sep = "")
          #seqnames(from)[seqnames(from)=="chrMT"] <- "chrM"
        }
      }
      #if(seqlevelsStyle(from)!="UCSC") seqlevelsStyle(from) <- "UCSC"
      seqlevels(from) <- sub("^chr", prefix, seqlevels(from),
                             ignore.case = TRUE)
    } else{
      ## remove chr
      seqlevels(from) <- sub("^chr", prefix, seqlevels(from),
                             ignore.case = TRUE)
    }
    from
  }
  seql <- seqlevelsStyle(to)
  seqf <- seqlevelsStyle(from)
  if (!seql[1] %in% seqf) {
    tried <- try({
      seqlevelsStyle(from) <- seql[1]
    })
    if (inherits(tried, "try-error")) {
      from <- forceFormatSeqnames(from, to)
    }
  }

  seql <- seqlevelsStyle(to)
  seqf <- seqlevelsStyle(from)
  if (length(intersect(seql, seqf)) == 0) {
    from <- forceFormatSeqnames(from, to)
  }
  from
}

### write chain block
#' @NoRd
.writeInBlock <- function(BW_Score,
                         Rep_Chr,
                         Chr_Size,
                         Rep_Start,
                         Rep_End,
                         Repeat,
                         Rep_Size,
                         Ref_Start,
                         Ref_End,
                         refLine,
                         repLine,
                         reverse,
                         RepeatID,
                         liftIn,
                         include = NULL) {
  if ((as.integer(Ref_End) - as.integer(Ref_Start)) <= 1)
    return()

  if (!is.null(include))
    if (!Repeat %in% include)
      return()

  Rep_End = Rep_Start + stringr::str_count(repLine, "[A-Z]")
  if(reverse==FALSE)
    paste(
      'chain',
      BW_Score,
      Rep_Chr,
      Chr_Size,
      '+',
      format(Rep_Start, scientific = FALSE),
      format(Rep_End,   scientific = FALSE),
      Repeat,
      Rep_Size,
      '+',
      format(Ref_Start, scientific = FALSE),
      format(Ref_End,   scientific = FALSE),
      RepeatID
    ) %>%
    writeLines(con = liftIn)
  else if(reverse==TRUE)
    paste(
      'chain',
      BW_Score,
      Rep_Chr,
      Chr_Size,
      '+',
      format(Rep_Start, scientific = FALSE),
      format(Rep_End,   scientific = FALSE),
      Repeat,
      Rep_Size,
      '-',
      format(Rep_Size-Ref_End, scientific = FALSE),
      format(Rep_Size-Ref_Start,scientific = FALSE),
      RepeatID
    ) %>%
    writeLines(con = liftIn)

  repGaps = stringi::stri_locate_all(repLine, regex = "--*")[[1]]
  refGaps = stringi::stri_locate_all(refLine, regex = "--*")[[1]]

  rbind(cbind(repGaps, -1),
        cbind(refGaps, +1)) -> Gaps

  Gaps[order(Gaps[, 1], na.last = NA),] %>%
    matrix(ncol = 3) -> Merged.Gaps

  if (nrow(Merged.Gaps) > 0) {
    for (i in 1:nrow(Merged.Gaps)) {
      if (i == 1) {
        size = Merged.Gaps[1, 1] - 1
        gap  = Merged.Gaps[1, 2] - Merged.Gaps[1, 1] + 1
      }
      else{
        size = Merged.Gaps[i, 1] - Merged.Gaps[i - 1, 2] - 1
        gap  = Merged.Gaps[i, 2] - Merged.Gaps[i, 1] + 1
      }
      if (Merged.Gaps[i, 3] == 1)
        write.table(list(size, gap, 0),
                    row.names = FALSE, col.names = FALSE,
                    append = TRUE, quote = FALSE,
                    sep = '\t', file = liftIn)
      else
        write.table(list(size, 0, gap),
                    row.names = FALSE, col.names = FALSE,
                    append = TRUE, quote = FALSE,
                    sep = '\t', file = liftIn)
    }
    writeLines(as.character(nchar(repLine) - Merged.Gaps[nrow(Merged.Gaps), 2]), con = liftIn)
  }
  else{
    writeLines(as.character(nchar(repLine)), con = liftIn)
  }

  writeLines("", con = liftIn)
}
