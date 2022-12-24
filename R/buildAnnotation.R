#' Function to convert RepeatMasker data to Genomic Range for annotating ChIP-seq peaks or binding sites
#'
#' @param file RepeatMasker file used for extracting repeat name and genomic region, e.g., hg38.fa.out
#' @examples
#'          ##Example:
#'          Repeats.Human.hg38 = buildAnnotation("hg38.fa.out")
#' @export
buildAnnotation <- function(file){
  Repeats <- readr::read_table(file, skip = 2,
                               col_names = c("swScore",	"percDiv",	"percDel",	"percIns",	"seqnames",	"Start",	"End",
                                             "genoLeft",	"strand",	"repName", "repFamily",	"repStart",	"repEnd",	"repLeft",	"id"),
                               col_types = cols(swScore = col_integer(),  percDiv = col_double(), percDel = col_double(), percIns = col_double())) %>%
    dplyr::select(seqnames, Start, End, strand, SW.Score = swScore,
                  Family = repFamily, Name = repName,
                  repeatStart = repStart, repeatEnd = repEnd, repeatLeft = repLeft, repeatID = id
    ) %>%
    dplyr::mutate(Family= as.factor(Family),
                  Name  = as.factor(Name),
                  strand= if_else(strand=="C", "-", "+")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  names(Repeats) <- Repeats$repeatID
  return(Repeats)
}
