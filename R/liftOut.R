#' Function to extract all genomic regions that maps to some specific locus within a reference repeat element
#'
#' @param alignment RepeatMasker alignment file used for extracting genomic region, e.g., hg38.fa.align
#' @param Repeat Name of Repeat to be extracted
#' @param start_pos Starting position of extracted locus within repeat consensus sequence, defined in Dfam database
#' @param end_pos   End position of extracted locus within repeat consensus sequence, defined in Dfam database
#' @param RepeatTagID  Extract regions only within defined RepeatTagIDs defined in alignment file, e.g., c_b1s251i0; If NULL, no such restriction.
#' @param RepeatID  Extract regions only within defined RepeatIDs, e.g., 1, 2, 3, etc; If NULL, no such restriction.
#' @examples
#'          ##Example 1: lift out all sites from THE1B repeats in hg38 genome
#'          sites = liftOut(alignment = "hg38.fa.align", Repeat = "THE1B", start_pos = 202, end_pos = 221)
#'          ##Example 2: lift out sites from L1P1_orf2 repeats with repeatID 1000 and 1005
#'          sites = liftOUt(alignment = "hg38.fa.align", Repeat = "L1P1_orf2", start_pos = 1000, end_pos = 1020, RepeatID = c(1000, 1005))
#' @export
liftOut <- function(alignment,
                    Repeat,
                    start_pos,
                    end_pos,
                    RepeatTagID = NULL,
                    RepeatID = NULL) {

  Target_Repeat   = Repeat
  Target_RepeatTagID = unique(RepeatTagID)
  Target_RepeatID = unique(RepeatID)

  infCon   <- file(alignment, open = "rt")

  result = tribble( ~seqnames, ~start, ~end, ~strand, ~Sequence, ~RepeatTagID, ~RepeatID)

  #The guideLine is used to find the blocks containing the specified repeat elements
  guideLine = paste0("[0-9].*", Target_Repeat, "#")

  remaining=TRUE
  while (remaining) {
    lines = readLines(infCon, n = 1000000)

    if (length(lines) < 1000000) remaining=FALSE

    while (!startsWith(lines[[length(lines) - 1]], "Gap_init") &
           remaining)
      lines = c(lines, readLines(infCon, n = 1))

    Block_positions = which(stringr::str_starts(lines, guideLine))

    Block_count = 1
    while (Block_count <= length(Block_positions)) {
      line = lines[Block_positions[Block_count]]

      words = strsplit(line, split = " |#",)[[1]]

      BW_Score = words[[1]]
      Rep_Chr  = words[[5]]
      Rep_Start = as.integer(words[[6]]) - 1
      Rep_End  = words[[7]]
      Chr_Left = strsplit(words[[8]], split = "\\(|\\)")[[1]][[2]]

      Chr_Size = as.integer(Rep_End) + as.integer(Chr_Left)

      if (words[[9]] == "C") {
        reverse  = TRUE
        Repeat   = words[[10]]
        Ref_Left = strsplit(words[[12]], split = "\\(|\\)")[[1]][[2]]
        Ref_End  = as.integer(words[[13]])
        Ref_Start = as.integer(words[[14]]) - 1
        RepeatTagID = words[[15]]
        RepeatID = words[[16]]
      }
      else{
        reverse  = FALSE
        Repeat   = words[[9]]
        Ref_Start = as.integer(words[[11]]) - 1
        Ref_End  = as.integer(words[[12]])
        Ref_Left = strsplit(words[[13]], split = "\\(|\\)")[[1]][[2]]
        RepeatTagID = words[[14]]
        RepeatID = words[[15]]
      }

      Rep_Size = as.integer(Ref_End) + as.integer(Ref_Left)
      refLine  = ""
      repLine  = ""

      i = Block_positions[Block_count] + 2
      while (!startsWith(lines[i], "Matrix")) {
        repLine = paste0(repLine, strsplit(lines[i], "  *")[[1]][[4]])
        refLine = paste0(refLine, strsplit(lines[i + 2], "  *")[[1]][[4]])
        i = i + 4
      }

      if (Ref_Start < start_pos &
          Ref_End >= end_pos &
          (is.null(Target_RepeatID) | (RepeatID%in%Target_RepeatID)) &
          (is.null(Target_RepeatTagID) | (RepeatTagID%in%Target_RepeatTagID))) {
        bases_position = stringr::str_locate_all(refLine, "[A-Z]")[[1]][, 'start']

        if (reverse == TRUE) {
          substr(repLine,
                 bases_position[[Ref_End - end_pos + 1]],
                 bases_position[[Ref_End - start_pos + 1]]) %>%
            Biostrings::DNAString() %>%
            Biostrings::reverseComplement() %>%
            as.character() -> Sequence

          stringr::str_count(substr(repLine, 1, bases_position[[Ref_End -
                                                                  end_pos + 1]]),
                             "[A-Z]") + Rep_Start -> start

          stringr::str_count(substr(repLine, 1, bases_position[[Ref_End -
                                                                  start_pos + 1]]),
                             "[A-Z]") + Rep_Start -> end

        }
        else if (!reverse) {
          substr(repLine,
                 bases_position[[start_pos - Ref_Start]],
                 bases_position[[end_pos - Ref_Start]]) -> Sequence

          stringr::str_count(substr(repLine, 1, bases_position[[start_pos -
                                                                  Ref_Start]]),
                             "[A-Z]") + Rep_Start -> start

          stringr::str_count(substr(repLine, 1, bases_position[[end_pos -
                                                                  Ref_Start]]),
                             "[A-Z]") + Rep_Start -> end
        }

        writeLines(paste(Rep_Chr, start, end, Sequence, RepeatTagID, RepeatID))
        result = add_row(
          result,
          seqnames = Rep_Chr,
          start = start,
          end = end,
          strand = ifelse(reverse, "-", "+"),
          Sequence = Sequence,
          RepeatTagID = RepeatTagID,
          RepeatID = RepeatID
        )
      }
      Block_count = Block_count+1
    }
  }

  close(infCon)
  return(GenomicRanges::makeGRangesFromDataFrame(result, keep.extra.columns = TRUE))
}
