#' Construct Chain file and Repeat sizes file for liftOver operation between genome coordinates and repeat coordinates
#'
#' @param alignment RepeatMasker alignment file used for extracting genomic region, e.g., hg38.fa.align
#' @param chainFile Name of liftOver chain file to be written
#' @param sizeFile Name of repeat sizes file to be written
#' @param include whether or not only build chain file for select repeat family; if NULL, output all repeat instances;
#' @examples
#'          ##Example 1: Construct chain file for all repeats in hg38 genome
#'          buildChain(alignment = "hg38.fa.align", chainFile = "Hg38ToRepeat.over.chain", sizeFile = "hg38.Repeat.sizes")
#'          ##Example 2: Construct chain file for THE1A repeats only in hg38 genome
#'          buildChain(alignment = "hg38.fa.align", chainFile = "Hg38ToTHE1A.over.chain", sizeFile = "THE1A.Repeat.sizes", include="THE1A")
#' @export
buildChain <- function(alignment   = "hg38.fa.align",
                       chainFile   = "Hg38ToRepeat.over.chain",
                       sizeFile    = "Repeat.sizes",
                       include     = NULL) {
  infCon   <- file(alignment, open = "rt")
  liftIn.Con  <- file(chainFile, open = "wt")

  sizeCon  <- file(sizeFile, open = "wt")
  sizes = list()

  while (TRUE) {
    lines = readLines(infCon, n = 5000000)

    if (length(lines) == 0)
      break

    while (!startsWith(lines[[length(lines) - 2]], "Gap_init"))
      lines = c(lines, readLines(infCon, n = 1))

    i = 1
    while (i <= length(lines)) {
      line = lines[[i]]

      if (grepl("^[0-9]", line)) {
        if (i > 1) {
          sizes[[Repeat]] = Rep_Size
          .writeInBlock(
            BW_Score = BW_Score,
            Rep_Chr = Rep_Chr,
            Chr_Size = Chr_Size,
            Rep_Start = Rep_Start,
            Rep_End = Rep_End,
            Repeat  = Repeat,
            Rep_Size = Rep_Size,
            Ref_Start = Ref_Start,
            Ref_End = Ref_End,
            refLine = refLine,
            repLine = repLine,
            reverse = reverse,
            RepeatID = RepeatID,
            liftIn  = liftIn.Con,
            include  = include
          )
        }

        words = strsplit(line, split = " |#",)[[1]]

        BW_Score = words[[1]]
        Rep_Chr = words[[5]]
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
          RepeatID = words[[16]]
        }
        else{
          reverse  = FALSE
          Repeat   = words[[9]]
          Ref_Start = as.integer(words[[11]]) - 1
          Ref_End  = as.integer(words[[12]])
          Ref_Left = strsplit(words[[13]], split = "\\(|\\)")[[1]][[2]]
          RepeatID = words[[15]]
        }

        Rep_Size = as.integer(Ref_End) + as.integer(Ref_Left)
        refLine  = ""
        repLine  = ""

        i = i + 2
      }
      else if (grepl("^  chr", line)) {
        repLine = paste0(repLine, strsplit(line, "  *")[[1]][[4]])
        refLine = paste0(refLine, strsplit(lines[i + 2], "  *")[[1]][[4]])
        i = i + 2
      }
      else
        i = i + 1

    }

  }

  close(liftIn.Con)

  for (name in names(sizes))
    write.table(list(name, sizes[[name]]),
                row.names = FALSE, col.names = FALSE,
                append = TRUE, quote = FALSE,
                sep = '\t', file = sizeCon)

  close(sizeCon)
}
