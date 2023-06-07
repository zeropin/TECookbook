# TECookbook: Annotating and parsing transposable-elements associated data

This package is developed to annotate and parse Transposabe-elements (TEs)-associated data. It has been successfully used to study the co-evolution relationship of ZNF10 and ZNF382 with LINE-1 family repeats, which can be accessed via [ZNF10](https://github.com/zeropin/ZFPCookbook/tree/master/ZNF10) and [ZNF382](https://github.com/zeropin/ZFPCookbook/tree/master/ZNF382) respectively.

If you are analyzing ChIP-seq data associated with human or mouse repeats, you can directly download [prebuilt annotation datasets](https://share.weiyun.com/wB9jqSaO) to save time. The easiest way to start using this package is to copy and modify existing workflow of [ZNF10](https://github.com/zeropin/ZFPCookbook/blob/master/ZNF10/htmls/Analysis-of-ZNF10-with-LINE1s.pdf) or [ZNF382](https://github.com/zeropin/ZFPCookbook/blob/master/ZNF382/htmls/Analysis-of-ZNF382-with-LIINE-1.pdf), which demostrates how to process raw ChIP-seq sequencing files into normalized signal tracks defined in reference repeat coordinates (bedgraph formats) for visualization and how to extract the putative binding sites sequence from any defined repeat locus for specificity analysis. Unix-like environment is recommended to run the [Makefile](https://github.com/zeropin/ZFPCookbook/blob/master/ZNF10/R/Makefile) for liftIn operation.

If you are analyzing some non-human/mouse repeats data, you need to download and process the standard RepeatMasker output (.out and .align files) into some chain file using the **buildChain** function of this package first, and then perform liftIn operation for signal visualization.

## Functions and code examples of TECookbook

**buildAnnotation**: Convert [RepeatMasker file (.out)](https://repeatmasker.org/species/hg.html) to Genomic Range for annotation of ChIP-seq peaks or binding sites

```r
Example: build annotation datasets from human or mouse RepeatMasker output
Repeats.Human.hg38 = buildAnnotation("hg38.fa.out")
```

**annotateSitesInRepeat**: Annotate a list of binding sites that are contained in certain repeat elements, defined by RepeatMasker

```r
Example: annotate mySiteList by Repeat.hg38, and output all sites
library(TECookbook)
data(Reoeats.Human.hg38)
annotateSite = annotateSitesInRepeat(mySiteList, AnnotationData=Repeats.Human.hg38, 
               SiteLocForDistance = "middle", RepeatLocForDistance = "middle",
               output = "all")
```
**annotatePeaksNearRepeat**: Annotate a list of ChIP-seq peaks that overlap with repeat elements

```r
Example: annotate myPeakList by Repeat.hg38, and output all peaks
library(TECookbook)
data(Repeats.Human.hg38)
annotatePeak = annotatePeaksNearRepeat(myPeakList, AnnotationData=Repeats.Human.hg38,
                                       minOverlap=1, output = "all")
```

**buildChain**: Construct a liftOver chain file and repeat sizes file based on [RepeatMasker alignment file (.align)](https://repeatmasker.org/species/hg.html) to map ChIP signals onto repeat coordinates

```r
Example: convert human RepeatMasker output into chain file for liftIn operation
TECookook::buildChain(alignment = "hg38.fa.align",
                      chainFile = "Hg38ToRepeat.over.chain",
                      sizeFile = "hg38.Repeat.sizes")
```

**liftOut**: Lift all sites out of a specific repeat family at defined locus based on [RepeatMasker alignment file (.align)](https://repeatmasker.org/species/hg.html)

```r
Example: extract all sequences located in 202 to 221 positions of human THE1B elements
sites = liftOut(alignment = "hg38.fa.align", Repeat = "THE1B", start_pos = 202, end_pos = 221)
```

## Installation instruction:

You can pull and install TECookbook package through R command:
```r
remotes::install_github("zeropin/TECookbook")
```
RStudio is recommended to install and use this package.

## Prebuilt datasets

There are some prebuilt datasets that you can download and use in conjuction with this package directly

[**Repeat.Human.hg38**](https://share.weiyun.com/3gXU6Chs): Annotation data based on RepeatMasker file hg38.fa.out

[**Repeat.Human.hg19**](https://share.weiyun.com/tIlSmg3m): Annotation data based on RepeatMasker file hg19.fa.out

[**Repeat.Mouse.mm10**](https://share.weiyun.com/TIYK2Q8s): Annotation data based on RepeatMasker file mm10.fa.out

[**Hg38ToRepeat.over.chain**](https://share.weiyun.com/wB9jqSaO): The liftOver file from hg38 genomic coordinates to reference repeat coordinates, constructed by buildChain function

[**Repeat.sizes**](https://share.weiyun.com/wB9jqSaO): The Repeat sizes file in Human genome, constructed by buildChain function


## Acknowledgement:

Some functions of this package, such as annotate\*\*, are derived and modified from other packages [ChIPpeakAnno](https://github.com/jianhong/ChIPpeakAnno). I want to thank their generous sharing of the source codes and permission for reuse.

Please contact me if you have any question.

Zheng Zuo [zeropin\@live.cn](mailto:zeropin@live.cn)
