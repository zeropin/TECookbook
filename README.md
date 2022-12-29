# TECookbook: Annotating and parsing transposable-elements associated data

This package is developed to annotate and parse Transposabe-elements (TEs)-associated data.

## Functions within TECookbook

**buildAnnotation**: Convert [RepeatMasker file (.out)](https://repeatmasker.org/species/hg.html) to Genomic Range for annotation of ChIP-seq peaks or binding sites

**annotateSitesInRepeat**: Annotate a list of binding sites that are contained in certain repeat elements, defined by RepeatMasker

**annotatePeaksNearRepeat**: Annotate a list of ChIP-seq peaks that overlap with repeat elements

**buildChain**: Construct a liftOver chain file and repeat sizes file based on [RepeatMasker alignment file (.align)](https://repeatmasker.org/species/hg.html) to map ChIP signals onto repeat coordinates

**liftOut**: Lift all sites out of a specific repeat family at defined locus

## Prebuilt datasets

There are some prebuilt datasets that you can download and use in conjuction with this package directly

[**Repeat.Human.hg38**](https://share.weiyun.com/3gXU6Chs): Annotation data based on RepeatMasker file hg38.fa.out

[**Repeat.Human.hg19**](https://share.weiyun.com/tIlSmg3m): Annotation data based on RepeatMasker file hg19.fa.out

[**Repeat.Mouse.mm10**](https://share.weiyun.com/TIYK2Q8s): Annotation data based on RepeatMasker file mm10.fa.out

[**Hg38ToRepeat.over.chain**](https://share.weiyun.com/wB9jqSaO): The liftOver file from hg38 genomic coordinates to reference repeat coordinates, constructed by buildChain function

[**Repeat.sizes**](https://share.weiyun.com/wB9jqSaO): The Repeat sizes file in Human genome, constructed by buildChain function

## Installation instruction:

You can pull and install TECookbook package through R command:
```r
remotes::install_github("zeropin/TECookbook")
```

## Application examples

TECookbook has been successfully used to study the relationship of ZNF10 and ZNF382 with LINE-1 family repeats, which can be accessed via
[ZNF10](https://github.com/zeropin/ZFPCookbook/ZNF10) and [ZNF382](https://github.com/zeropin/ZFPCookbook/ZNF382) respectively.

## Acknowledgement:

Some functions of this package, such as annotate\*\*, are derived and modified from other packages [ChIPpeakAnno](https://github.com/jianhong/ChIPpeakAnno). I want to thank their generous sharing of the source codes and permission for reuse.

Please contact me if you have any question.

Zheng Zuo [zeropin\@live.cn](mailto:zeropin@live.cn)
