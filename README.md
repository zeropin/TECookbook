# TECookbook: Annotating and parsing transposable-elements associated data

This package is developed to annotate and parse Transposabe-elements (TEs)-associated data.

## Currently it contains four functions:

**annotateSitesInRepeat**: Annotate a list of binding sites that are contained in certain repeat elements, defined by RepeatMasker

**annotatePeaksNearRepeat**: Annotate a list of ChIP-seq peaks that overlap with repeat elements

**buildChain**: Construct a liftOver chain file and repeat sizes file based on [RepeatMasker alignment file (.align)](https://repeatmasker.org/species/hg.html) to map ChIP signals onto repeat coordinates

**liftOut**: Lift all sites out of a specific repeat family at defined locus

## There are some prebuilt datasets that you can download and use in conjuction with this package directly

[**Repeat.Human.hg38**](https://share.weiyun.com/3gXU6Chs): Annotation data based on RepeatMasker file hg38.fa.out

**Repeat.Human.hg19**: Annotation data based on RepeatMasker file hg19.fa.out

**Repeat.Mouse.mm10**: Annotation data based on RepeatMasker file mm10.fa.out

## Installation instruction:

You can pull and install TECookbook package through R command: `remotes::install_github("zeropin/TECookbook")`

## Acknowledgement:

Some functions of this package, such as annotate\*\*, are derived and modified from other packages [ChIPpeakAnno](https://github.com/jianhong/ChIPpeakAnno). I want to thank their generous sharing of the source codes and permission for reuse.

Please contact me if you have any question.

Zheng Zuo [zeropin\@live.cn](mailto:zeropin@live.cn)
