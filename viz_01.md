VIZ\_01
================
Faith Igomodu
June 9, 2024

# Explatory Data Analysis of Gene Expression Dataset

Graphs are used in data visualization of gene expression to uncovers
patterns that answer questons concerning variability in gene expression.
The following graphs are used;

-   Box plot
-   Violin plot
-   Stripchart
-   Dot plot
-   Histogram and density plots
-   ECDF plot
-   Q-Q plot
-   Volcano plot
-   Heatmaps
-   Dendograms

For this practice set, will generate box plot using gene expression data
from the Cancer Genome Atlas Dataset.

## Load Libraries

``` r
library(ggpubr)
```

    ## Loading required package: ggplot2

``` r
library(ggplot2)
```

## The Dataset

The Cancer Genome Atlas data is a R package and is free to the public.
The data has clinical and genomic data from gene expression to CNV
profiling, SNP genotyping, DNA methylation, miRNA profiling, exome
sequencing and more. The RTCGA R package is made by Marcin Kosinski et
al., and is a easy way to work with genomic data made available by The
Cancer Genome Atlas Data.

``` r
#Install RTCGA package
#BiocManager::install("RTCGA")
```

mRNA expression data for five genes GATA3, PTEN, XBP1, ESR1 and MUC1
from three different cancer datasets;

-Breast Invasive Carcinoma -Ovarian Serous Cytadenocarcinoma -Lung
Squamous Cell Carcinoma

will be used for the following analysis.

``` r
#load mRNA expression dataset 
library(RTCGA)
```

    ## Welcome to the RTCGA (version: 1.34.0). Read more about the project under https://rtcga.github.io/RTCGA/

``` r
infoTCGA()
```

    ##                   Cohort  BCR Clinical   CN LowP Methylation mRNA mRNASeq miR
    ## ACC-counts           ACC   92       92   90    0          80    0      79   0
    ## BLCA-counts         BLCA  412      412  410  112         412    0     408   0
    ## BRCA-counts         BRCA 1098     1097 1089   19        1097  526    1093   0
    ## CESC-counts         CESC  307      307  295   50         307    0     304   0
    ## CHOL-counts         CHOL   51       45   36    0          36    0      36   0
    ## COAD-counts         COAD  460      458  451   69         457  153     457   0
    ## COADREAD-counts COADREAD  631      629  616  104         622  222     623   0
    ## DLBC-counts         DLBC   58       48   48    0          48    0      48   0
    ## ESCA-counts         ESCA  185      185  184   51         185    0     184   0
    ## FPPP-counts         FPPP   38       38    0    0           0    0       0   0
    ## GBM-counts           GBM  613      595  577    0         420  540     160 565
    ## GBMLGG-counts     GBMLGG 1129     1110 1090   52         936  567     676 565
    ## HNSC-counts         HNSC  528      528  522  108         528    0     520   0
    ## KICH-counts         KICH  113      113   66    0          66    0      66   0
    ## KIPAN-counts       KIPAN  973      941  883    0         892   88     889   0
    ## KIRC-counts         KIRC  537      537  528    0         535   72     533   0
    ## KIRP-counts         KIRP  323      291  289    0         291   16     290   0
    ## LAML-counts         LAML  200      200  197    0         194    0     179   0
    ## LGG-counts           LGG  516      515  513   52         516   27     516   0
    ## LIHC-counts         LIHC  377      377  370    0         377    0     371   0
    ## LUAD-counts         LUAD  585      522  516  120         578   32     515   0
    ## LUSC-counts         LUSC  504      504  501    0         503  154     501   0
    ## MESO-counts         MESO   87       87   87    0          87    0      87   0
    ## OV-counts             OV  602      591  586    0         594  574     304 570
    ## PAAD-counts         PAAD  185      185  184    0         184    0     178   0
    ## PCPG-counts         PCPG  179      179  175    0         179    0     179   0
    ## PRAD-counts         PRAD  499      499  492  115         498    0     497   0
    ## READ-counts         READ  171      171  165   35         165   69     166   0
    ## SARC-counts         SARC  261      261  257    0         261    0     259   0
    ## SKCM-counts         SKCM  470      470  469  118         470    0     469   0
    ## STAD-counts         STAD  443      443  442  107         443    0     415   0
    ## STES-counts         STES  628      628  626  158         628    0     599   0
    ## TGCT-counts         TGCT  150      134  150    0         150    0     150   0
    ## THCA-counts         THCA  503      503  499   98         503    0     501   0
    ## THYM-counts         THYM  124      124  123    0         124    0     120   0
    ## UCEC-counts         UCEC  560      548  540  106         547   54     545   0
    ## UCS-counts           UCS   57       57   56    0          57    0      57   0
    ## UVM-counts           UVM   80       80   80   51          80    0      80   0
    ##                 miRSeq RPPA MAF rawMAF
    ## ACC-counts          80   46  90      0
    ## BLCA-counts        409  344 130    395
    ## BRCA-counts       1078  887 977      0
    ## CESC-counts        307  173 194      0
    ## CHOL-counts         36   30  35      0
    ## COAD-counts        406  360 154    367
    ## COADREAD-counts    549  491 223    489
    ## DLBC-counts         47   33  48      0
    ## ESCA-counts        184  126 185      0
    ## FPPP-counts         23    0   0      0
    ## GBM-counts           0  238 290    290
    ## GBMLGG-counts      512  668 576    806
    ## HNSC-counts        523  212 279    510
    ## KICH-counts         66   63  66     66
    ## KIPAN-counts       873  756 644    799
    ## KIRC-counts        516  478 417    451
    ## KIRP-counts        291  215 161    282
    ## LAML-counts        188    0 197      0
    ## LGG-counts         512  430 286    516
    ## LIHC-counts        372   63 198    373
    ## LUAD-counts        513  365 230    542
    ## LUSC-counts        478  328 178      0
    ## MESO-counts         87   63   0      0
    ## OV-counts          453  426 316    469
    ## PAAD-counts        178  123 150    184
    ## PCPG-counts        179   80 179      0
    ## PRAD-counts        494  352 332    498
    ## READ-counts        143  131  69    122
    ## SARC-counts        259  223 247      0
    ## SKCM-counts        448  353 343    366
    ## STAD-counts        436  357 289    395
    ## STES-counts        620  483 474    395
    ## TGCT-counts        150  118 149      0
    ## THCA-counts        502  222 402    496
    ## THYM-counts        124   90 123      0
    ## UCEC-counts        538  440 248      0
    ## UCS-counts          56   48  57      0
    ## UVM-counts          80   12  80      0

``` r
# Install mRNA dataset 
#installTCGA(packages = "RTCGA.mRNA" )
```

``` r
#Install clinical dataset 
#installTCGA(packages = "RTCGA.clinical" )
```

``` r
#Extract data from TCGA dataset 
library(RTCGA.mRNA)
library(RTCGA.clinical)
expr_data <- expressionsTCGA(BRCA.mRNA, OV.mRNA, LUSC.mRNA,
                        extract.cols = c("GATA3", "PTEN", "XBP1","ESR1", "MUC1"))
expr_data
```

    ## # A tibble: 1,305 × 7
    ##    bcr_patient_barcode          dataset    GATA3   PTEN   XBP1   ESR1    MUC1
    ##    <chr>                        <chr>      <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ##  1 TCGA-A1-A0SD-01A-11R-A115-07 BRCA.mRNA  2.87   1.36   2.98   3.08   1.65  
    ##  2 TCGA-A1-A0SE-01A-11R-A084-07 BRCA.mRNA  2.17   0.428  2.55   2.39   3.08  
    ##  3 TCGA-A1-A0SH-01A-11R-A084-07 BRCA.mRNA  1.32   1.31   3.02   0.791  2.99  
    ##  4 TCGA-A1-A0SJ-01A-11R-A084-07 BRCA.mRNA  1.84   0.810  3.13   2.50  -1.92  
    ##  5 TCGA-A1-A0SK-01A-12R-A084-07 BRCA.mRNA -6.03   0.251 -1.45  -4.86  -1.17  
    ##  6 TCGA-A1-A0SM-01A-11R-A084-07 BRCA.mRNA  1.80   1.31   4.04   2.80   3.53  
    ##  7 TCGA-A1-A0SO-01A-22R-A084-07 BRCA.mRNA -4.88  -0.237 -0.725 -4.49  -1.46  
    ##  8 TCGA-A1-A0SP-01A-11R-A084-07 BRCA.mRNA -3.14  -1.24  -1.19  -1.63  -0.987 
    ##  9 TCGA-A2-A04N-01A-11R-A115-07 BRCA.mRNA  2.03   1.21   2.28   4.12   0.668 
    ## 10 TCGA-A2-A04P-01A-31R-A034-07 BRCA.mRNA -0.293  0.288 -1.61   0.473  0.0115
    ## # ℹ 1,295 more rows

``` r
#Determine number of samples in each dataset 

samples_num <- table(expr_data$dataset)
samples_num
```

    ## 
    ## BRCA.mRNA LUSC.mRNA   OV.mRNA 
    ##       590       154       561

``` r
#Simplify name of dataset by removing mRNA 
expr_data$dataset <- gsub(pattern = ".mRNA", replacement = "",  expr_data$dataset)
```

``` r
#Simplify patient barcode column 
expr_data$bcr_patient_barcode <- paste0(expr_data$dataset, c(1:590, 1:561, 1:154))
expr_data
```

    ## # A tibble: 1,305 × 7
    ##    bcr_patient_barcode dataset  GATA3   PTEN   XBP1   ESR1    MUC1
    ##    <chr>               <chr>    <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ##  1 BRCA1               BRCA     2.87   1.36   2.98   3.08   1.65  
    ##  2 BRCA2               BRCA     2.17   0.428  2.55   2.39   3.08  
    ##  3 BRCA3               BRCA     1.32   1.31   3.02   0.791  2.99  
    ##  4 BRCA4               BRCA     1.84   0.810  3.13   2.50  -1.92  
    ##  5 BRCA5               BRCA    -6.03   0.251 -1.45  -4.86  -1.17  
    ##  6 BRCA6               BRCA     1.80   1.31   4.04   2.80   3.53  
    ##  7 BRCA7               BRCA    -4.88  -0.237 -0.725 -4.49  -1.46  
    ##  8 BRCA8               BRCA    -3.14  -1.24  -1.19  -1.63  -0.987 
    ##  9 BRCA9               BRCA     2.03   1.21   2.28   4.12   0.668 
    ## 10 BRCA10              BRCA    -0.293  0.288 -1.61   0.473  0.0115
    ## # ℹ 1,295 more rows

Genes GATA3, PTEN, and XBP1 are cancer genes and are used during this
tutorial. Moreover,

-   GATA3 - a bio-marker of lunminal tumors, a type of breast cancer
    that start in the lining of the mammary ducts. The expectation is
    higher expression of GATA3 genens in BRCA gene expression data than
    in OV and LUSC.

-   PTEN (Phosphatase and tensin homolog) encodes a phosphatase enzyme
    in humans and is a tumor supressor gene. The loss and or varation in
    the gene expression level of PTEN is associated with a broad
    spectrum of human cancers. The expectation is PTEN will be expressed
    at different rates in all three tissue ; breast, lung and ovaries.

-   XBP1 (X-box binding protein 1) - encodes the XBP1 transcription
    factor and is upregulated in cancer cells and promotes tumor growth
    and metastasis. The expression of XBP1 should be high in all three
    datasets.

The objective is to compare the gene expression level of PTEN, XBP1, and
GATA3 in each cancer cell type.Then compare expression level in
different Cell subpopulations.

``` r
#Plot the GATA3, PTEN, and XBP1 without creating a list. 

library(ggpubr)

# GATA3
ggboxplot(expr_data, x = "dataset", y = "GATA3",
          title = "GATA3", ylab = "Expression",
          color = "dataset", palette = "lancet")
```

![](viz_01_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# PTEN
ggboxplot(expr_data, x = "dataset", y = "PTEN",
          title = "PTEN", ylab = "Expression",
          color = "dataset", palette = "lancet")
```

![](viz_01_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
# XBP1
ggboxplot(expr_data, x = "dataset", y = "XBP1",
          title = "XBP1", ylab = "Expression",
          color = "dataset", palette = "lancet")
```

![](viz_01_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
#Plot the GATA3, PTEN, and XBP1 by creating a list. 


list_ofplots <- ggboxplot(expr_data, x = "dataset", 
               y = c("GATA3", "PTEN", "XBP1"),
               title = c("GATA3", "PTEN", "XBP1"),
               ylab = "Expression", 
               color = "dataset", palette = "lancet")
# View GATA3
list_ofplots$GATA3
```

![](viz_01_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# View PTEN
list_ofplots$PTEN
```

![](viz_01_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
# View XBP1
list_ofplots$XBP1
```

![](viz_01_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
#Combine all three plots onto one planel. 

ggboxplot(expr_data, x = "dataset",
          y = c("GATA3", "PTEN", "XBP1"),
          combine = TRUE,
          ylab = "Expression",
          color = "dataset", palette = "lancet")
```

![](viz_01_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

From graph below, GATA3 and XBP1 expression level is higher in BRCA
dataset than in OV, and LUSC. This is consistent with the hypothesis
mentioned earlier.PTEN expression level is about the same in all three
datasets.Next compare the gene expression of cell subpopulations.

``` r
#Merge the plots 

ggboxplot(expr_data, x = "dataset",
          y = c("GATA3", "PTEN", "XBP1"),
          merge = TRUE,
          ylab = "Expression", 
          palette = "jco")
```

![](viz_01_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

The plot below shows difference in gene expression among cell types.
From graph, BRCA, OV and LUSC cells each have a unique spread in gene
expression.

``` r
#Generate  plot that compares expression level in different cell subpopulations 

ggboxplot(expr_data, x = "dataset",
          y = c("GATA3", "PTEN", "XBP1"),
          merge = "flip",
          ylab = "Expression", 
          palette = "jco")
```

![](viz_01_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Jitter is a scatterplot laid over box plot to identify genes with
highest and lowest gene expression values. Each jitter point is a
datapoint in the dataset.

``` r
ggboxplot(expr_data, x = "dataset",
          y = c("GATA3", "PTEN", "XBP1"),
          combine = TRUE,
          color = "dataset", palette = "jco",
          ylab = "Expression", 
          # Add jittered points
          add = "jitter",  
          # Point size and the amount of jittering
          add.params = list(size = 0.1, jitter = 0.2), 
          # column with point labels
          label = "bcr_patient_barcode",                
           # Select  labels to display
          label.select = list(top.up = 2, top.down = 2),
         # label font
          font.label = list(size = 8, face = "italic"), 
           # Stop label text overplotting
          repel = TRUE                                 
          
          )
```

![](viz_01_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## References

Facilitating exploratory data visualization: application to TCGA genomic
Data - Articles - STHDA. (n.d.).
<http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/77-facilitating-exploratory-data-visualization-application-to-tcga-genomic-data/>
