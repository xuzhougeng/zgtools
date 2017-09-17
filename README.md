# 垃圾R包：zgtools使用指北

## 安装方法：
```r
# 这个必须装，感谢Y叔提醒
install.package("devtools")
source("https://bioconductor.org/biocLite.R")
# 感谢Y叔提醒
biocLite("xuzhougeng/zgtools")
```

## 功能：
顾名思义，zgtools, 就是我的常用工具，比如说我会把我自己会经常用到的流程写成函数。目前的功能如下：
- auto_ge_analyzer.R: 自动化基因表达分析流程。 说是自动，其实你还是需要自己提供文件路径和试验设计。其中文件必须得是用[salmon](https://combine-lab.github.io/salmon/)得到的count matrix file。你看这个R包是不是很垃圾，完全不能实现我心里那种，大喊一声"给我分析RNA-Seq数据"， 过了几分钟就有结果的那种效果。
- enrich_analyzer.R: auto_ge_analyzer的其中一个组件，把Y叔的[clusterProfiler](https://github.com/GuangchuangYu/clusterProfiler) 里面的富集分析工具放到一个函数里面。 没有任何的创意，我唯一的工作就是写了一个函数，我自己都看不过去了。
- eda_plot.R : auto_ge_analyzer的其中一个组件。 把一些探索性数据分析(explorary data analysis)的图用一个函数的方式画了出来。 我的工作就是写了一个函数，做了code的搬运工，一点都不geek。
- 其他： 我没有想好写啥，慢慢更新吧。

## 使用方法
> 这个主要是写个自己看的，难得写了多那么函数的R包，而不是一个函数，并且有那么多参数，不写一个使用方法，以后根本用不来，

```r
col_data <- data.frame(sample= paste0('sample', seq(1,12,1)),
                          line=c(103,106,111,113,128,130,131,133,142,148,154,158),
                          types=factor(rep(c('col-0',191),each=3,times=2)),
                          treat_time=factor(rep(c('0h','4h'),each=6)),
                       row.names = paste0('sample', seq(1,12,1)) )
files <- file.path("F:/Data/YBB/counts/salmon_smem/",paste0('sample',col_data$line), "quant.sf")
names(files) <- paste0("sample", 1:12)
# txdb and orgdb
library(AnnotationHub)
suppressMessages(ah <- AnnotationHub())
suppressMessages(txdb <- ah[['AH52247']])
suppressMessages(org <- ah[["AH53758"]])

# analysis
library(zgtools)
dds <- file2dds(files, "salmon", colData = col_data, design = ~ treat_time, txDb = txdb)
enrich_list <- auto_ge_analyzer(filepath = files, type = 'salmon', colData = col_data, design = ~ treat_time, txDb = txdb, orgDb = org)

```


## 版本更新
v0.13: 增加了`file2dds`函数. 目前是根据`salmon`的结果路径和实验设计构建`DESeqDataSeq`对象。以后希望支持多种格式。


## 注意事项
由于目前仅跑过拟南芥，所以仅仅测试过拟南芥，其他物种，我就当不存在。

## 缘起
这是我在学习 Hadley Wickham 的 [R package](http://r-pkgs.had.co.nz/r.html) 后的第一个R包作品。

我一直想把自己常用的代码封装起来，这样就能非常方便的使用了。于是我没事就去Y叔的[GitHub](https://github.com/GuangchuangYu/) 观摩他的代码，然后我总能在他的代码部分看到如下内容
```r
##' drawing phylogenetic tree from phylo object
##'
##'
##' @title ggtree
##' @param tr phylo object
##' @param mapping aes mapping
##' @param layout one of 'rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle' or 'daylight'
##' @param open.angle open angle, only for 'fan' layout
##' @param mrsd most recent sampling date
```

不明觉厉，于是我就默默地去搜索这些神秘代码是什么出处，于是就被我搜到Hadley写的[R package](http://r-pkgs.had.co.nz/r.html) 这本书。最好的学习方式就是践行，于是在稍微看了这本书后，我就动手写这个R包了。


这个R包，不出意外，一定会是一个不太好用的R包，里面充斥了大量的代码问题，比如说没有考虑到系统变量，没有考虑到用户的使用习惯，没有太多的创新性。 但是，至少我把他写出来了，后面就是不断改进了。
