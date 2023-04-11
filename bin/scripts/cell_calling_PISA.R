#cell_calling, method in knee and emptydrops

parser = argparse::ArgumentParser(description="")
parser$add_argument('--matrix', help='MEX matrix dir')
parser$add_argument('--forcecells',help='Force pipeline to use this number of cells, bypassing cell calling algorithm.')
parser$add_argument('--expectcells', help='Expected number of recovered cells, used as input to cell calling algorithm. [default: 3000]')
parser$add_argument('--method', help='Cell calling method, Choose from auto and emptydrops, default: auto')
parser$add_argument('--output',help='Output dir, default: current dir')
parser$add_argument('--minumi',help='min umi, default: 1000')
args = parser$parse_args()

if (is.null(args$output)) {
    args$output<-getwd()
}
if (is.null(args$method)) {
    args$method<- "auto"
}
if (is.null(args$expectcells)) {
    args$expectcells<- 3000
}
if (is.null(args$minumi)) {
    args$minumi<- 1000
}

ReadPISA <- function(mex_dir=NULL,
                     barcode.path = NULL,
                     feature.path = NULL,
                     matrix.path=NULL,
                     use_10X=FALSE) {
  if (is.null(mex_dir) && is.null(barcode.path)  && is.null(feature.path) &&
      is.null(matrix.path)) {
    stop("No matrix set.")
  }
  if (!is.null(mex_dir) && !file.exists(mex_dir) ) {
    stop(paste0(mex_dir, " does not exist."))
  }
  if (is.null(barcode.path)  && is.null(feature.path) && is.null(matrix.path)) {
    barcode.path <- paste0(mex_dir, "/barcodes.tsv.gz")
    feature.path <- paste0(mex_dir, "/features.tsv.gz")
    matrix.path <- paste0(mex_dir, "/matrix.mtx.gz")
  }
  spliced.path <- paste0(mex_dir, "/spliced.mtx.gz")
  unspliced.path <- paste0(mex_dir, "/unspliced.mtx.gz")
  spanning.path <- paste0(mex_dir, "/spanning.mtx.gz")
  
  if (!file.exists(barcode.path) || !file.exists(feature.path)) {
    stop(paste0("No expression file found at ", mex_dir))
  }
  
  .ReadPISA0 <- function(barcode.path, feature.path, matrix.path, use_10X) {
    mat <- Matrix::readMM(file = matrix.path)
    feature.names <- read.delim(feature.path,
                                header = FALSE,
                                stringsAsFactors = FALSE
    )
    barcode.names <- read.delim(barcode.path,
                                header = FALSE,
                                stringsAsFactors = FALSE
    )
    colnames(mat) <- barcode.names$V1
    if (use_10X == TRUE) {
      rownames(mat) <- make.unique(feature.names$V2)
    } else {
      rownames(mat) <- make.unique(feature.names$V1)
    }
    mat
  }
  
  if (!file.exists(spliced.path) && file.exists(matrix.path)) {
    return(.ReadPISA0(barcode.path, feature.path, matrix.path, use_10X))
  }
  mat <- list()
  cat("Load spliced matrix ...\n")
  mat$spliced <- Matrix::readMM(file = spliced.path)
  cat("Load unspliced matrix ...\n")
  mat$unspliced <- Matrix::readMM(file = unspliced.path)
  cat("Load spanning matrix ...\n")
  mat$spanning <- Matrix::readMM(file = spanning.path)
  
  feature.names <- read.delim(feature.path,
                              header = FALSE,
                              stringsAsFactors = FALSE
  )
  barcode.names <- read.delim(barcode.path,
                              header = FALSE,
                              stringsAsFactors = FALSE
  )
  colnames(mat$spliced) <- barcode.names$V1
  rownames(mat$spliced) <- make.unique(feature.names$V1)
  colnames(mat$unspliced) <- barcode.names$V1
  rownames(mat$unspliced) <- make.unique(feature.names$V1)
  colnames(mat$spanning) <- barcode.names$V1
  rownames(mat$spanning) <- make.unique(feature.names$V1)
  mat
}

library(DropletUtils)
library(dplyr)
library(ggplot2)
library(cowplot)

mtx <- ReadPISA(args$matrix)
br.out <- barcodeRanks(mtx, lower = 100)
br.out.sort <- br.out[order(br.out$total, decreasing=T),]
len <- nrow(br.out.sort)
UMIsor <- br.out.sort$total

if(!is.null(args$forcecells) & length(numeric(args$forcecells))){
  forcecell <- as.numeric(args$forcecells)
  cutoff <- forcecell
  CellNum <- cutoff
  tmp<-data.frame(barcodes=1:len,UMI=UMIsor,Beads=c(rep("true",cutoff),rep("noise",len-cutoff)))
  beads_barcodes <- rownames(br.out.sort)[1:cutoff]
}else if(args$method=="auto"){
  cutoff <- nrow(subset(br.out.sort,br.out.sort$total>=metadata(br.out)$inflection))
  CellNum <- cutoff
  tmp<-data.frame(barcodes=1:len,UMI=UMIsor,Beads=c(rep("true",cutoff),rep("noise",len-cutoff)))
  beads_barcodes <- rownames(br.out.sort)[1:cutoff]
}else{
  set.seed(123)
  e.out <- emptyDropsCellRanger(mtx,n.expected.cells=as.numeric(args$expectcells), 
                    max.percentile=0.99, max.min.ratio=10,
                    umi.min=as.numeric(args$minumi),umi.min.frac.median=0.01,cand.max.n=20000,
                    ind.min=45000,ind.max=90000,round=TRUE,niters=10000)
  is.cell <- e.out$FDR <= 0.01
  CellNum <- sum(is.cell, na.rm=TRUE)
  cell.counts  <- mtx[,which (is.cell), drop = FALSE ]
  list <- colnames(cell.counts)
  br.out.df <- as.data.frame(br.out.sort)
  br.out.df$barcode <- rownames(br.out.df)
  br.out.df$Beads <- br.out.df$barcode %in% list
  rownames(br.out.df)<-NULL
  br.out.df <- br.out.df%>%arrange(desc(total))%>%mutate(ranks=rownames(br.out.df))
  br.out.df$Beads[which(br.out.df$Beads =='FALSE')] <- 'noise'
  br.out.df$Beads[which(br.out.df$Beads =='TRUE')] <- 'true'
  tmp <- data.frame(br.out.df$ranks,br.out.df$total,br.out.df$Beads)
  colnames(tmp) <- c('barcodes','UMI','Beads')
  beads_barcodes <- subset(br.out.df,Beads=="true")$barcode

}

write.csv(tmp,file = paste(args$output,"/cutoff.csv",sep=""), quote = FALSE, row.names = FALSE)
write.table(beads_barcodes,file = paste(args$output,"/beads_barcodes.txt",sep=""), quote = FALSE, row.names=FALSE, col.names=FALSE)

called <- sum(tmp[which(tmp$Beads %in% "true"),"UMI"])
in_cell <- round(called/sum(tmp$UMI),digits=4)
UMI_mean <-mean(tmp[which(tmp$Beads %in% "true"),"UMI"])
Gene_mean <- mean(colSums(mtx[,beads_barcodes] != 0))




png(file=paste(args$output,"/beads_count_summary.png",sep=""), width=1200,height=400,res=80)
p1 = ggplot(tmp,aes(x=barcodes,y=UMI)) + xlim(0,10) + ylim(0,9)
p1 = p1 +annotate("text",x=0.2,y=9,label="Estimated Number of beads:",size=10,hjust=0)
p1 = p1 +annotate("text",x=3,y=6.5,label=CellNum,size=20,hjust=0, colour="red")
p1 = p1 +annotate("text",x=0.2,y=3.5,label="Reads in beads:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=3.5,label=round(in_cell,3),size=6,hjust=1)
p1 = p1 +annotate("segment",x = 0, xend = 10, y = 3, yend = 3,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=2.5,label="Mean UMI counts per beads:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=2.5,label=round(UMI_mean),size=6,hjust=1)
p1 = p1 +annotate("segment", x = 0, xend = 10, y = 2, yend = 2,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=1.5,label="Mean Genes per beads:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=1.5,label=round(Gene_mean),size=6,hjust=1)
p1 = p1 +annotate("segment",x  = 0, xend  = 10, y  = 1, yend = 1,colour = "blue")
p1 = p1 +theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major=element_blank(),panel.grid.minor = element_blank())

p = ggplot(tmp,aes(x=as.numeric(barcodes),y=as.numeric(UMI)))
p = p + geom_line(aes(color=Beads),size=2) +scale_color_manual(values=c("#999999","blue"))
p = p + scale_x_log10(name="Barcodes",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + scale_y_log10(name="UB",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + theme_bw() + geom_vline(xintercept =CellNum)
p = p + theme(legend.position = "none")


plot_grid(p1, p,  rel_widths = c(2,1.5),ncol=2)

dev.off()
