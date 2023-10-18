library(ggplot2)
library(ggrepel)
##' correlaiton plot for pan tissue analysis
##'
##' @title pancorplot
##' @param font.size font size
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 scale_x_log10
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 theme
##' @importFrom ggrepel geom_text_repel
##' @importFrom ggplot2 aes
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_gpsa
##' @export
##' 
pancorplot <- function(data=plotdf,
                       label=c("all","none","BRCA"),
                       shape = 21,size=5,alpha=0.6,
                       colour="black",stroke = 1.5,
                       labelcolor = "black",
                       box.padding = 0.6, max.overlaps = Inf,
                       legend.position="top"){
  
  data$group = ifelse(data$p.value <= 0.05,ifelse(data$cor>=0,"1","2"),"3")
  data$group = factor(data$group,levels = c("1","2","3"),labels = c("pos","neg","none"))
  
  options(scipen = 2)
  p=ggplot(data,aes(-log10(p.value),cor,fill=group))+
    geom_point(size=size,alpha=alpha,shape =shape,colour=colour,stroke =stroke)+
    scale_y_continuous(expand = c(0,0),limits = c(-1.1,1.1),breaks = seq(-1,1,0.2))+
    scale_x_log10(limits = c(0.01, 1000),breaks = c(0.01,0.1,10,1000))+
    geom_hline(yintercept = 0,size=1.2)+
    geom_vline(xintercept = -log10(0.05),size=1.2)+
    labs(x=bquote(-log[10]~italic("P")),y=paste0(data$method[1]," correlation (r)"))+
    theme(legend.position=legend.position)+
    theme_gpsa()
  
  if(label[1]=="all"){
    label = "all"
    p + geom_text_repel(aes(label=type),col=labelcolor,box.padding = box.padding, max.overlaps = max.overlaps)
  }else if(label[1]=="none"){
    p
  }else if(label[1]=="sig"){
    p + geom_text_repel(data=subset(data, data$p.value < 0.05),
                        aes(label=type),col=labelcolor,box.padding = box.padding, max.overlaps = max.overlaps)
  }else{
    p + geom_text_repel(data=subset(data, data$type %in% label),
                        aes(label=type),col=labelcolor,box.padding = box.padding, max.overlaps = max.overlaps)
  }
}


##' correlaiton plot for pan tissue analysis
##'
##' @title getpancordata
##' @param gene1
##' @param gene2
##' @param data
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_gpsa
##' @export
##' 
getpancordata <- function(gene1,gene2,data=tcga_splitdata,
                          alternative = c("two.sided", "less", "greater"),
                          method = c("pearson", "kendall", "spearman"),
                          exact = NULL, conf.level = 0.95, continuity = FALSE, ...){
  do.call(rbind,lapply(data, function(x){
    dd = cor.test(as.numeric(x[,gene1]),as.numeric(x[,gene2]),
                  alternative = alternative,
                  method = method,
                  exact =exact, conf.level = conf.level, continuity = continuity, ...)
    data.frame(type=x$type[1],cor=dd$estimate,p.value=dd$p.value,method = method)
  }))
}


##' correlaiton plot for two element by guozi
##'
##' @title sincorplot_gz
##' @param data
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_smooth
##' @importFrom ggplot2 geom_rug
##' @importFrom ggplot2 theme_minimal
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 theme
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_gpsa
##' @export
##' 
##' 
sincorplot <- function(a,b,type="ALL",data= tcga_mRNA_exprSet,method = "spearman",title=T,source = "TCGA"){
  if(type=="ALL"){
    plot_df <- data[,c(a,b)]
    type = "Data"
  }else{
    plot_df <- data[data$type %in% type,c(a,b)]
  }
  
  if(title){title = paste0(type," from ",source)}else{title=NULL}
  
  names(plot_df) <- c("geneA","geneB")
  
  ggplot(plot_df,aes(geneA,geneB))+
    geom_point(col="#984ea3")+
    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
    geom_rug(col="#7fc97f")+
    stat_cor(method = method,digits = 3,size=5)+
    xlab(a)+
    ylab(b)+
    labs(x=a,y=b,title = title)+
    theme_bw()+
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))
}

##' ggplot theme of gpsa
##'
##' @title theme_gpsa
##' @param font.size font size
##' @return ggplot theme
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 element_line
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_rect
##' @importFrom ggplot2 margin
##' @importFrom ggplot2 unit
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_gpsa
##' @export
theme_gpsa <- function() {
  theme(axis.title=element_text(size=20),
        axis.text = element_text(face = "bold",size = 16),
        axis.ticks.length=unit(.4, "cm"),
        axis.ticks = element_line(colour = "black", size = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
}



pancorBatch <- function(gene1,gene2,splitdata){
  pancordata = do.call(rbind,lapply(splitdata, function(x){
    dd  <- cor.test(as.numeric(x[,gene1]),as.numeric(x[,gene2]),method="spearman")
    data.frame(type=x$type[1],cor=dd$estimate,p.value=dd$p.value )
  }))
  pancordata = pancordata[!is.na(pancordata$p.value),]
  subdata = pancordata[pancordata$p.value < 0.05,]
  data.frame(gene1=gene1,
             count= nrow(subdata),
             avercor= mean(subdata$cor),
             positive= sum(subdata$cor > 0),
             negative= sum(subdata$cor < 0))
}

