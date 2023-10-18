# 计算相关系数的函数
genecor.parallel <- function(data,gene,cl){
  cl <- makeCluster(cl)
  y <- as.numeric(data[gene,])
  rownames <- rownames(data)
  dataframes <- do.call(rbind, parLapply(cl=cl,rownames, function(x){
    dd  <- cor.test(as.numeric(data[x,]), y, type="spearman")
    data.frame(Gene_1=gene, Gene_2=x, cor=dd$estimate, p.value=dd$p.value)
  }))
  stopCluster(cl)
  return(dataframes)
}

# 画图的函数
genecor_circleplot <- function(x){
  Corr <- data.frame(rbind(data.frame(Gene=x[,1], Correlation=x[,3]), 
                           data.frame(Gene=x[,2], Correlation=x[,3])), stringsAsFactors = F)      
  Corr$Index <- seq(1,nrow(Corr),1) #记录基因的原始排序，记录到Index列
  Corr <- Corr[order(Corr[,1]),] #按照基因名排序
  corrsp <- split(Corr,Corr$Gene)
  corrspe <- lapply(corrsp, function(x){x$Gene_Start<-0
  
  #依次计算每个基因的相关系数总和，作为基因终止位点
  if (nrow(x)==1){x$Gene_End<-1}else{
    x$Gene_End<-sum(abs(x$Correlation))} 
  x})
  GeneID <- do.call(rbind,corrspe)
  GeneID <- GeneID[!duplicated(GeneID$Gene),]
  
  #基因配色
  mycol <- pal_d3("category20c")(20)
  n <- nrow(GeneID)
  GeneID$Color <- mycol[1:n]
  
  #连线的宽度是相关系数的绝对值
  Corr[,2] <- abs(Corr[,2]) 
  corrsl <- split(Corr,Corr$Gene)
  aaaaa <- c()
  corrspl <- lapply(corrsl,function(x){nn<-nrow(x)
  for (i in 1:nn){
    aaaaa[1] <- 0
    aaaaa[i+1] <- x$Correlation[i]+aaaaa[i]}
  bbbbb <- data.frame(V4=aaaaa[1:nn],V5=aaaaa[2:(nn+1)])
  bbbbbb <- cbind(x,bbbbb)
  bbbbbb
  })
  Corr <- do.call(rbind,corrspl)
  
  #根据Index列，把基因恢复到原始排序
  Corr <- Corr[order(Corr$Index),]
  
  #V4是起始位置，V5是终止位置
  #把它写入Links里，start_1和end_1对应Gene_1，start_2和end_2对应Gene_2
  x$start_1 <- Corr$V4[1:(nrow(Corr)/2)]
  x$end_1 <- Corr$V5[1:(nrow(Corr)/2)]
  x$start_2 <- Corr$V4[(nrow(Corr)/2 + 1):nrow(Corr)]
  x$end_2 <- Corr$V5[(nrow(Corr)/2 + 1):nrow(Corr)]
  
  #连线（相关系数）的配色
  #相关系数最大为1，最小-1，此处设置201个颜色
  #-1到0就是前100，0到1就是后100
  color <- data.frame(colorRampPalette(c("#67BE54", "#FFFFFF", "#F82C2B"))(201))
  #根据相关系数的数值，给出相应的颜色
  for (i in 1:nrow(x)){
    x[i,8] <- substring(color[x[i,3] * 100 + 101, 1], 1, 7)
  }
  names(x)[8] <- "color"
  
  #绘图区设置
  #par(mar=rep(0,4))
  circos.clear()
  circos.par(start.degree = 90, #从哪里开始画，沿着逆时针顺序
             gap.degree = 5, #基因bar之间的间隔大小
             track.margin = c(0,0.23), #值越大，基因跟连线的间隔越小
             cell.padding = c(0,0,0,0)
  )
  circos.initialize(factors = GeneID$Gene,
                    xlim = cbind(GeneID$Gene_Start, GeneID$Gene_End))
  
  #先画基因
  circos.trackPlotRegion(ylim = c(0, 1), factors = GeneID$Gene, 
                         track.height = 0.05, #基因线条的胖瘦
                         panel.fun = function(x, y) {
                           name = get.cell.meta.data("sector.index") 
                           i = get.cell.meta.data("sector.numeric.index") 
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(x = mean(xlim), y = 1,
                                       labels = name,
                                       cex = 1, #基因ID文字大小
                                       niceFacing = TRUE, #保持基因名的头朝上
                                       facing = "bending", #基因名沿着圆弧方向，还可以是reverse.clockwise
                                       adj = c(0.5, -2.8), #基因名所在位置，分别控制左右和上下
                                       font = 2 #加粗
                           )
                           circos.rect(xleft = xlim[1], 
                                       ybottom = ylim[1],
                                       xright = xlim[2], 
                                       ytop = ylim[2],
                                       col = GeneID$Color[i],
                                       border = GeneID$Color[i])
                           
                           circos.axis(labels.cex = 0.7, 
                                       direction = "outside"
                           )})
  
  #画连线
  for(i in 1:nrow(x)){
    circos.link(sector.index1 = x$Gene_1[i], 
                point1 = c(x[i, 4], x[i, 5]),
                sector.index2 = x$Gene_2[i], 
                point2 = c(x[i, 6], x[i, 7]),
                col = paste(x$color[i], "C9", sep = ""), 
                border = FALSE, 
                rou = 0.7
    )}
  
  #画图例
  i <- seq(0,0.995,0.005)
  rect(-1+i/2, #xleft
       -1, #ybottom
       -0.9975+i/2, #xright
       -0.96, #ytop
       col = paste(as.character(color[,1]), "FF", sep = ""),
       border = paste(as.character(color[,1]), "FF", sep = ""))
  text(-0.97, -1.03, "-1")
  text(-0.51, -1.03, "1")
}