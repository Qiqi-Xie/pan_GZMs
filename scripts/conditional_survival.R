conditional_surv_est <- function(basekm, t1, t2) {
  if (class(basekm) != "survfit") {
    stop(
      "Argument to basekm must be of class survfit"
    )
  }
  
  if (max(t1) > max(basekm$time)) {
    stop(
      paste(
        "Argument to t1 specifies a value outside the range of observed times;", "the maximum observed time is", round(max(basekm$time), 2)
      )
    )
  }
  
  if (max(t2) > max(basekm$time)) {
    stop(paste(
      "Argument to t2 specifies a value outside the range of observed times;",
      "the maximum observed time is", round(max(basekm$time), 2)
    ))
  }
  
  cs <- summary(basekm, times = c(t1, t2))$surv[2] /
    summary(basekm, times = c(t1, t2))$surv[1]
  cs.sq <- cs^2
  
  d <- basekm$n.event[basekm$time >= t1 &
                        basekm$time <= t2 &
                        basekm$n.event > 0]
  
  r <- basekm$n.risk[basekm$time >= t1 &
                       basekm$time <= t2 &
                       basekm$n.event > 0]
  
  dr <- d / (r * (r - d))
  var.cs <- 1 / (log(cs)^2) * sum(dr)
  ci <- cs^(exp(c(1, -1) * stats::qnorm(0.975) * sqrt(var.cs)))
  ci.cs <- round(ci, 2)
  
  return(
    list(
      cs_est = round(cs, 2),
      cs_lci = ci.cs[1],
      cs_uci = ci.cs[2]
    )
  )
}

con_prob <- function(i){
  y <- as.numeric(i)
  table=purrr::map_df(
    prob_times,
    ~conditional_surv_est(
      basekm = myfit,
      t1 = y,
      t2 = .x)) %>%
    mutate(years = prob_times) %>%
    select(years, everything())
  table$con=i
  table$label=ifelse(table$years<y,NA,paste0((table$cs_est)*100,"%"))
  data.frame(table)
}

plot_condsurv <- function (survdt = NULL,
                           at, 
                           surv.cut = NULL,
                           main = "", 
                           xlab = "Time (Years)", 
                           ylab = "Survival probability", 
                           curv.col = NULL,
                           lwd = 1.2,
                           legend.pos = "top") 
{
  # survdat: 生存数据，至少包含两列，生存时间（futime，年为单位）和生存状态（fustat）
  # at：一个数值序列，表示条件生存时间中的附加时间，如在第三年存活的基础上额外存活at年的概率
  # surv.cut：KM曲线在surv.cut处截断
  # main：图像主题，默认无
  # xlab：图像x轴名称，默认以年为单位
  # ylab：图像y轴名称，默认以25%为间隔
  # curv.col：各条件生存曲线的颜色
  # lwd：线条宽度
  # legend.pos：图例位置，默认为顶部
  
  library(ggplot2)
  library(survival)
  library(tidyverse)
  
  basekm <- survfit(Surv(futime, fustat)~ 1, data=survdt, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
  
  if (max(at) > max(basekm$time)) {
    stop(paste("Argument to at specifies value(s) outside the range of observed times;", 
               "the maximum observed time is", round(max(basekm$time), 2)))
  }
  if(surv.cut <= max(at)) {
    stop("Time cutoff should be greater than the maximal conditional survival point\n")
  }
  if(!all(is.element(c("fustat","futime"),colnames(survdt)))) {
    stop("Make sure the survival data has two columns of fustat and futime\n")
  }
  
  nt <- length(at)
  fitkm <- list()
  fitkmdat <- list()
  for (i in 1:nt) {
    fitkm[[i]] <- survival::survfit(formula = stats::as.formula(basekm$call$formula), 
                                    data = eval(basekm$call$data), start.time = at[i])
    fitkmdat[[i]] <- tibble::tibble(timept = fitkm[[i]]$time, 
                                    prob = fitkm[[i]]$surv)
  }
  condsurvdat <- fitkmdat %>% purrr::map_df(`[`, .id = "which_at") %>% 
    dplyr::mutate(condtime = factor(which_at, levels = seq(1, nt), labels = at))
  
  condsurvdat$condtime <- paste0(condsurvdat$condtime," Year")
  
  # 第一幅条件生存概率曲线
  if(is.null(surv.cut)) {
    condsurv_plot <- ggplot(condsurvdat, aes(x = timept, y = prob, color = condtime)) + 
      scale_color_manual(values = mycol) + 
      geom_step(lwd = lwd) + 
      ylim(0, 1) + 
      labs(x = xlab, y = ylab, title = main, color = "Given conditional survival") + 
      geom_segment(aes(x = 0, y = 1, xend = at[length(at)], yend = 1),lwd = lwd, color = curv.col[at[length(at)] + 1]) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
            panel.border = element_rect(colour = "black", fill=NA, size=0.2),
            axis.ticks = element_line(size=0.2, color="black"),
            axis.ticks.length = unit(0.2, "cm"),
            legend.position = legend.pos,
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 10))
    ggsave("conditional_survival_curve.pdf")
    
  } else {
    condsurv_plot <- ggplot(condsurvdat, aes(x = timept, y = prob, color = condtime)) + 
      scale_color_manual(values = mycol) + 
      geom_step(lwd = lwd) + 
      ylim(0, 1) + 
      xlim(0, surv.cut) + 
      labs(x = xlab, y = ylab, title = main, color = "Given conditional survival") + 
      geom_segment(aes(x = 0, y = 1, xend = at[length(at)], yend = 1),lwd = lwd, color = curv.col[at[length(at)] + 1]) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
            panel.border = element_rect(colour = "black", fill=NA, size=0.2),
            axis.ticks = element_line(size=0.2, color="black"),
            axis.ticks.length = unit(0.2, "cm"),
            legend.position = legend.pos,
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 10))
  }
  
  # 输出条件生存概率表
  condsurvmat <- matrix("",
                        nrow = max.condsurv + 1,
                        ncol = floor(max(survdt$futime)) + 1,
                        dimnames = list(paste0("y",0:max.condsurv),
                                        paste0("y",0:floor(max(survdt$futime)))))
  condsurvmat <- as.data.frame(condsurvmat)
  
  for (given_year in 0:max.condsurv) {
    for (reach_year in given_year:floor(max(survdt$futime))) {
      condsurvmat[paste0("y",given_year),paste0("y",reach_year)] <- 
        paste0(conditional_surv_est(basekm, given_year, reach_year)$cs_est * 100,"%")
    }
  }
  write.csv(condsurvmat,file = "output/conditional_survival_matrix.csv",row.names = T,col.names = NA,quote = F)
  
  return(list(condsurvdat = condsurvdat, condsurvcurve = condsurv_plot, condsurvmat = condsurvmat, cond = at, prob_times = seq(0,surv.cut,1), basekm = basekm))
}