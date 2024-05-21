
mr_analysis <- function(exp_data,out_data){
#----harmonise data---
Mydf<- harmonise_data(exposure_dat=data.frame(exp_data),
                      outcome_dat=data.frame(out_data),
                      action= 2)

Mydf$exposure <- exposure
Mydf$outcome <- outcome
write.xlsx(Mydf,"analysis results/mr_data.xlsx")

#接下来就可以进行MR分析了，在这里作者定义了5种方法，包括固定效应和随机效应模型
res<-mr(Mydf)
res
write.xlsx(res,"analysis results/mr_analysis_result.xlsx")

#生成了结果，结论是一样的，精神病对前臂骨质疏松没关联。

#进行了一个（MR-PRESSO）检验，这个也是多水平效应检验，P值应该要大于0.05
pres <- mr_presso(BetaOutcome="beta.outcome",
          BetaExposure ="beta.exposure",
          SdOutcome ="se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest =TRUE,DISTORTIONtest = TRUE,
          data =Mydf, NbDistribution = 1000,SignifThreshold = 0.05)

pres
write.xlsx(pres,"analysis results/mr_presso.xlsx")

# 查看离群值
Outlier <- pres$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`

if(!is.null(Outlier))
  Mydf <- Mydf[-Outlier,] # 剔除离群值61和70

#异质性检验，没有异质性的。
mr_hete <- mr_heterogeneity(Mydf,
                 method_list=c("mr_egger_regression", "mr_ivw"))

write.xlsx(mr_hete,"analysis results/mr_heterogeneity.xlsx")

#多水平校验，这里是没有多水平效应的
pleio<- mr_pleiotropy_test(Mydf)
pleio

write.xlsx(pleio,"analysis results/mr_pleiotropy_test.xlsx")

#生成OR和可信区间
OR<-generate_odds_ratios(res)
write.xlsx(OR,"analysis results/OR value.xlsx")

#Leave-one-out analysis是指逐步剔除SNP后观察剩余的稳定性，理想的是剔除后变化不大
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6)
par(mar = c(0,4,2,0))
png(filename = "analysis results/1.mr_leaveoneout.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

single<- mr_leaveoneout(Mydf)
p1<- mr_leaveoneout_plot(single)
print(p1)

dev.off()

#散点图
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/2.mr_scatter_plot.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

p2 <- mr_scatter_plot(res,Mydf)
print(p2)

dev.off()


#绘制森林图
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/3.mr_forest_plot.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

res_single<- mr_singlesnp(Mydf)
p3 <- mr_forest_plot(res_single)
print(p3)

dev.off()


#---forestploter---

library(grid)
#install.packages("forestploter")
library(forestploter)

fore_df <- data.frame(Exposure=exposure,Outcome=outcome,method="",Sample_size=265431,
                      OR=NA,p95CI=NA,pvalue="",low=NA,high=NA,se=NA,OR_95DI="")
or_sele <- dplyr::select(.data = OR,exposure,outcome,method,nsnp,or,
                         exposure,pval,or_lci95,or_uci95,se,lo_ci,up_ci)

names(or_sele) <- names(fore_df)
or_sele$Sample_size <- ""
or_sele$OR <- round(OR$or,3) %>% as.numeric()
or_sele$p95CI <- paste0(round(OR$or_lci95,3),"-",round(OR$or_uci95,3))
or_sele$pvalue <- round(OR$pval,4) %>% as.numeric()
or_sele$low <- round(OR$or_lci95,3) %>% as.numeric()
or_sele$high <- round(OR$or_uci95,3) %>% as.numeric()
or_sele$se <- round(OR$se,3) %>% as.numeric()
or_sele$OR_95DI <- paste0(round(or_sele$OR,3)," (",or_sele$p95CI,")")

df_forest <- rbind(fore_df,or_sele)

min(df_forest$low[-1])
max(df_forest$high[-1])

spac=round(max(df_forest$high[-1]),0)+3
df_forest$" " <- c("               ")# 决定森林线图的宽度


colnames(df_forest)[c(6,11)] <- c("95%CI","OR (95%CI)")

tm <- forest_theme(base_size = 11,  #文本的大小
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,   #可信区间点的形状
                   ci_col = "#762a83",    #CI的颜色
                   ci_fill = "blue",     #ci颜色填充
                   ci_alpha = 0.8,        #ci透明度
                   ci_lty = 1,            #CI的线型
                   ci_lwd = 2,          #CI的线宽
                   ci_Theight = 0.2, # Set an T end at the end of CI  ci的高度，默认是NULL
                   # Reference line width/type/color   参考线默认的参数，中间的竖的虚线
                   refline_lwd = 1,       #中间的竖的虚线
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color  垂直线宽/类型/颜色   可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lwd = 2,              #可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders   更改填充和边框的摘要颜色
                   summary_fill = "yellow",       #汇总部分大菱形的颜色
                   summary_col = "#4575b4",
                   # Footnote font size/face/color  脚注字体大小/字体/颜色
                   footnote_cex = 1,
                   footnote_fontface = "italic",
                   footnote_col = "red")

while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/4.forest_plot.png",
    width=1000, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)


a1=round(min(df_forest$low[-1]),0)
a3=round(max(df_forest$high[-1]),0)
a2 <- round((a1+a3)/2,0)
a4=a3+round((a3-a2)/2,0)
if(a4==a3)
  a4=a4+1

p4 <- forestploter::forest(data=df_forest[,c(1:4,11,12,7)],
                           est = as.numeric(df_forest$OR), #效应值
                           lower = df_forest$low,   #可信区间下限
                           upper = df_forest$high,   #可信区间上限
                           sizes = 1.5*df_forest$se,
                           ci_column = 6,   #在那一列画森林图，要选空的那一列
                           ref_line = 1,
                           arrow_lab = c("protective factor", "risk factor"),
                           xlim = c(0, 6),
                           ticks_at = c(a1,a2,a3,a4),
                           footnote = "",
                           theme = tm)
p4

print(p4)

dev.off()#关闭Plots

#绘制漏斗图，主要是看蓝线周围的散点是否对称
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/5.funnel_plot.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

p5 <- mr_funnel_plot(res_single)
print(p5)

dev.off()#关闭Plots

mr_pict <-list(pict=p5)
return(mr_pict)

}



