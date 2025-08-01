# dat is your harmonised file 

library(TwoSampleMR)
library(ggplot2)
#scatter plot
results<-mr(dat)
scatter_plot<-mr_scatter_plot(results,dat)


png('test_scatter.png')
scatter_plot
dev.off()

ggsave(scatter_plot[[1]], file = "scatter_plot.png", width = 7, height = 7)



#forest plot 

res_single <- mr_singlesnp(dat)
forest_plot <- mr_forest_plot(res_single)
forest_plot[[1]]

ggsave(forest_plot[[1]], file = "forest_plot.png", width = 7, height = 7)


#LOO analysis 


res_loo <- mr_leaveoneout(dat)
LOO_plot <- mr_leaveoneout_plot(res_loo)
LOO_plot[[1]]

ggsave(LOO_plot[[1]], file = "LOO_plot.png", width = 7, height = 7)