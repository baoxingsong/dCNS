library(ggplot2)
data = read.table("maize_sorghum5p_MSA_plot")

changetoM <- function ( position ){
  100000 - position
}

p = ggplot(data, aes(x = V1, y = V2)) +  geom_polygon(aes( group = V3, fill=V4), alpha=0.4) +
    scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(labels=changetoM) + labs(x="", y="", title="")+labs(fill = "length") +
  theme_bw() +theme_grey(base_size = 36)  + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
p

png("maize_sorghum5p_MSA_plot.png" , width=1332, height=450)
p
dev.off()
pdf("maize_sorghum5p_MSA_plot.pdf" , width=22, height=7.5)
p
dev.off()
