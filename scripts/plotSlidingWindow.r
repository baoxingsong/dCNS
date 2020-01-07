library(ggplot2)
data = read.table("SORBI_3009G024600_5_prime_slidingwindow_1")


p = ggplot(data, aes(x = V1, y = V3, color=V2)) +  geom_line() + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 12) + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))

png("SORBI_3009G024600_5_prime_slidingwindow_1.png" , width=2000, height=450)
p
dev.off()
pdf("SORBI_3009G024600_5_prime_slidingwindow_1.pdf" , width=33, height=7.5)
p
dev.off()
