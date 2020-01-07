library(ggplot2)
data = read.table("maize_sorghum5p")

mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))), 
  x=c(data$V1,data$V2,data$V4,data$V3), 
  score=c(data$V5,data$V5,data$V5,data$V5),
  y=c(rep("sorghum", nrow(data)),rep("sorghum", nrow(data)),rep("maize", nrow(data)),rep("maize", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}


p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=score), alpha=0.4) +scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36) + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
p
png("maize_sorghum5p.png" , width=2000, height=450)
p
dev.off()
pdf("maize_sorghum5p.pdf" , width=33, height=7.5)
p
dev.off()


data = read.table("maize_sorghum5p_syntenic")

mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))), 
  x=c(data$V1,data$V2,data$V4,data$V3), 
  score=c(data$V5,data$V5,data$V5,data$V5),
  y=c(rep("sorghum", nrow(data)),rep("sorghum", nrow(data)),rep("maize", nrow(data)),rep("maize", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}


p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=(score^2)), alpha=0.4) +scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36) + theme(axis.line = element_blank(),
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(),
                                                 panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                                 panel.background = element_blank(),
                                                 axis.text.y = element_text( colour = "black"),
                                                 axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
p
png("maize_sorghum5p_syntenic.png" , width=2000, height=450)
p
dev.off()
pdf("maize_sorghum5p_syntenic.pdf" , width=33, height=7.5)
p
dev.off()





library(ggplot2)
data = read.table("93983.5queryo")

mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))),
  x=c(data$V1,data$V2,data$V4,data$V3),
  length=c(data$V5,data$V5,data$V5,data$V5),
  pvalue=c(data$V6,data$V6,data$V6,data$V6),
  y=c(rep("maize", nrow(data)),rep("maize", nrow(data)),rep("sorghum", nrow(data)),rep("sorghum", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}

p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=-log10(pvalue)), alpha=0.4) +
    scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36)  + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))


png("SORBI_3009G024600_5_prime_rev_o.png" , width=2000, height=450)
p
dev.off()
pdf("SORBI_3009G024600_5_prime_rev_o.pdf" , width=33, height=7.5)
p
dev.off()




library(ggplot2)
data = read.table("SORBI_3009G024600_5_prime_rev_o_weighted")

mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))),
  x=c(data$V1,data$V2,data$V4,data$V3),
  length=c(data$V5,data$V5,data$V5,data$V5),
  pvalue=c(data$V6,data$V6,data$V6,data$V6),
  y=c(rep("maize", nrow(data)),rep("maize", nrow(data)),rep("sorghum", nrow(data)),rep("sorghum", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}




library(ggplot2)
data = read.table("SORBI_3001G121600_5_prime_o")
data = data[which(data$V6<0.01),]
mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))),
  x=c(data$V1,data$V2,data$V4,data$V3),
  length=c(data$V5,data$V5,data$V5,data$V5),
  pvalue=c(data$V6,data$V6,data$V6,data$V6),
  y=c(rep("sorghum", nrow(data)),rep("sorghum", nrow(data)),rep("maize", nrow(data)),rep("maize", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}

p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=-log10(pvalue)), alpha=0.4) +
    scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(limits = c(18000, 100000), labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36)  + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))

png("SORBI_3001G121600_5_prime_o.png" , width=2000, height=450)
p
dev.off()
pdf("SORBI_3001G121600_5_prime_o.pdf" , width=33, height=7.5)
p
dev.off()








library(ggplot2)
data = read.table("SORBI_3009G024600_5_prime_many_to_many_seed_xtend")
data = data[which(data$V6<0.01),]
mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))),
  x=c(data$V1,data$V2,data$V4,data$V3),
  length=c(data$V5,data$V5,data$V5,data$V5),
  pvalue=c(data$V6,data$V6,data$V6,data$V6),
  y=c(rep("maize", nrow(data)),rep("maize", nrow(data)),rep("sorghum", nrow(data)),rep("sorghum", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}

p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=-log10(pvalue)), alpha=0.4) +
    scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(limits = c(0, 100000), labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36)  + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))

png("1o.png" , width=2000, height=450)
p
dev.off()
pdf("1o.pdf" , width=33, height=7.5)
p
dev.off()






library(ggplot2)
data = read.table("SORBI_3009G024600_5_prime_many_to_many_2gaps_o")

mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))),
  x=c(data$V1,data$V2,data$V4,data$V3),
  length=c(data$V5,data$V5,data$V5,data$V5),
  pvalue=c(data$V6,data$V6,data$V6,data$V6),
  y=c(rep("maize", nrow(data)),rep("maize", nrow(data)),rep("sorghum", nrow(data)),rep("sorghum", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}

p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=-log10(pvalue)), alpha=0.4) +
    scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(limits = c(0, 100000), labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36)  + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))


png("2.png" , width=2000, height=450)
p
dev.off()
pdf("2.pdf" , width=33, height=7.5)
p
dev.off()







library(ggplot2)
data = read.table("3sorghumo")
mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))),
  x=c(data$V1,data$V2,data$V4,data$V3),
  length=c(data$V5,data$V5,data$V5,data$V5),
  pvalue=c(data$V6,data$V6,data$V6,data$V6),
  y=c(rep("maize", nrow(data)),rep("maize", nrow(data)),rep("sorghum", nrow(data)),rep("sorghum", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}

p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=-log10(pvalue)), alpha=0.4) +
    scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(limits = c(0, 100000), labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36)  + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))

png("3o.png" , width=2000, height=450)
p
dev.off()
pdf("3o.pdf" , width=33, height=7.5)
p
dev.off()







library(ggplot2)
data = read.table("4sorghumo")
mypositions1 <- data.frame(
  id = c(c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data)), c(1:nrow(data))),
  x=c(data$V1,data$V2,data$V4,data$V3),
  length=c(data$V5,data$V5,data$V5,data$V5),
  pvalue=c(data$V6,data$V6,data$V6,data$V6),
  y=c(rep("maize", nrow(data)),rep("maize", nrow(data)),rep("sorghum", nrow(data)),rep("sorghum", nrow(data)))
)

changetoM <- function ( position ){
  100000 - position
}

p = ggplot(mypositions1, aes(x = x, y = y)) +  geom_polygon(aes( group = id, fill=-log10(pvalue)), alpha=0.4) +
    scale_fill_gradient(low="blue", high="red")+ scale_x_continuous(limits = c(0, 100000), labels=changetoM) + labs(x="", y="", title="")+
  theme_bw() +theme_grey(base_size = 36)  + theme(axis.line = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                           panel.background = element_blank(),
                                           axis.text.y = element_text( colour = "black"),
                                           axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))

png("4o.png" , width=2000, height=450)
p
dev.off()
pdf("4o.pdf" , width=33, height=7.5)
p
dev.off()


