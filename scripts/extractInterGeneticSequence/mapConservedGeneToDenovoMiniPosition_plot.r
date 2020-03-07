library("ggplot2")

data <- read.table("miniDistance")

changetoK <- function ( position ){
  position=(10^position)/1000;
  paste(position, "K", sep="")
}

changetoK2 <- function ( position ){
  position=position/1000;
  paste(position, "K", sep="")
}

pdf(file="conserved_hist.pdf")
ggplot(data=data, aes(x=log10(data$V2))) +
  geom_histogram(color="black", fill="white", binwidth = 0.25)+
  scale_x_continuous(labels=changetoK)+
  labs(y="Frequency", x="minimum extent of flanking",  title="")+
  theme_bw() +theme_grey(base_size = 31)+
  theme(axis.line = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()

pdf(file="conserved_accumulation.pdf")
ggplot(data=data, aes(x=log10(data$V2))) +
  geom_step(aes(y=..y..),stat="ecdf")+
  scale_x_continuous(labels=changetoK)+
  labs(y="Accumulation proportion\nof mapped genes", x="minimum extent of flanking",  title="")+
  theme_bw() +theme_grey(base_size = 31)+
  theme(axis.line = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()


pdf(file="conserved_small_accumulation.pdf")

ggplot(data=data, aes(x=(data$V2))) +
  geom_step(aes(y=..y..),stat="ecdf")+
  scale_x_continuous(labels=changetoK2)+
  coord_cartesian(xlim=c(5000, 50000))+
  geom_hline(yintercept = 0.5) + 
  labs(y="Accumulation proportion\nof mapped genes", x="minimum extent of flanking",  title="")+
  theme_bw() +theme_grey(base_size = 31)+
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()

