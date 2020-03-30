library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

dat = read.table(file=args[1], sep="\t", header=F)

dat$V2 = (1-dat$V2)*100

fig1 = ggplot(dat, aes(x=V1, y=V2)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x='Position within alignment',
       y='% Coverage in alignment') +
  theme(text=element_text(size=20))

fig2 = fig1 +
  scale_x_continuous(limits=c(0,as.numeric(args[3]))) +
  geom_hline(yintercept=95, colour='red')

fig3 = fig1 +
  scale_x_continuous(limits=c(nrow(dat)-as.numeric(args[3]),nrow(dat))) +
  geom_hline(yintercept=95, colour='red')

file1 = paste0(args[2],'gapFrequencies.fullGenome.pdf')
file2 = paste0(args[2],'gapFrequencies.5prime.pdf')
file3 = paste0(args[2],'gapFrequencies.3prime.pdf')

pdf(file1)
fig1
dev.off()

pdf(file2)
fig2
dev.off()

pdf(file3)
fig3
dev.off()


