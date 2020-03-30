library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

dat1 = read.table(file=args[1], sep='\t', header=F)
dat2 = read.table(file=args[2], sep='\t', header=F)

dat1$Position = c(1:nrow(dat1))
dat2$Position = c(1:nrow(dat2))

dat1$V1 = dat1$V1 * 100
dat2$V1 = dat2$V1 * 100

fig1 = ggplot(dat1, aes(x=Position, y=V1)) + 
    geom_bar(stat='identity') +
    theme_bw() +
    labs(x='Position along reference genome',
         y='Frequncy of Mutation (%)')

fig2 = ggplot(dat2, aes(x=Position, y=V1)) +
    geom_bar(stat='identity') +
    theme_bw() +
    labs(x='Position along reference genome',
         y='Frequency of Mutation (%)')

file1 = paste0(args[3], '.mutationFrequency.uncondensed.pdf')
file2 = paste0(args[3], '.mutationFrequency.condensed.pdf')

pdf(file1)
fig1
dev.off()

pdf(file2)
fig2
dev.off()
