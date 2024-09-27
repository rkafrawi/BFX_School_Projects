out_file = 'inner_dist_RK_real'
pdf('inner_dist_RK_real.inner_distance_plot.pdf')
fragsize=rep(c(-247.5,-242.5,-237.5,-232.5,-227.5,-222.5,-217.5,-212.5,-207.5,-202.5,-197.5,-192.5,-187.5,-182.5,-177.5,-172.5,-167.5,-162.5,-157.5,-152.5,-147.5,-142.5,-137.5,-132.5,-127.5,-122.5,-117.5,-112.5,-107.5,-102.5,-97.5,-92.5,-87.5,-82.5,-77.5,-72.5,-67.5,-62.5,-57.5,-52.5,-47.5,-42.5,-37.5,-32.5,-27.5,-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5,82.5,87.5,92.5,97.5,102.5,107.5,112.5,117.5,122.5,127.5,132.5,137.5,142.5,147.5,152.5,157.5,162.5,167.5,172.5,177.5,182.5,187.5,192.5,197.5,202.5,207.5,212.5,217.5,222.5,227.5,232.5,237.5,242.5,247.5),times=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,49,631,549,591,635,692,778,1071,1614,2467,3810,5820,8914,13177,18980,24179,30826,36146,40412,44101,48215,48389,49792,49007,48154,45683,43616,40804,37813,34090,30615,28166,24895,22351,20048,17184,15222,13515,11981,10325,9036,7915,6941,5967,5365,4584,3964,3560,2956,2727,2362,2092,1954,1691,1592,1445,1344,1189,1081))
frag_sd = sd(fragsize)
frag_mean = mean(fragsize)
frag_median = median(fragsize)
write(x=c("Name","Mean","Median","sd"), sep="	", file=stdout(),ncolumns=4)
write(c(out_file,frag_mean,frag_median,frag_sd),sep="	", file=stdout(),ncolumns=4)
hist(fragsize,probability=T,breaks=100,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")
lines(density(fragsize,bw=10),col='red')
dev.off()
