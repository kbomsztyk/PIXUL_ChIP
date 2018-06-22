#BAM file 
#MACS2 calling, same parameters. Both used input/control
ctcf<-read.table("AB14_HCT116_CTCF_closest_ENCODE.bed",header=F)
mean(ctcf$V7==0)

ctcf_rev<-read.table("CTCF_ENCODE_closest_to_ours.bed",header=F)
mean(ctcf_rev$V7==0)
sum(ctcf_rev$V7==0)

h3k4me1<-read.table("AB2_HCT116_H3K4me1_closest_ENCODE.bed",header=F)
sum(h3k4me1$V7==0)
mean(h3k4me1$V7==0)

h3k27me3<-read.table("NC5_HTON_H3K27me3_closest_ENCODE.bed",header=F)
sum(h3k27me3$V7==0)
mean(h3k27me3$V7==0)

h3k9me3<-read.table("AB6_HCT116_H3K9me3_closest_ENCODE.bed",header=F)
sum(h3k9me3$V7==0)
mean(h3k9me3$V7==0)

k36me3<-read.table("AB13_H3K36me3_closest_ENCODE.bed",header=F)
sum(k36me3$V7==0)
mean(k36me3$V7==0)

#show that peaks that overlap are high quality peaks in both datasets, so it is fine for other peaks to not overlap#
#K27me3
h3k27me3<-read.table("NC5_HTON_H3K27me3_closest_ENCODE.bed",header=F)
h3k27me3$id<-paste(h3k27me3$V1,h3k27me3$V2,h3k27me3$V3,sep="_")
our_k27<-read.table("NC5_HTON_K27me3_A018_peaks.broadPeak",header=F,sep="\t")
our_k27$id<-paste(our_k27$V1,our_k27$V2,our_k27$V3,sep="_")
peaks_overlap<- our_k27$V8[our_k27$id %in% h3k27me3$id[h3k27me3$V7==0]]
peaks_only_ours<- our_k27$V8[!(our_k27$id %in% h3k27me3$id[h3k27me3$V7==0])]

temp<-stack(list(overlapping=peaks_overlap,unique_to_us=peaks_only_ours))
colnames(temp)[2]<-"Category"
library(ggplot2)
pdf("H3K27me3_our_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<10,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("H3K27me3")
dev.off()

encode_k27<-read.table("H3K27me3_with_control_ENCFF330JMX_peaks.broadPeak",header=F,sep="\t")
encode_k27$id<-paste(encode_k27$V1,encode_k27$V2,encode_k27$V3,sep="_")
h3k27me3<-read.table("NC5_HTON_H3K27me3_closest_ENCODE.bed",header=F)
h3k27me3$id<-paste(h3k27me3$V4,h3k27me3$V5,h3k27me3$V6,sep="_")
peaks_overlap<- encode_k27$V8[encode_k27$id %in% h3k27me3$id[h3k27me3$V7==0]]
peaks_only_encode<- encode_k27$V8[!(encode_k27$id %in% h3k27me3$id[h3k27me3$V7==0])]

temp<-stack(list(overlapping=peaks_overlap,unique_to_ENCODE=peaks_only_encode))
colnames(temp)[2]<-"Category"
pdf("H3K27me3_ENCODE_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<10,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("H3K27me3")
dev.off()

#CTCF
ctcf<-read.table("AB14_HCT116_CTCF_closest_ENCODE.bed",header=F)
ctcf$id<-paste(ctcf$V1,ctcf$V2,ctcf$V3,sep="_")
our_ctcf<-read.table("AB14_HCT116_CTCF_peaks.broadPeak",header=F,sep="\t")
our_ctcf$id<-paste(our_ctcf$V1,our_ctcf$V2,our_ctcf$V3,sep="_")
peaks_overlap<- our_ctcf$V8[our_ctcf$id %in% ctcf$id[ctcf$V7==0]]
peaks_only_ours<- our_ctcf$V8[!(our_ctcf$id %in% ctcf$id[ctcf$V7==0])]

temp<-stack(list(overlapping=peaks_overlap,unique_to_us=peaks_only_ours))
colnames(temp)[2]<-"Category"
pdf("CTCF_our_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<50,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("CTCF")
dev.off()

encode_ctcf<-read.table("CTCF_with_control_ENCFF467JND_peaks.broadPeak",header=F,sep="\t")
encode_ctcf$id<-paste(encode_ctcf$V1,encode_ctcf$V2,encode_ctcf$V3,sep="_")
ctcf<-read.table("AB14_HCT116_CTCF_closest_ENCODE.bed",header=F)
ctcf$id<-paste(ctcf$V4,ctcf$V5,ctcf$V6,sep="_")
peaks_overlap<- encode_ctcf$V8[encode_ctcf$id %in% ctcf$id[ctcf$V7==0]]
peaks_only_encode<- encode_ctcf$V8[!(encode_ctcf$id %in% ctcf$id[ctcf$V7==0])]

temp<-stack(list(overlapping=peaks_overlap,unique_to_ENCODE=peaks_only_encode))
colnames(temp)[2]<-"Category"
pdf("CTCF_ENCODE_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<50,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("CTCF")
dev.off()

#K36me3
k36me3<-read.table("AB13_H3K36me3_closest_ENCODE.bed",header=F)
our_k36me3<-read.table("AB13_K36me3_PIXUL_HCT116_peaks.broadPeak",header=F,sep="\t")
our_k36me3$id<-paste(our_k36me3$V1,our_k36me3$V2,our_k36me3$V3,sep="_")
k36me3$id<-paste(k36me3$V1,k36me3$V2,k36me3$V3,sep="_")
peaks_overlap<- our_k36me3$V8[our_k36me3$id %in% k36me3$id[k36me3$V7==0]]
peaks_only_ours<- our_k36me3$V8[!(our_k36me3$id %in% k36me3$id[k36me3$V7==0])]
temp<-stack(list(overlapping=peaks_overlap,unique_to_us=peaks_only_ours))
colnames(temp)[2]<-"Category"
pdf("H3K36me3_our_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<10,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("H3K36me3")
dev.off()

k36me3<-read.table("AB13_H3K36me3_closest_ENCODE.bed",header=F)
encode_k36me3<-read.table("H3K36me3_with_control_ENCFF341WQZ_peaks.broadPeak",header=F,sep="\t")
encode_k36me3$id<-paste(encode_k36me3$V1,encode_k36me3$V2,encode_k36me3$V3,sep="_")
k36me3$id<-paste(k36me3$V4,k36me3$V5,k36me3$V6,sep="_")
peaks_overlap<- encode_k36me3$V8[encode_k36me3$id %in% k36me3$id[k36me3$V7==0]]
peaks_only_encode<- encode_k36me3$V8[!(encode_k36me3$id %in% k36me3$id[k36me3$V7==0])]
temp<-stack(list(overlapping=peaks_overlap,unique_to_ENCODE=peaks_only_encode))
colnames(temp)[2]<-"Category"
pdf("H3K36me3_ENCODE_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<10,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("H3K36me3")
dev.off()

#k4me1
k4me1<-read.table("AB2_HCT116_H3K4me1_closest_ENCODE.bed",header=F)
our_k4me1<-read.table("AB2_HCT116_H3K4me1_peaks.broadPeak",header=F,sep="\t")
our_k4me1$id<-paste(our_k4me1$V1,our_k4me1$V2,our_k4me1$V3,sep="_")
k4me1$id<-paste(k4me1$V1,k4me1$V2,k4me1$V3,sep="_")
peaks_overlap<- our_k4me1$V8[our_k4me1$id %in% k4me1$id[k4me1$V7==0]]
peaks_only_ours<- our_k4me1$V8[!(our_k4me1$id %in% k4me1$id[k4me1$V7==0])]
temp<-stack(list(overlapping=peaks_overlap,unique_to_us=peaks_only_ours))
colnames(temp)[2]<-"Category"
pdf("H3K4me1_our_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<20,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("H3K4me1")
dev.off()

k4me1<-read.table("AB2_HCT116_H3K4me1_closest_ENCODE.bed",header=F)
encode_k4me1<-read.table("H3K4me1_with_control_ENCFF115MLX_peaks.broadPeak",header=F,sep="\t")
encode_k4me1$id<-paste(encode_k4me1$V1,encode_k4me1$V2,encode_k4me1$V3,sep="_")
k4me1$id<-paste(k4me1$V4,k4me1$V5,k4me1$V6,sep="_")
peaks_overlap<- encode_k4me1$V8[encode_k4me1$id %in% k4me1$id[k4me1$V7==0]]
peaks_only_encode<- encode_k4me1$V8[!(encode_k4me1$id %in% k4me1$id[k4me1$V7==0])]
temp<-stack(list(overlapping=peaks_overlap,unique_to_ENCODE=peaks_only_encode))
colnames(temp)[2]<-"Category"
pdf("H3K4me1_ENCODE_peaks_quality_overlap.pdf",width=4,height=3)
ggplot(temp[temp$values<10,],aes(x=values)) +geom_density(aes(fill=Category),alpha=0.3) + theme(legend.position = c(1,1),legend.justification = c(1,1)) + xlab("-log10(peak p-value)") + ylab("Density") + ggtitle("H3k4me1")
dev.off()

#Motif enrichment
library(dplyr)
ctcf<-read.table("AB14_HCT116_CTCF_closest_ENCODE.bed",header=F)
ctcf$id<-paste(ctcf$V1,ctcf$V2,ctcf$V3,sep="_")
our_ctcf<-read.table("AB14_HCT116_CTCF_peaks.broadPeak",header=F,sep="\t")
our_ctcf$id<-paste(our_ctcf$V1,our_ctcf$V2,our_ctcf$V3,sep="_")
peaks_overlap<- our_ctcf[our_ctcf$id %in% ctcf$id[ctcf$V7==0],]
write.table(peaks_overlap[,1:3],file="CTCF_overlap.bed",sep="\t",quote=F,row.names=F,col.names=F)
our_ctcf<-our_ctcf[our_ctcf$V9>3 & our_ctcf$V7>5,]
write.table(our_ctcf[,1:3],file="our_CTCF_top.bed",sep="\t",quote=F,row.names=F,col.names=F)

h3k27me3<-read.table("NC5_HTON_H3K27me3_closest_ENCODE.bed",header=F)
h3k27me3$id<-paste(h3k27me3$V1,h3k27me3$V2,h3k27me3$V3,sep="_")
our_k27<-read.table("NC5_HTON_K27me3_A018_peaks.broadPeak",header=F,sep="\t")
our_k27$id<-paste(our_k27$V1,our_k27$V2,our_k27$V3,sep="_")
peaks_overlap<- our_k27[our_k27$id %in% h3k27me3$id[h3k27me3$V7==0],]
write.table(peaks_overlap[,1:3],file="H3K27me3_overlap.bed",sep="\t",quote=F,row.names=F,col.names=F)

k36me3<-read.table("AB13_H3K36me3_closest_ENCODE.bed",header=F)
our_k36me3<-read.table("AB13_K36me3_PIXUL_HCT116_peaks.broadPeak",header=F,sep="\t")
our_k36me3$id<-paste(our_k36me3$V1,our_k36me3$V2,our_k36me3$V3,sep="_")
k36me3$id<-paste(k36me3$V1,k36me3$V2,k36me3$V3,sep="_")

peaks_overlap<-our_k36me3[our_k36me3$id %in% k36me3$id[k36me3$V7==0],]
write.table(peaks_overlap[,1:3],file="H3K36me3_overlap.bed",sep="\t",quote=F,row.names=F,col.names=F)


k4me1<-read.table("AB2_HCT116_H3K4me1_closest_ENCODE.bed",header=F)
our_k4me1<-read.table("AB2_HCT116_H3K4me1_peaks.broadPeak",header=F,sep="\t")
our_k4me1$id<-paste(our_k4me1$V1,our_k4me1$V2,our_k4me1$V3,sep="_")
k4me1$id<-paste(k4me1$V1,k4me1$V2,k4me1$V3,sep="_")
peaks_overlap<- our_k4me1[our_k4me1$id %in% k4me1$id[k4me1$V7==0],]
write.table(peaks_overlap[,1:3],file="H3K4me1_overlap.bed",sep="\t",quote=F,row.names=F,col.names=F)

