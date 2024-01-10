chl<-read.table("Desktop/CMIP6/chl.txt", header = FALSE)
chl2100<-chl[grep("2100",chl$V1),]
write.table(chl2100,"Desktop/CMIP6/chl_2100.txt", row.names = FALSE, quote = FALSE, )

nh4<-read.table("Desktop/CMIP6/nh4.txt", header = FALSE)
nh42100<-nh4[grep("2100",nh4$V1),]
write.table(nh42100,"Desktop/CMIP6/nh4_2100.txt", row.names = FALSE, quote = FALSE, )

o2<-read.table("Desktop/CMIP6/o2.txt", header = FALSE)
o22100<-o2[grep("2100",o2$V1),]
write.table(o22100,"Desktop/CMIP6/o2_2100.txt", row.names = FALSE, quote = FALSE, )

ph<-read.table("Desktop/CMIP6/ph.txt", header = FALSE)
ph2100<-ph[grep("2100",ph$V1),]
write.table(ph2100,"Desktop/CMIP6/ph_2100.txt", row.names = FALSE, quote = FALSE, )

fe<-read.table("Desktop/CMIP6/fe.txt", header = FALSE)
fe2100<-fe[grep("2100",fe$V1),]
write.table(fe2100,"Desktop/CMIP6/fe_2100.txt", row.names = FALSE, quote = FALSE, )

sal<-read.table("Desktop/CMIP6/sal.txt", header = FALSE)
sal2100<-sal[grep("2100",sal$V1),]
write.table(sal2100,"Desktop/CMIP6/sal_2100.txt", row.names = FALSE, quote = FALSE, )

temp<-read.table("Desktop/CMIP6/temp.txt", header = FALSE)
temp2100<-temp[grep("2100",temp$V1),]
write.table(temp2100,"Desktop/CMIP6/temp_2100.txt", row.names = FALSE, quote = FALSE, )

