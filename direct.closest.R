#!/software/R-3.3.0/bin/Rscript --vanilla
args=commandArgs(trailingOnly = T)
freal=args[1]
chrdp=args[2]
out=args[3]
chr=args[4]
start=args[5]
end=args[6]


genotype_closest=function(alldepth, chrdp, stat){
  ### This function performs genotyping according to mean depth in the region, ignoring depth segments.
    ## say hello
    cat(" [from script (genotype_closest):] Genotyping using interval-wide means\n")
    ## prepare
	print(stat)
        chr=unique(alldepth$chr)
        alldepth$chr=NULL
        pos=alldepth$pos
        alldepth$pos=NULL
        nms=colnames(alldepth)
        #alldepth=alldepth[,(nms[length(nms)]):=NULL]
    ## then do the plot
    allout=NULL
    for(i in 1:nrow(stat)){
        chr=stat[i,]$chr
        start=stat[i,]$start
        stop=stat[i,]$stop
        
        dmin=alldepth[pos>start & pos<stop,];
        numsnp=nrow(dmin)
        obase=paste(chr, start, stop, numsnp, sep=".")
	print(dim(dmin))
        u=apply(dmin, 2, function(x) mean(as.numeric(as.character(x)), na.rm=T));
        topcp=ceiling(max(u/chrdp$V2, na.rm=T)*2)/2;
        print(topcp)
	print(max(u/chrdp$V2, na.rm=T)*2)
	genos=seq(0, topcp, by=0.5);
        attgeno=sapply(u/chrdp$V2, function(x){
      testv=abs(x-genos);
      ret=(0:(length(genos)-1))[testv==min(testv)];
      if(length(ret)==1){return(ret)}
      else if (1 %in% ret){return(1)}
      else if(ret[1]>1){return(min(ret, na.rm=T))}
      else{return(max(ret, na.rm=T))}});
        
        # outfile
        dgeno=sapply(u/chrdp$V2, function(x){testv=abs(x-genos);return(min(testv))})
        outdf=cbind(chrdp, rep(chr, nrow(chrdp)), rep(start, nrow(chrdp)),rep(stop, nrow(chrdp)), u, attgeno, dgeno)
        colnames(outdf)=c("sample","chrwide_depth","chr","start", "end" ,  "region_depth", "copy_number", "inverse_confidence")
        #write.table(outdf, paste(obase, "out", sep="."), quote=F, row.names=F, sep="\t")
        allout=rbind(allout, outdf)
    }
    return(allout)

}


plot_cnv_qc=function(outdf, alldepth){
                alldepth$chr=NULL
        pos=alldepth$pos
        alldepth$pos=NULL
          chr=unique(outdf$chr)
          start=unique(outdf$start)
          stop=unique(outdf$end)
  # general plot cfg
  #png(paste(obase, "png", sep="."), width=1500, height=800)
  par(mfrow=c(1,2), mar = c(2,3,1,2))

  # plot 

  # plot mean-based
  # needed: attgeno, u/chrdp, dmin, topcp, start, stop
  cat(" [from script (plot_cnv_qc):] Generating plots for interval ", chr, ":", start,"-", stop,"\n")

  attgeno=outdf$copy_number
  u=outdf$region_depth
  chrdp=outdf$chrwide_depth
  topcp=ceiling(max(u/chrdp, na.rm=T)*2)/2;
  dmin=alldepth[pos>start & pos<stop,];
        
  plot(u/chrdp, col=attgeno+1, pch=19, xaxt="n", xlab="samples", ylab="normalised depth", 
       main=paste(chr, paste(start, stop, sep="-"), sep=":"), ylim=c(0, topcp+0.1));
  lal=apply(dmin, 1, function(x){as.numeric(as.character(x))/chrdp});
  nsnps=ncol(lal)
  k=ceiling(nsnps/10)
  plot(0, xlim=c(1, nsnps), ylim=c(0, topcp+0.1), type="n", xaxt="n", xlab="position", ylab="normalised depth", main=paste(ncol(lal), "SNPs"));
    cat(" [from script (plot_cnv_qc) (DEBUG):]      Drawing lines with n=", nsnps, " and k=", k, " for event ",chr, ":",start, "-", stop, " ", nrow(dmin)/9, "\n")
  # lal=lal[1:3,]
  apply(lal, 1, function(x) {lines(rollmean(x, k))})
  #dev.off()
}


library(data.table)
library(zoo)
d=fread(freal, header=T)
alld=fread(chrdp, header=F)
stat=data.table(chr=chr, start=start, stop=end)
closest=genotype_closest(d, alld, stat)
print(head(closest))

towrite=closest
towrite$norm_depth=towrite$region_depth/towrite$chrwide_depth
write.table(towrite, paste(out, "tsv", sep="."), quote=F, col.names=T, row.names=F, sep="\t")

png(paste(out, "png", sep="."), width=1000)
plot_cnv_qc(closest, d)
dev.off()

dm=data.table(FID=closest$sample, IID=closest$sample, FI=0, MI=0, SEX=0, PH=0)
dm$copy1=2
dm$copy2=2
dm$copy1[closest$copy_number<2]=1
dm$copy2[closest$copy_number==0]=1
write.table(dm, paste(out, "ped", sep="."), quote=F, col.names=F, row.names=F)
write.table(data.frame(chr=chr, id=paste0(chr, ":", start, "-", end), lol=0, start=start), paste(out, "map", sep="."), quote=F, col.names=F, row.names=F)


