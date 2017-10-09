#!/software/R-3.3.0/bin/Rscript --vanilla
args=commandArgs(trailingOnly = T)
freal=args[1]
chr=args[2]
out=args[3]
include_dups=as.logical(as.integer(args[4]))

cat(paste("Opening", freal, "...\n"))
mar.default <- par('mar')

library(data.table)
pdf(paste0(out, ".pdf"), width=15	)
        allbreaks=fread(freal)
cat("Plotting...\n")

alldepth=allbreaks
        confidence=0.05
        
mind=min(allbreaks$depth)
maxd=max(allbreaks$depth)
breaks=seq(0, ceiling(max(allbreaks$depth)), by=0.02)
histo=cut(allbreaks$depth, breaks)
        par(new=F, mar=mar.default)
    mar.default <- par('mar')
    mar.left <- mar.default
    mar.right <- mar.default
    
    mar.left[4] <- 0
    mar.right[2] <- 0
    mar.one=c(0,2,1.5,0)
mar.two=c(0,0,1.5,0)
mar.three=c(2, 2, 0, 0)
par(fig=c(0,0.8, 0.3,1), mar=mar.one)
plot(c(allbreaks$min, allbreaks$max), c(allbreaks$depth, allbreaks$depth), type="n", ylim=c(0, ceiling(max(allbreaks$depth))), xaxt="n", main=freal)
segments(x0=allbreaks$min, y0=allbreaks$depth, x1=allbreaks$max, lwd=3, col=adjustcolor("darkslategray", alpha=0.2))


## Calling
library(mixtools)
mvect=seq(0, ceiling(max(allbreaks$depth)), by=0.5)
mix=normalmixEM(allbreaks$depth, k=length(mvect), mu=mvect, sigma=1, mean.constr=mvect, sd.constr=rep('a', length(mvect)))
parm=cbind(mix$mu, mix$sigma)
apply(parm, 1, function(x){
    abline(h=qnorm(1-confidence, mean=x[1], sd=x[2], lower.tail=T), lty=2, col="gray")
    abline(h=qnorm(1-confidence, mean=x[1], sd=x[2], lower.tail=F), lty=2, col="gray")
    
})
discr=seq(min(alldepth$min), max(alldepth$max), by=5000)
callz=t(sapply(allbreaks$depth, function(meas) {
    botparm=parm[parm[,1]==floor(meas*2)/2,]
    topparm=parm[parm[,1]==ceiling(meas*2)/2,]
    topp=max(pnorm(meas, mean=botparm[1], sd=botparm[2], lower.tail=F), pnorm(meas, topparm[1], sd=topparm[2], lower.tail=T))
    topn=ifelse(pnorm(meas, mean=botparm[1], sd=botparm[2], lower.tail=F)>pnorm(meas, topparm[1], sd=topparm[2], lower.tail=T), botparm[1], topparm[1])
    return(c(meas, 2*topp, topn))
}))

allbreaks$assigned_depth=callz[,3]
allbreaks$p_value=callz[,2]
tocall=allbreaks[allbreaks$assigned_depth<1, ] # replace this with != to call dups and dels
discr=seq(min(alldepth$min), max(alldepth$max), by=5000)

numhethom=sapply(discr, function(x){
    in_region=tocall[tocall$min<=x & tocall$max>=x,]
    nrow(in_region[in_region$p_value>confidence,])
    })
numhethom2=sapply(discr, function(x){
    in_region=tocall[tocall$min<=x & tocall$max>=x,]
    if(nrow(in_region)==0){return(0)}
    adjust_factor=nrow(in_region[in_region$p_value<confidence,])
    if(!is.finite(adjust_factor)){adjust_factor=1}
    nrow(in_region[in_region$p_value>confidence,])/adjust_factor
    })

a=rle(numhethom>0) ## arbitrary threshold here we come
svs=t(rbind(a$values, a$lengths, min(discr)+cumsum(a$lengths)*5000))
svs=svs[svs[,1]==1,,drop=F]
if(!all(is.na(svs))){
rect(xleft=svs[,3]-svs[,2]*5000, xright=svs[,3], ybottom=-1e6, ytop=1e6, border=NA, col=adjustcolor("blue3", alpha=0.4))
}

a=rle(numhethom2>1) ## arbitrary threshold here we come
svs=t(rbind(a$values, a$lengths, min(discr)+cumsum(a$lengths)*5000))
svs=svs[svs[,1]==1,,drop=F]
if(!all(is.na(svs))){
rect(xleft=svs[,3]-svs[,2]*5000, xright=svs[,3], ybottom=-1e6, ytop=1e6, border=NA, col=adjustcolor("deeppink", alpha=0.4))
}

## Plot right histo
par(fig=c(0.8,1, 0.3,1), new=T, mar=mar.two)
plot(0, type="n", ylim=c(0,ceiling(max(allbreaks$depth))), xlim=c(0, max(table(histo))), xaxt="n", yaxt="n", xlab="", ylab="")
grid()
segments(x0=rep(0, length(breaks)), x1=table(histo), y0=breaks, y1=breaks, lwd=4)
par(fig=c(0,0.8, 0,0.3), new=T, mar=mar.three)
plot(discr, numhethom,type="l", lty=2, col="gray")
points(discr, numhethom, pch=".", cex=2, col="blue3")
points(discr, numhethom2, pch=".", cex=2, col="deeppink")
points(discr, numhethom2, type="l", lty=2, col="deeppink")


dev.off()
cat("Generating diagnostic dataset...\n")

calls=svs
calls[,2]=svs[,3]-svs[,2]*5000

cnv_stats=data.frame(t(apply(calls, 1, function(x){
    relevant=allbreaks[allbreaks$min<x[3] & allbreaks$max>x[2],]
    avgdelp=mean(relevant[relevant$assigned_depth<1,]$p_value)
    avgdupp=mean(relevant[relevant$assigned_depth>1,]$p_value)
    del_in_conf=nrow(relevant[relevant$assigned_depth<1 & relevant$p_value>confidence,])
    del_out_conf=nrow(relevant[relevant$assigned_depth<1 & relevant$p_value<confidence,])
    dup_in_conf=nrow(relevant[relevant$assigned_dept>1 & relevant$p_value>confidence,])
    dup_out_conf=nrow(relevant[relevant$assigned_depth>1 & relevant$p_value<confidence,])
    deldup=nrow(relevant[relevant$assigned_depth<1,])/nrow(relevant[relevant$assigned_depth>1,])
    delaf=(nrow(relevant[relevant$assigned_depth==0,])*2+nrow(relevant[relevant$assigned_depth==0.5,]))/nrow(relevant)
    ret=c(x, 
          table(cut(relevant$assigned_depth, c(-1, 0.25, 0.75, 1.25, 1e6), labels=c("hom", "het", "normal", "dup"))),
          delaf,
          deldup,
          avgdelp,
          del_in_conf,
          del_out_conf,
          del_in_conf/(del_out_conf+del_in_conf),
          avgdupp,
          dup_in_conf,
          dup_out_conf,
          dup_in_conf/(dup_out_conf+dup_in_conf)
         )
    names(ret)[1:3]=c("chr", "start", "stop")
    names(ret)[8:9]=c("del_af","del_ratio")
    names(ret)[10:17]=c("avg_del_p", "del_hc", "del_lc", "del_hc_total_ratio", "avg_dup_p", "dup_hc", "dup_lc", "dup_hc_total_ratio")
    ret
})))
cnv_stats$chr=chr
write.table(cnv_stats, paste0(out, ".stats"), row.names=F, quote=F)

cat("Generating PLINK dataset (")
cat(ifelse(include_dups, "including dups, writing .cnv/.fam)\n", "excluding dups, writing .ped/.map and coercing dups to 1)\n"))

if(include_dups){
cnv_df=apply(calls, 1, function(x){
    relevant=allbreaks[allbreaks$min<x[3] & allbreaks$max>x[2],]
    cbind(relevant$sample, relevant$sample, rep(chr, nrow(relevant)), rep(x[2], nrow(relevant)), rep(x[3], nrow(relevant)), relevant$assigned_depth*2, relevant$p_value, rep(0, nrow(relevant)))
    
 })

cnv_df=do.call(rbind, cnv_df)
 write.table(cnv_df, paste0(out, ".cnv"), row.names=F, col.names=F, quote=F)

}else{
# exclude dups, plink dataset
alldepth=allbreaks
samples=unique(alldepth$sample)
geno=NULL
forgetme=apply(as.matrix(calls), 1, function(x){
    copy1=rep(1, length(samples))
    copy2=rep(1, length(samples))
    names(copy1)=samples
    names(copy2)=samples
    a=alldepth[alldepth$min<=x[3] & alldepth$max>=x[2],];
    homs=a$sample[a$assigned_depth==0]			# This is really quite permissive. Whatever they're closer to, 0 or 0.5, they get assigned to.
    hets=a$sample[a$assigned_depth==0.5]
    copy1[c(homs,hets)]=2
    copy2[homs]=2
    geno<<-cbind(geno, copy1,copy2)
})
geno=cbind(samples, samples, rep(0, length(samples)), rep(0, length(samples)), rep(0, length(samples)), 
           rep(-9, length(samples)), geno)
calls[,1]=chr
calls=cbind(calls,calls[,3])
calls[,3]=0

write.table(calls, paste0(out, ".map"), quote=F, row.names=F, col.names=F)
write.table(geno, paste0(out, ".ped"), quote=F, row.names=F, col.names=F)
}

cat("End of script.\n")
