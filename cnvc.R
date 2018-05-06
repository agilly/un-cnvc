#!/software/R-3.3.0/bin/Rscript --vanilla
args=commandArgs(trailingOnly = T)
freal=args[1]
chrdp=args[2]
out=args[3]
include_dups=as.logical(as.integer(args[4]))
complexity=args[5]
method=args[6]
nogeno=as.logical(as.integer(args[7]))

cat(" [from script (main):] Sourcing functions.\n")


write_plink=function(m, outf, which_method, include_dups){
  cat(" [from script (write_plink):] Writing PLINK ", ifelse(include_dups, ".cnv (including dups)", ".ped/map (excluding dups)"), " using method ", method, ".\n")

  if(which_method %in% c("segment", "both")){
    if(include_dups){
        ## Write .cnv file
        # FID IID CHR BP1 BP2 COPYNUMBER SCORE DUMMY_0
        ped=data.frame(FID=m$sample, IID=m$sample, CHR=m$chr, BP1=m$start, BP2=m$stop, TYPE=2*m$segment_total_assigned_depth, SCORE=m$segment_weighted_pvalue, SITES=0)
        write.table(ped, paste(outf, "segment", "cnv", sep="."), quote=F, col.names=F, row.names=F, sep="\t")
        ## + FAM file
        ## FID IID IF IM SEX PHENO
        splv=unique(m$sample)
        fam=data.frame(FID=splv, IID=splv, IF=0, IM=0, SEX=0, PHENO=0)
        write.table(fam, paste(outf, "segment", "fam", sep="."), quote=F, col.names=F, row.names=F, sep="\t")
      }else{
        ## Write .ped/.map
        out=data.frame(FID=unique(m$sample), IID=unique(m$sample), IF=0, IM=0, SEX=0, PHENO=0)
        out$IID=as.character(out$IID)
        map=NULL
        for(sv_start in unique(m$start)){
          ## for each SV, we add the double genotype vector
          tmpdf=m[m$start==sv_start,]
          rdf=data.frame(IID=as.character(tmpdf$sample), geno=tmpdf$segment_total_assigned_depth)
          rdf$copy1=2
          rdf$copy2=2
          rdf$copy1[rdf$geno<1]=1
          rdf$copy2[rdf$geno==0]=1
          rdf$geno=NULL
          out=merge(out, rdf, by="IID")


          ## we create the MAP field entry
          chr=unique(tmpdf$chr)[1]
          stop=unique(tmpdf$stop)[1]
          map=rbind(map, c(chr, paste0(chr, ":", sv_start, "-", stop),0,sv_start))
        }
        write.table(out, paste(outf, "segment", "ped", sep="."), quote=F, col.names=F, row.names=F, sep="\t")
        write.table(map, paste(outf, "segment", "map", sep="."), quote=F, col.names=F, row.names=F, sep="\t")

      }
  }
  if(which_method %in% c("means", "both")){
    if(include_dups){
        ## Write .cnv file
        # FID IID CHR BP1 BP2 COPYNUMBER SCORE DUMMY_0
        ped=data.frame(FID=m$sample, IID=m$sample, CHR=m$chr, BP1=m$start, BP2=m$stop, TYPE=2*m$segment_total_assigned_depth, SCORE=m$segment_weighted_pvalue, SITES=0)
        write.table(ped, paste(outf, "means", "cnv", sep="."), quote=F, col.names=F, row.names=F, sep="\t")
        ## + FAM file
        ## FID IID IF IM SEX PHENO
        splv=unique(m$sample)
        fam=data.frame(FID=splv, IID=splv, IF=0, IM=0, SEX=0, PHENO=0)
        write.table(fam, paste(outf, "means", "fam", sep="."), quote=F, col.names=F, row.names=F, sep="\t")
      }else{
        ## Write .ped/.map
        out=data.frame(FID=unique(m$sample), IID=unique(m$sample), IF=0, IM=0, SEX=0, PHENO=0)
        out$IID=as.character(out$IID)
        map=NULL
        for(sv_start in unique(m$start)){
          ## for each SV, we add the double genotype vector
          tmpdf=m[m$start==sv_start,]
          rdf=data.frame(IID=as.character(tmpdf$sample), geno=tmpdf$means_assigned_depth)
          rdf$copy1=2
          rdf$copy2=2
          rdf$copy1[rdf$geno<1]=1
          rdf$copy2[rdf$geno==0]=1
          rdf$geno=NULL
          out=merge(out, rdf, by="IID")

          ## we create the MAP field entry
          chr=unique(tmpdf$chr)[1]
          stop=unique(tmpdf$stop)[1]
          map=rbind(map, c(chr, paste0(chr, ":", sv_start, "-", stop),0,sv_start))
        }

      }
      write.table(out, paste(outf, "means", "ped", sep="."), quote=F, col.names=F, row.names=F, sep="\t")
      write.table(map, paste(outf, "means", "map", sep="."), quote=F, col.names=F, row.names=F, sep="\t")

  }
}

plot_window = function(allbreaks, confidence, sv1, sv2, discr, numhethom, numhethom2){
  #pdf(paste0(out, ".pdf"), width=15 )
  cat(" [from script (plot_window):] Plotting region...\n")
  # plotting parameters
  mar.default <- par('mar')
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

  # generating histo
  mind=min(allbreaks$depth)
  maxd=max(allbreaks$depth)
  breaks=seq(0, ceiling(max(allbreaks$depth)), by=0.02)
  histo=cut(allbreaks$depth, breaks)

  # Plotting depth segments
  plot(c(allbreaks$min, allbreaks$max), c(allbreaks$depth, allbreaks$depth), type="n", ylim=c(0, ceiling(max(allbreaks$depth))), xaxt="n", main=freal)
  segments(x0=allbreaks$min, y0=allbreaks$depth, x1=allbreaks$max, lwd=2, col=adjustcolor("darkslategray", alpha=0.5))

  # draws confidence lines on plot
  apply(parm, 1, function(x){
      abline(h=qnorm(1-confidence, mean=x[1], sd=x[2], lower.tail=T), lty=2, col="gray")
      abline(h=qnorm(1-confidence, mean=x[1], sd=x[2], lower.tail=F), lty=2, col="gray")
  })

  # first calling method
  if(!all(is.na(sv1))){
  rect(xleft=sv1[,3]-sv1[,2]*5000, xright=sv1[,3], ybottom=-1e6, ytop=1e6, border=NA, col=adjustcolor("blue3", alpha=0.4))
  }

  # second
  if(!all(is.na(sv2))){
  rect(xleft=sv2[,3]-sv2[,2]*5000, xright=sv2[,3], ybottom=-1e6, ytop=1e6, border=NA, col=adjustcolor("deeppink", alpha=0.4))
  }

  ## Plot right histo
  par(fig=c(0.8,1, 0.3,1), new=T, mar=mar.two)
  plot(0, type="n", ylim=c(0,ceiling(max(allbreaks$depth))), xlim=c(0, max(table(histo))), xaxt="n", yaxt="n", xlab="", ylab="")
  grid()
  segments(x0=rep(0, length(breaks)), x1=table(histo), y0=breaks, y1=breaks, lwd=4)

  ## plot bottom part
  par(fig=c(0,0.8, 0,0.3), new=T, mar=mar.three)
  plot(discr, numhethom,type="l", lty=2, col="gray")
  points(discr, numhethom, pch=".", cex=2, col="blue3")
  points(discr, numhethom2, pch=".", cex=2, col="deeppink")
  points(discr, numhethom2, type="l", lty=2, col="deeppink")
  #dev.off()
}


get_stats=function(calls, allbreaks){
  cat(" [from script (get_stats):] Producing segment statistics...\n")

  chr=unique(calls$chr)[1]
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
  if(nrow(cnv_stats)>1){
  cnv_stats$chr=chr
  #write.table(cnv_stats, paste0(out, ".stats"), row.names=F, quote=F)
  }
  return(cnv_stats)
}


genotype_closest=function(alldepth, chrdp, stat){
  ### This function performs genotyping according to mean depth in the region, ignoring depth segments.
    ## say hello
    cat(" [from script (genotype_closest):] Genotyping using interval-wide means\n")
    ## prepare
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
        u=apply(dmin, 2, function(x) mean(as.numeric(as.character(x)), na.rm=T));
        topcp=ceiling(max(u/chrdp$V2, na.rm=T)*2)/2;
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


plot_cnv_qc=function(stats, outdf, plinkcall, alldepth, segments){
          chr=unique(alldepth$chr)
        alldepth$chr=NULL
        pos=alldepth$pos
        alldepth$pos=NULL
  # general plot cfg
  #png(paste(obase, "png", sep="."), width=1500, height=800)
  par(mfrow=c(2,2), mar = c(2,2,1,2))

  # plot 

  # plot mean-based
  # needed: attgeno, u/chrdp, dmin, topcp, start, stop
  apply(stats, 1, function(x){
  chr=x[1]
  start=x[2]
  stop=x[3]
  cat(" [from script (plot_cnv_qc):] Generating plots for interval ", chr, ":", start,"-", stop,"\n")

  outdf=outdf[outdf$start==start,]
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
  #lal=lal[1:3,]
  apply(lal, 1, function(x) {lines(rollmean(x, k))})
  #dev.off()

  ## third plot : comparison
  pmin=plinkcall[plinkcall$start==start,]
  outdf$ndp=outdf$region_depth/outdf$chrwide_depth
  m=merge(pmin, outdf, by="sample")
  relevant=segments[segments$min<x[3] & segments$max>x[2],]
  plot(m$ndp, m$adp, pch=m$copy_number, col=m$acall*2+1, xlab="average depth(method 1)", ylab="segment depth(method 2)",
  main="color=segment call, symbol=depth call", ylim=c(0, max(relevant$depth)+0.2))

  ## last plot: segments
  
  plot(0, type="n", xlim=c(start-(stop-start)*0.1, stop+(stop-start)*0.1), ylim=c(0, max(relevant$depth)+0.2), xlab="position", ylab="depth")
  rect(start, -1, stop, 2*max(relevant$depth), border=NA, col="lightgray")
  badsamples=pmin$sample[pmin$n>1]
  relevant$color="darkslategray"
  relevant$lwd=1
  relevant$color[relevant$sample %in% badsamples]="firebrick"
  relevant$lwd[relevant$sample %in% badsamples]=2
  relevant=relevant[order(relevant$color),]
  segments(y0=relevant$depth, x0=relevant$min, x1=relevant$max, col=relevant$color, lwd=relevant$lwd)
  })
}


genotype_segment=function(allbreaks, calls, include_dups){

  chr=unique(calls[,1])
  cat(" [from script (genotype_segment):] ")
  cat("Genotyping using constant segments (")
  cat(ifelse(include_dups, "including dups)\n", "coercing dups to 1)\n"))

  if(include_dups){
    cnv_df=NULL
    apply(calls, 1, function(x){
    relevant=allbreaks[allbreaks$min<x[3] & allbreaks$max>x[2],]
    # aggregation: for each sample, we compute the following vote:
    # relative length of depth segment x confidence x assigned depth
    relevant$umin=ifelse(relevant$min<x[2], x[2], relevant$min)
    relevant$umax=ifelse(relevant$max>x[3], x[3], relevant$max)
    relevant$length=abs(relevant$umax-relevant$umin)
    event_length=abs(x[3]-x[2])
    
    sm=NULL
    allsamples=unique(relevant$sample)
    for(s in allsamples){
    r=relevant[relevant$sample==s,]
    nr=nrow(r)
    quot=sum(r$length*r$assigned_depth)/(event_length)
    avgp=sum(r$length*r$p_value)/(event_length)
    #print(r)
    #print(c(nr, event_length, quot, avgp))
    quot=round(2*quot)/2
    sm=rbind(sm, c(s, quot, avgp))
    }
    sm=as.data.table(sm)
    colnames(sm)=c("sample", "depth", "confidence")
    sm$depth=as.numeric(sm$depth)
    sm$confidence=as.numeric(sm$confidence)
    u=data.frame(FID=allsamples, IID=allsamples, CHR=rep(chr, length(allsamples)), START=rep(x[2], length(allsamples)), END=rep(x[3], 
    length(allsamples)))
    u=merge(u, sm, all.x=T, by.x="FID", by.y="sample")
    u$DUMMY=rep(0, length(allsamples))
    cnv_df<<-rbind(cnv_df, u)
    })
    #cnv_df=do.call(rbind, cnv_df)
    write.table(cnv_df, paste0(out, ".cnv"), row.names=F, col.names=F, quote=F)
    fam_file=unique(cnv_df[,1])
    fam_file=cbind(fam_file, fam_file)
    fam_file=cbind(fam_file, rep(0, nrow(fam_file)), rep(0, nrow(fam_file)), rep(0, nrow(fam_file)), rep(-9, nrow(fam_file)))
    write.table(fam_file, paste0(out, ".fam"), row.names=F, col.names=F, quote=F)
  }else{
  # Exclude duplicates, force dups as 1
  # exclude dups, plink dataset
  alldepth=allbreaks
  alldepth$sample=as.character(alldepth$sample)
  samples=unique(alldepth$sample)
  geno=NULL
  # for every call
  forgetme=apply(as.matrix(calls), 1, function(x){

      # select only the depth segments that overlap with the region
      a=alldepth[alldepth$min<=x[3] & alldepth$max>=x[2],];
      a$min[a$min<x[2]]=x[2]
      a$min[a$max>x[3]]=x[3]
      # get list of samples that carry hom/het deletions
      o=NULL
      for(s in samples){
        sampledat=a[a$sample==s,]
        n=nrow(sampledat)
        sampledat$l=abs(sampledat$max-sampledat$min)
        dp=mean(sampledat$depth*sampledat$l)/sum(sampledat$l)
        p=mean(sampledat$p_value*sampledat$l)/sum(sampledat$l)
        adp=mean(sampledat$assigned_depth*sampledat$l)/sum(sampledat$l)
        topcp=ceiling(max(adp, na.rm=T)*2)/2;
        possibles=seq(0, topcp, by=0.5)
        acall=sapply(adp, function(x, possibles){v=abs(x-possibles);return(possibles[v==min(v)])}, possibles=possibles)[1]
        o=rbind(o, data.frame(sample=s, chr=x[1], start=x[2], stop=x[3], n=n, dp=dp, adp=adp, acall=acall,p=p))
      }
      return(o)
  })
  forgetme=do.call(rbind, forgetme)
  forgetme$sample=as.character(forgetme$sample)
  return(forgetme)
}

}


call_sv=function(dfcall, chrwide_depth, cparm){
      
    #lastcol=ncol(dfcall)
    #dfcall[,lastcol]=NULL
    nms=colnames(dfcall)
    #dfcall=dfcall[,(nms[length(nms)]):=NULL]
    #dfcall=as.data.frame(dfcall)
    cat(" [from script (call_sv):] Read Genotypes for ", ncol(dfcall)-2, "samples.\n")
    avgmeans_region=as.numeric(sapply(dfcall[,3:ncol(dfcall)], mean))
    
    snames=chrwide_depth$V1
    chrwide_depth=chrwide_depth$V2
    names(chrwide_depth)=snames
    colnames(dfcall)=c("chr", "pos", snames)
    cat(" [from script (call_sv):] Read per-sample depths.\n") 
    if(length(avgmeans_region)!=length(chrwide_depth)){
  stop("Error : number of samples different from number of depth columns.")
    }
    reldepth2 = avgmeans_region/chrwide_depth
    dfcall=as.data.frame(dfcall)
    returndf=NULL
    i=1
  cat(paste(" [from script (call_sv):] Using complexity parameter", cparm, "\n"))
    rpc=rpart.control( minbucket = 5, cp = cparm, 
              maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10,
              surrogatestyle = 0, maxdepth = 30)
    for (selected_sample in snames){
        if(selected_sample %in% colnames(dfcall)){
            cat(paste0(selected_sample, " (", i, "/", length(snames), ")\r"))
                        i=i+1
            flush.console()
        regiondata=data.frame(pos=dfcall$pos, dp=dfcall[,selected_sample]/chrwide_depth[selected_sample])

        tree <- rpart(dp ~ pos, data=regiondata, control=rpc)
        x=regiondata$pos
        s <- seq(min(x), max(x), by=100)
        pred=predict(tree, data.frame(pos=s))
        downsampled=data.frame(pos=seq(min(x), max(x), by=100), pred=pred)
        ret=cbind(aggregate(downsampled$pos, by=list(downsampled$pred), min), aggregate(downsampled$pos, by=list(downsampled$pred), max)$x)
        colnames(ret)=c("depth", "min","max")
        ret$sample=selected_sample
        ret=ret[,c("sample", "depth", "min","max")]
        returndf=rbind(returndf, ret)
        }
    }
    #write.table(returndf, paste(out, cparm, '.segments', sep="."), quote=F, row.names=F)
    return(returndf)
}


cat(" [from script (main):] Reading files.\n")

#freal="1.2-3M.alldp";chrdp="../average_chrdps/chr1.avgdp";out="test"
suppressMessages(library(data.table));
suppressMessages(library(mixtools)); 
suppressMessages(library(zoo))
suppressMessages(library(rpart))

## main depth per genotype file (longest to load)
dfcall = fread(freal, header=T)
#colnames(dfcall)[c(1,2)]=c("chr", "pos")
#print(colnames(dfcall))
## per-chromosome depth per sample
chrwide_depth=read.table(chrdp, stringsAsFactors=F, header=F)

## perform piecewise constant regression
allbreaks=call_sv(dfcall, chrwide_depth, complexity)
alldepth=allbreaks
confidence=0.05


## Calling
mvect=seq(0, ceiling(max(allbreaks$depth, na.rm=T)), by=0.5)
mix=normalmixEM(allbreaks$depth, k=length(mvect), mu=mvect, sigma=1, mean.constr=mvect, sd.constr=rep('a', length(mvect)))
parm=cbind(mix$mu, mix$sigma)


# discretisation
discr=seq(min(alldepth$min), max(alldepth$max), by=5000)
# for every depth measure, assign it a depth c(depth, prob, n)
callz=t(sapply(allbreaks$depth, function(meas) {
    # closest ideal depths
    botparm=parm[parm[,1]==floor(meas*2)/2,]
    topparm=parm[parm[,1]==ceiling(meas*2)/2,]
    # compute highest p
    plo=pnorm(meas, mean=botparm[1], sd=botparm[2], lower.tail=F)
    phi=pnorm(meas, mean=topparm[1], sd=topparm[2], lower.tail=T)
    topp=max(plo, phi)
    # deduce best ideal depth
    topn=ifelse(plo>phi, botparm[1], topparm[1])
    return(c(meas, 2*topp, topn))
}))

allbreaks$assigned_depth=callz[,3]
allbreaks$p_value=callz[,2]

# we don't want to call breaks that are greater (or close to) one in depth
if(include_dups){
tocall=allbreaks[allbreaks$assigned_depth!=1, ] 
}else{
tocall=allbreaks[allbreaks$assigned_depth<1, ] 
}

#this is our discretisation grid on the x axis
discr=seq(min(alldepth$min), max(alldepth$max), by=5000)

# these two functions are indicators of presence, along
# the grid, of depths different than one

# we use two indicators: numhethom
# calculates the number of high-confidence depth segments in the region
numhethom=sapply(discr, function(x){
    in_region=tocall[tocall$min<=x & tocall$max>=x,]
    nrow(in_region[in_region$p_value>confidence,])
    })

# numhethom2 is the ratio of high vs low confidence in the region
numhethom2=sapply(discr, function(x){
    in_region=tocall[tocall$min<=x & tocall$max>=x,]
    if(nrow(in_region)==0){return(0)}
    adjust_factor=nrow(in_region[in_region$p_value<confidence,])
    if(!is.finite(adjust_factor)){adjust_factor=1}
    nrow(in_region[in_region$p_value>confidence,])/adjust_factor
    })

##### REGION CALLING 1 : 
# set of longest continuous lengths where count of depths > 0
a=rle(numhethom>0) ## arbitrary threshold here we come
svs=t(rbind(a$values, a$lengths, min(discr)+cumsum(a$lengths)*5000))
svs=svs[svs[,1]==1,,drop=F]
sv1=svs

##### REGION CALLING 2 :
a=rle(numhethom2>1) ## arbitrary threshold here we come
svs=t(rbind(a$values, a$lengths, min(discr)+cumsum(a$lengths)*5000))
svs=svs[svs[,1]==1,,drop=F]
sv2=svs


calls=svs
calls[,2]=svs[,3]-svs[,2]*5000

## At this point, the SVs are called using the two methods, but not genotyped.

## We plot the CNVs in the region
pdf(paste(out, "window.pdf", sep="."), width=15)
plot_window(allbreaks, confidence, sv1, sv2, discr, numhethom, numhethom2)
dev.off()

if(nogeno){
cat(" [from script (main):] No genotyping requested, exiting.\n")
}else{

## Compute statistics for SVs
colnames(calls)=c("chr", "start", "stop")
calls=as.data.frame(calls)
cnv_stats=get_stats(calls, allbreaks)

## genotyping method 1
segm=genotype_segment(allbreaks, calls, include_dups)

## genotyping method 2 (means)
closest=genotype_closest(dfcall, chrwide_depth, cnv_stats)


## generating pdf of stats for all called CNVs
pdf(paste(out, "QC.pdf", sep="."), width=10)
plot_cnv_qc(cnv_stats, closest, segm, dfcall, allbreaks)
dev.off()

## merging and writing
m=merge(segm, closest, by=c("sample", "start"))
m$chr.y=NULL;m$end=NULL;
m$chr=m$chr.x;m$chr.x=NULL;m=m[,c("sample", "chr", "start", "stop", "n", "dp", "adp", "acall", "p", "chrwide_depth", "region_depth", "copy_number", "inverse_confidence")]
m$copy_number=m$copy_number/2
colnames(m)[5:13]=c("segments", "segment_weighted_depth", "segment_weighted_assigned_depth", "segment_total_assigned_depth", "segment_weighted_pvalue", "chrwide_depth", 
  "region_depth", "means_assigned_depth", "means_inverse_confidence")
write.table(m, paste(out, "genotypes.csv", sep="."), sep=",", quote=F, row.names=F)

write_plink(m, out, method, include_dups)
}
cat(" [from script (main):] End of script.\n")

