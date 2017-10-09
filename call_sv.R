#!/software/bin/Rscript --vanilla

argv <- commandArgs(TRUE)
library(data.table)
call_sv=function(filename, sampledepth, out){
    library(rpart)
    dfcall = fread(filename, na.strings=".")
    lastcol=ncol(dfcall)
    dfcall[,lastcol]=NULL
    print("Read Genotypes.")
    avgmeans_region=as.numeric(sapply(dfcall[,3:ncol(dfcall)], mean))
    chrwide_depth=read.table(sampledepth, stringsAsFactors=F)
    snames=chrwide_depth$V1
    chrwide_depth=chrwide_depth$V2
    names(chrwide_depth)=snames
    colnames(dfcall)=c("chr", "pos", snames)
    print("Read per-sample depths.") 
    reldepth2 = avgmeans_region/chrwide_depth
    dfcall=as.data.frame(dfcall)
    returndf=NULL
    i=1
    for (selected_sample in snames){
        if(selected_sample %in% colnames(dfcall)){
            cat(paste0(selected_sample, " (", i, "/", length(snames), ")\r"))
                        i=i+1
            flush.console()
        regiondata=data.frame(pos=dfcall$pos, dp=dfcall[,selected_sample]/chrwide_depth[selected_sample])

        tree <- rpart(dp ~ pos, data=regiondata)
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
    write.table(returndf, out, quote=F, row.names=F)
}

call_sv(argv[1], argv[2], argv[3])
