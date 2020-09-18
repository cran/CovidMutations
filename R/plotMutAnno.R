#' Plot the mutation statistics after annotating the "nucmer" object by "indelSNP" function
#'
#' @description Basic descriptions for the mutational events.
#'
#' @param covid_annot The mutation effects provided by "indelSNP" function.
#' @param figureType Figure type for: "MostMut", "MutPerSample", "VarClasses", "VarType", "NucleoEvents", "ProEvents".
#' @param outdir The output directory.
#'
#' @return Plot the selected figure type as output.
#' @export
#' @importFrom graphics par barplot
#' @importFrom grDevices png dev.off heat.colors terrain.colors cm.colors topo.colors rainbow
#' @examples
#' data("covid_annot")
#' # make sure the covid_annot is a dataframe
#' covid_annot <- as.data.frame(covid_annot)
#' #outdir <- tempdir() specify your output directory
#' plotMutAnno(covid_annot = covid_annot, figureType = "MostMut", outdir = NULL)
plotMutAnno <- function(covid_annot = covid_annot, figureType = "MostMut", outdir = NULL){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))# resetting the par when exit function
    if(figureType == "MostMut"){
      if(is.null(outdir) == FALSE){
        figurePath<- file.path(outdir,paste("MostMut_","covid_annot.png", sep = ""),sep = "")
        png(figurePath,width = 3000,height = 2000,res=300)
        # Most mutated samples
        occ<-sort(table(covid_annot$sample),decreasing=TRUE)[1:20]
        par(las=2,mar=c(15,5,5,1))
        barplot(occ,ylab="nr of mutations",main="Most mutated samples",col=topo.colors(length(occ)))
        dev.off()
      }
      if(is.null(outdir) == TRUE){
        occ<-sort(table(covid_annot$sample),decreasing=TRUE)[1:20]
        barplot(occ,ylab="nr of mutations",main="Most mutated samples",col=topo.colors(length(occ)))
      }

    }
    if(figureType == "MutPerSample"){
      if(is.null(outdir) == FALSE){
        figurePath<- file.path(outdir,paste("MutPerSample_","covid_annot.png",sep = ""),sep = "")
        png(figurePath,width = 3000,height = 2000,res=300)
        occ<-table(table(covid_annot$sample))
        par(las=2,mar=c(5,5,5,1))
        barplot(occ,xlab="nr of mutations",main="Overall mutations per sample",col=topo.colors(length(occ)))
        dev.off()
      }
      if(is.null(outdir) == TRUE){
        occ<-table(table(covid_annot$sample))
        barplot(occ,xlab="nr of mutations",main="Overall mutations per sample",col=topo.colors(length(occ)))
      }

    }
    if(figureType == "VarClasses"){
      if(is.null(outdir) == FALSE){
        figurePath<- file.path(outdir,paste("VarClasses_","covid_annot.png",sep = ""),sep = "")
        png(figurePath,width = 3000,height = 2000,res=300)
        occ<-sort(table(covid_annot$varclass),decreasing=TRUE)
        par(las=2,mar=c(8,5,5,1))
        barplot(occ,ylab="nr of events",main="Most frequent events per class",col=cm.colors(length(occ)))
        dev.off()
      }
      if(is.null(outdir) == TRUE){
        occ<-sort(table(covid_annot$varclass),decreasing=TRUE)
        barplot(occ,ylab="nr of events",main="Most frequent events per class",col=cm.colors(length(occ)))
      }

    }
    if(figureType == "VarType"){
      if(is.null(outdir) == FALSE){
        figurePath<- file.path(outdir,paste("VarType_","covid_annot.png",sep = ""),sep = "")
        png(figurePath,width = 3000,height = 2000,res=300)
        occ<-sort(table(apply(covid_annot[,c("refvar","qvar")],1,paste0,collapse=">")),decreasing=TRUE)[1:20]
        par(las=2,mar=c(5,5,5,1))
        barplot(occ,ylab="nr of samples",main="Most frequent per type",col=cm.colors(length(occ)))
        dev.off()
      }
      if(is.null(outdir) == TRUE){
        occ<-sort(table(apply(covid_annot[,c("refvar","qvar")],1,paste0,collapse=">")),decreasing=TRUE)[1:20]
        barplot(occ,ylab="nr of samples",main="Most frequent per type",col=cm.colors(length(occ)))
      }

    }
    if(figureType == "NucleoEvents"){
      if(is.null(outdir) == FALSE){
        figurePath<- file.path(outdir,paste("NucleoEvents_","covid_annot.png",sep = ""),sep = "")
        png(figurePath,width = 3000,height = 2000,res=300)
        occ<-sort(table(apply(covid_annot[,c("refvar","refpos","qvar")],1,paste0,collapse="")),decreasing=TRUE)[1:10]
        par(las=2,mar=c(8,5,5,1))
        barplot(occ,ylab="nr of samples",main="Most frequent events (nucleotide)",col=heat.colors(length(occ)))
        dev.off()
      }
      if(is.null(outdir) == TRUE){
        occ<-sort(table(apply(covid_annot[,c("refvar","refpos","qvar")],1,paste0,collapse="")),decreasing=TRUE)[1:10]
        barplot(occ,ylab="nr of samples",main="Most frequent events (nucleotide)",col=heat.colors(length(occ)))
      }

    }
    if(figureType == "ProEvents"){
      if(is.null(outdir) == FALSE){
        figurePath<- file.path(outdir,paste("ProEvents_","covid_annot.png",sep = ""),sep = "")
        png(figurePath,width = 3000,height = 2000,res=300)
        occ<-sort(table(apply(covid_annot[,c("protein","variant")],1,paste0,collapse=":")),decreasing=TRUE)[1:10]
        par(las=2,mar=c(8,5,5,1))
        barplot(occ,ylab="nr of samples",main="Most frequent events (protein)",col=rainbow(length(occ)))
        dev.off()
      }
      if(is.null(outdir) == TRUE){
        occ<-sort(table(apply(covid_annot[,c("protein","variant")],1,paste0,collapse=":")),decreasing=TRUE)[1:10]
        barplot(occ,ylab="nr of samples",main="Most frequent events (protein)",col=rainbow(length(occ)))
      }

    }


}
