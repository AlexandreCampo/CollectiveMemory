
# read a results file
#m = read.table ("res.txt")

# get the second part of the exp and round time
#m = subset (m, m[, 2] >= 2500)
#m[,2] = round (m[,2], 0)


library(gridBase) ;
library(gplots) ;
library(plotrix) ;
library (boot) ;


medianbootstrap = function (x, idx) {return(median(x[idx] + rnorm(1,sd=1/sqrt(length(x)))))}
bootR = 999
bootCIL = 25
bootCIU = 975

myBlue = rgb(red=0,blue=1,green=0,alpha=0.35)
myBlueLight = rgb(red=0,blue=1,green=0,alpha=0.17)
myRed = rgb(red=1,blue=0,green=0,alpha=0.99)
myRedLight = rgb(red=1,blue=0,green=0,alpha=0.25)

Bootstrap = function (m)
{
  med = vector("numeric", ncol(m))
  cil = vector("numeric", ncol(m))
  ciu = vector("numeric", ncol(m))
  for (i in 1:ncol(m))
    {
      print (paste ("Bootstraping time", i, "/", ncol(m)))
      
      # do a bootstrap
      x = m[,i]
      b = boot (x, medianbootstrap, R = bootR)
      o = order(b$t)
      d = b$t[o]
    
      # record median and CI
      med[i] = median (d)
      cil[i] = d[bootCIL]
      ciu[i] = d[bootCIU]
    }
  return (list(med,cil,ciu))
}



PlotAllChoices = function ()
  {
    resfiles = system(paste("ls results*.txt", sep=""), T)

    synth = NULL
    
    # for each resfile
    for (r in resfiles)
      {
        print (paste("Taking care of ", r))        

        # plot in time        
        output = system(paste("basename ", r, " .txt"), T)
        PlotChoice (r, output)
      }
  }


PlotChoice = function (filename, output = "fig")
{ 
    # format : replication / time / opinion of each robot
    m = read.table (filename)    

    m = subset (m, m[, 2] >= 2500)
    m[,2] = round (m[,2], 0) - 2500
    
    totalpop = 30
    
    # in quorum, gather the proportion of pop that believes recalling
    n = cbind(m[,1:3], totalpop - m[,3])
    
    # separate curves for 2 choices
    time = unique(n[,2])
    maxtime = time[length(time)]

    ends = subset(n, n[,2] == maxtime)

    repsA = subset(ends, ends[,3] >= 0.5 * totalpop)[,1]
    repsB = subset(ends, ends[,3] < 0.5 * totalpop)[,1]
            

    C = NULL
    NC = NULL
    for (i in repsA)      
    {
      C = rbind (C, subset(n, n[,1] == i)[,3])
      NC = rbind (NC, subset(n, n[,1] == i)[,4])
    }
    for (i in repsB)      
    {
      C = rbind (C, subset(n, n[,1] == i)[,4])
      NC = rbind (NC, subset(n, n[,1] == i)[,3])
    }
    
  
  # subsample, otherwise it is ugly
  x = length(time)
  step = x / 15
    
  ns = c(seq(1,x,step), x)
  time = time[ns]
  C = C[,ns]
  NC = NC[,ns]
    
    medC = vector("numeric", length(time))
    cilC = vector("numeric", length(time))
    ciuC = vector("numeric", length(time))
    
    if (nrow(C) > 0)
      {
        bs = Bootstrap (C)
        medC = bs[[1]] / totalpop
        cilC = bs[[2]] / totalpop
        ciuC = bs[[3]] / totalpop
      }
    
    medNC = vector("numeric", length(time))
    cilNC = vector("numeric", length(time))
    ciuNC = vector("numeric", length(time))
    
    if (nrow(NC) > 0)
      {
        bs = Bootstrap (NC)
        medNC = bs[[1]] / totalpop
        cilNC = bs[[2]] / totalpop
        ciuNC = bs[[3]] / totalpop
      }

    pdf(paste(output, "Dynamics.pdf", sep=""),width=6,height=6) ;    
    par(mar=c(5,5,2,2)+0.1) ;
    
    plot( medC~time, type="l", ylim = c(0,1), xlim=c(0, 2500), col=myRed, lwd=3, cex.axis=1.5, xlab="Time (s)", ylab="Fraction of population", cex.lab=1.5, main=output) ;
    polygon( c(cilC,rev(ciuC)) ~ c(time, rev(time) + 1), border = NA, col = myRedLight) ;

    lines( medNC ~ time, type="l", ylim = c(0,1), xlim=c(0, 2500), col=myBlue, lwd=3, cex.axis=1.5, cex.lab=1.5) ;
    polygon( c(cilNC,rev(ciuNC)) ~ c(time, rev(time) + 1), border = NA, col = myBlueLight) ;
    
    legend(x="right",legend=c("Winning opinion","Discarded opinion"),fill=c(myRed,myBlue),bg="white",cex=1.5,bty="n") ;
        
    dev.off() ;  


    png(paste(output, "Dynamics.png", sep=""),width=800,height=800)
    par(mar=c(7,7.5,1,1)+0.1)
    par(mgp = c(5,2,0))
    
    plot( medC~time, type="l", ylim = c(0,1), xlim=c(0, 2500), col=myRed, lwd=3, cex.axis=3, xlab="Time (s)", ylab="Fraction of population", cex.lab=3, main=output) ;
    polygon( c(cilC,rev(ciuC)) ~ c(time, rev(time) + 1), border = NA, col = myRedLight) ;

    lines( medNC ~ time, type="l", ylim = c(0,1), xlim=c(0, 2500), col=myBlue, lwd=3) ;
    polygon( c(cilNC,rev(ciuNC)) ~ c(time, rev(time) + 1), border = NA, col = myBlueLight) ;
    
    legend(x="right",legend=c("Winning opinion","Discarded opinion"),fill=c(myRed,myBlue),bg="white",cex=3,bty="n") ;
        
    dev.off() ;  

  }


PlotAllHists = function ()
  {
    resfiles = system(paste("ls results*.txt", sep=""), T)

    props = NULL
    cil = NULL
    ciu = NULL
    for (f in resfiles)
      {
        print (paste("Taking care of ", f))
        output = system(paste("basename ", f, " .txt"), T)
        r = PlotHist (f, output)        
        props = rbind(props, r[[1]])
        cil = rbind(cil, r[[2]])
        ciu = rbind(ciu, r[[3]])
      }

    # plot differences
    pdf ("threshold.pdf", useDingbats=F)
    len = nrow(props)
    halflen = len / 2
    idx1 = 1:halflen
    idx2 = (halflen+1):len

    x = seq(0, 0.60, 0.05)
    
    plot(x, props[idx1,1], type = "l", xlab = "Recall threshold", ylab = "Proportion of false positives", ylim = c(0.0,1), col = myRed, cex.lab = 1.6, lwd=2.5, las=1, cex.axis=1.6)
    polygon( c(x, rev(x)), c(cil[idx1,1], rev(ciu[idx1,1])), col = myRedLight, border = NA)
    lines(x, props[idx2,2], type = "l", col = myBlue, lwd=2.5)
    polygon( c(x, rev(x)), c(cil[idx2,2], rev(ciu[idx2,2])), col = myBlueLight, border = NA)
    legend(x="topright",legend=c("Failed recalls", "Wrong recalls"),fill=c(myRed,myBlue),bg="white",cex=3,bty="n") ;
    dev.off()

    png ("threshold.png",  width=800, height=800)
    par(mar=c(7.5,8.8,1,1)+0.1)
    par(mgp = c(6.5,2,0))

    len = nrow(props)
    plot(x, props[idx1,1], type = "l", xlab = "Recall threshold", ylab = "Proportion of false positives", ylim = c(0.0,1), col = myRed, cex.lab = 3, lwd=3, las=1, cex.axis=3)
    polygon( c(x, rev(x)), c(cil[idx1,1], rev(ciu[idx1,1])), col = myRedLight, border = NA)
    lines(x, props[idx2,2], type = "l", col = myBlue, lwd=2.5)
    polygon( c(x, rev(x)), c(cil[idx2,2], rev(ciu[idx2,2])), col = myBlueLight, border = NA)
    legend(x="topright",legend=c("Failed recalls", "Wrong recalls"),fill=c(myRed,myBlue),bg="white",cex=3,bty="n") ;
    dev.off()
    
return(list(props, idx1, idx2))
  }



PlotHist = function (filename, output = "fig")
{
    # format : replication / time / opinion of each robot
    m = read.table (filename)    

    m = subset (m, m[, 2] >= 2500)
    m[,2] = round (m[,2], 0) - 2500
    
    totalpop = 30
    
    n = cbind(m[,1:3], totalpop - m[,3])

    
    # now separate curves for 2 choices
    time = unique(n[,2])
    maxtime = time[length(time)]
    ends = subset(n, n[,2] == maxtime)

    repsA = subset(ends, ends[,3] >= 0.5 * totalpop)[,1]
    repsB = subset(ends, ends[,4] >= 0.5 * totalpop)[,1]

    S1 = length(repsA)
    S2 = length(repsB)
    
  res1 = binom.test(S1,S1+S2) ;
  res2 = binom.test(S2,S1+S2) ;

  print (res1)
  
    tmp = c(S1 / (S1+S2), S2 / (S1+S2))
  
  # ci formula : sqrt ( p * (1 - p) / n) * 1.96
  tmpci = sqrt (tmp * (1 - tmp) / (S1+S2)) * 1.96
  tmpciu = sapply (tmp + tmpci, function(x){return (min(x,1))})
  tmpcil = sapply (tmp - tmpci, function(x){return (max(x,0))})
  
  res = matrix(tmp, ncol=2, byrow=FALSE);
  ci.l = matrix(tmpcil,ncol=2,byrow=FALSE) ;
  ci.u = matrix(tmpciu,ncol=2,byrow=FALSE) ;
  

  pdf(paste(output, "Hist.pdf", sep=""),width=6,height=6) ;  
  par(xpd=T,mar=c(5,5,2,2)+0.1) ;  
  barplot2(res,plot.ci=TRUE,ci.l=ci.l,ci.u=ci.u,ci.lwd=3,beside=TRUE,cex.axis=1.5,cex.lab=1.5,
           ylim=c(-0.01,1.01),plot.grid=T,cex.names=1.5,grid.lwd=3,lwd=3, las=1,
           col=c(myRed,myBlue),
           ylab="Proportion of trials", xlab="Collective choice", names.arg=c("Opinion 1", "Opinion 2"));
  dev.off() ;

  png(paste(output, "Hist.png", sep=""),width=800,height=800) 
    par(mar=c(7,7.5,1,1)+0.1)
    par(mgp = c(5,2,0))
  barplot2(res,plot.ci=TRUE,ci.l=ci.l,ci.u=ci.u,ci.lwd=3,beside=TRUE,cex.axis=3,cex.lab=3,
           ylim=c(-0.01,1.01),plot.grid=T,cex.names=3,grid.lwd=3,lwd=3, las=1,
           col=c(myRed,myBlue),
           ylab="Proportion of trials", xlab="Collective choice", names.arg=c("Opinion 1", "Opinion 2"));
  dev.off() ;

  return (list(res, ci.l, ci.u))
}

