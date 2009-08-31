`bdoc` <-
function(traindata,testdata,delta=9.7e-8,epsilon=.2,priors="equal",stoppingrule=TRUE,impute=1,plot.file="pdf"){
    
    date1 <- Sys.time()
    if(delta<=0|delta>0.1) stop("delta must be strictly greater than 0 and less than or equal to 0.1")
    if(epsilon<0|epsilon>1) stop("epsilon must be greater than or equal to 0 and less than or equal to 1")
    testk = nrow(testdata)
    natodash<-function(x){
y<-ifelse(is.na(x)==TRUE,"-",x)
    }
    traindata<-apply(traindata,2,natodash)
    testk<-nrow(testdata)
    testdata<-cbind(rep("test",testk),rep("test",testk),testdata)
    colnames(testdata)<-colnames(traindata)

    testdata<-apply(testdata,2,natodash)
   #testp<-ncol(testdata)-2
    ifelse(is.vector(testdata)==TRUE,testp<-(length(testdata)-2),testp<-(ncol(testdata)-2))

    bases = 4
    set <- by(as.data.frame(traindata), as.data.frame(traindata)$species, subset)
    set2 <- list()
    s.num <- nrow(set)
    spec <- list()
    likelihood <- array(0, c(dim(set), bases, testp))
    dirichlet.prior<-function(alpha,scale=1){
        x<-rgamma(length(alpha),alpha,scale) #produces B x1-xk variates from gamma with common scale and shape=alpha1,...,alphak
        theta<-x/sum(x)         #computes theta=x/sum(x)  
        return<-list(theta=theta,x=x,alpha=alpha,scale=scale)
        }
    majrule <- function(x) {
        y1 <- names(which.max(table(ifelse(x != "A", ifelse(x != 
            "T", ifelse(x != "C", ifelse(x != "G", NA, "G"), 
            "C"), "T"), "A"))))
        if (identical(is.null(y1), TRUE) == TRUE) {
            y <- rep(NA, length(x))
        }
        else {
            y <- ifelse(x != "A", ifelse(x != "T", ifelse(x != 
                "C", ifelse(x != "G", y1, "G"), "C"), "T"), "A")
            y <- ifelse(is.na(x) == TRUE, y1, y)
        }
    }
    proprule <- function(x) {
        y1 <- ifelse(x != "A", ifelse(x != "T", ifelse(x != "C", 
            ifelse(x != "G", NA, "G"), "C"), "T"), "A")
        tab <- table(y1)
        if (identical(is.null(names(tab)), TRUE) == TRUE) {
            y <- rep(NA, length(x))
        }
        else {
            y <- y1
            for (i in 1:length(x)) {
                if (is.na(y1[i]) == TRUE) {
                  y[i] <- sample(names(tab), 1, prob = tab/sum(tab))
                }
            }
        }
        y
    }
    date15<-Sys.time()
    for (i in 1:s.num) {
        if (impute == 2) {
            if (nrow(set[[i]]) == 1) {
                set2[[i]] <- cbind(set[[i]][, 1:2], t(as.matrix(apply(set[[i]][, 
                  3:(testp + 2)], 2, majrule))))
            }
            else {
                set2[[i]] <- cbind(set[[i]][, 1:2], apply(set[[i]][, 
                  3:(testp + 2)], 2, majrule))
            }
        }
        else {
            if (nrow(set[[i]]) == 1) {
                set2[[i]] <- cbind(set[[i]][, 1:2], t(as.matrix(apply(set[[i]][, 
                  3:(testp + 2)], 2, proprule))))
            }
            else {
                set2[[i]] <- cbind(set[[i]][, 1:2], apply(set[[i]][, 
                  3:(testp + 2)], 2, proprule))
            }
        }
    }
    date16<-Sys.time()
    for (i in 1:s.num) {
        spec[[i]] <- as.list(apply(set2[[i]], 2, table))
        if (is.null(dim(spec[[i]][[1]])) == TRUE) {
            for (k in 3:(testp + 2)) {
                names(spec[[i]][[k]]) <- set2[[i]][[k]][[1]]
            }
        }
        for (j in 3:(testp + 2)) {
            likelihood[i, 1, j - 2] <- sort(ifelse(names(spec[[i]][[j]][1:length(spec[[i]][[j]])]) == 
                "A", spec[[i]][[j]][1:length(spec[[i]][[j]])]/sum(spec[[i]][[j]]), 
                0), decreasing = TRUE)[1]
            likelihood[i, 2, j - 2] <- sort(ifelse(names(spec[[i]][[j]][1:length(spec[[i]][[j]])]) == 
                "T", spec[[i]][[j]][1:length(spec[[i]][[j]])]/sum(spec[[i]][[j]]), 
                0), decreasing = TRUE)[1]
            likelihood[i, 3, j - 2] <- sort(ifelse(names(spec[[i]][[j]][1:length(spec[[i]][[j]])]) == 
                "C", spec[[i]][[j]][1:length(spec[[i]][[j]])]/sum(spec[[i]][[j]]), 
                0), decreasing = TRUE)[1]
            likelihood[i, 4, j - 2] <- sort(ifelse(names(spec[[i]][[j]][1:length(spec[[i]][[j]])]) == 
                "G", spec[[i]][[j]][1:length(spec[[i]][[j]])]/sum(spec[[i]][[j]]), 
                0), decreasing = TRUE)[1]
        }
    }
    date17<-Sys.time()
    adjlikelihood <- array(0, c(s.num, bases, testp))
    for (j in 1:s.num) {
        for (i in 1:testp) {
            adjlikelihood[j, , i] <- ifelse(likelihood[j, , i] == 
                0, delta, likelihood[j, , i] - (sum(ifelse(likelihood[j, 
                , i] == 0, delta, likelihood[j, , i])) - 1)/(bases - 
                (1/delta) * (sum(ifelse(likelihood[j, , i] == 
                  0, delta, likelihood[j, , i])) - 1)))
        }
    }
adjlikelihood2 <- array(0, c(s.num, bases, testp))
    for (j in 1:s.num) {
        for (i in 1:testp) {
            adjlikelihood2[j, , i] <- ifelse(is.na(adjlikelihood[j,,i])==TRUE,1/bases, adjlikelihood[j, , i])
        }
    }
    date18<-Sys.time()
    speciestots <- by(as.data.frame(traindata), as.data.frame(traindata)$species, nrow)
    if(length(priors)!=1){
       if (all.equal(sum(priors), 1) != TRUE) 
            warning("Improper priors specified.  Priors have been rescaled to sum to 1")
        priors<-priors/sum(priors)
        names(priors) <- names(speciestots)
        priors <- by(priors, names(priors), sort)
    }
    else{
     if (priors == "equal") {
        priors<-rep(1/s.num,s.num)
         names(priors)<-names(speciestots)
        priors<-by(priors,names(priors),sort)
    }
    else if(priors=="data"){
        priors <- speciestots/nrow(traindata)

    }
    else if(priors=="dir"){
priors <- dirichlet.prior(rep(1,s.num))$theta
      names(priors)<-names(speciestots)
     priors<-by(priors,names(priors),sort)

    }
    else{
       priors<-rep(1/s.num,s.num)
         names(priors)<-names(speciestots)
        priors<-by(priors,names(priors),sort)

    }
    }
    date3 <- Sys.time()
    post.mat <- matrix(data = 0, s.num, testp)
    post <- matrix(data = 0, s.num, bases)
    tot <- matrix(data = 0, bases, s.num)
    testpriors <- priors
    testlikelihood <- adjlikelihood2
    posteriors <- list()
    species.class <- matrix(data = 0, 2, testk)
    rownames(species.class) <- c("Species", "Prob")
    colnames(species.class) <- c(1:testk)
    unit <- vector()
    unit2 <- vector()
    trC<-as.numeric(ifelse(as.vector(as.matrix(traindata[,-c(1:2)]))=="A",65,ifelse(as.vector(as.matrix(traindata[,-c(1:2)]))=="T",84,ifelse(as.vector(as.matrix(traindata[,-c(1:2)]))=="C",67,ifelse(as.vector(as.matrix(traindata[,-c(1:2)]))=="G",71,0)))))
    identical.indicator<-vector()
    date14<-Sys.time()

date10<-Sys.time()
    dyn.load(paste("bdoc2", .Platform$dynlib.ext, sep = ""))
    dyn.load(paste("bdoc1", .Platform$dynlib.ext, sep = ""))
    for (k in 1:testk) {
        #bp <- testdata[k, 3:ncol(testdata)]
        ifelse(is.vector(testdata)==TRUE,bp<-testdata[3:length(testdata)],bp <- testdata[k, 3:ncol(testdata)])
        bpC<-as.numeric(ifelse(bp=="A",65,ifelse(bp=="T",84,ifelse(bp=="C",67,ifelse(bp=="G",71,0)))))
        
       
       if(stoppingrule==TRUE){
       res<-.C("bdoc",as.integer(s.num),as.integer(bases),as.integer(testp),
              as.integer(nrow(traindata)),as.double(epsilon),as.double(testpriors),
              as.numeric(adjlikelihood2),as.integer(bpC),
              as.integer(trC),
              posteriors=numeric(s.num*testp),
              q=as.integer(0),
              bar2=numeric(testp),PACKAGE="bdoc2")
          

          post.mat<-matrix(res$posteriors,s.num)
          j<-res$q
        }
       else{
         res<-.C("bdoc1",as.integer(s.num),as.integer(bases),as.integer(testp),
              as.integer(nrow(traindata)),as.double(epsilon),as.double(testpriors),
              as.numeric(adjlikelihood2),as.integer(bpC),
              as.integer(trC),
              posteriors=numeric(s.num*testp),
              q=as.integer(0),
              bar2=numeric(testp),PACKAGE="bdoc1")
          

          post.mat<-matrix(res$posteriors,s.num)
          j<-res$q
        }
        
        
        identical.post<-which(post.mat[,j]==max(post.mat[,j]))
        identical.indicator[k]<-ifelse(length(identical.post)==1,FALSE,TRUE)
        species.id <- cbind(c(rownames(testpriors)[which.max(post.mat[, j])], max(post.mat[, j])))
        post.mat.names<-vector()
        for(q in 1:j){
	     post.mat.names[q]<-paste("Position",q)
        }
        posteriors[[k]] <- list(species.id = species.id, post = post.mat[, 
            1:j])
        rownames(posteriors[[k]]$post)<-names(speciestots)
        colnames(posteriors[[k]]$post)<-post.mat.names
        species.class[, k] <- posteriors[[k]]$species.id
        testpriors <- priors
    }

date11<-Sys.time()

    dyn.unload(paste("bdoc2", .Platform$dynlib.ext, sep = ""))
    dyn.unload(paste("bdoc1", .Platform$dynlib.ext, sep = ""))

    date2 <- Sys.time()
    tottime <- date2 - date1
    imptime<-date16-date15
    liketime<-date18-date16
    classtime<-date11-date14
    
    ######################################################################################################
    ################ Creates and saves plots of the posteriors as pdf, jpg, png, wmf, or ps.
    ################
    ######################################################################################################
par(oma=c(1,0,0,0))   #This expands the bottom margin a bit to get the message about what noisy plots indicate
posteriors2<-matrix() #############################################################################################
      jitter2<-function(x){ ####These are to plot and warn in the event of identical sequences in the reference data#####
      x<-x+runif(1,0,1)/70  ############################################################################################# 
      }                     #############################################################################################  
#for(j in 1:nrow(testdata)){
for(j in 1:testk){
if(identical.indicator[j]==TRUE){
                  posteriors2<- t(apply(posteriors[[j]]$post,1,jitter2))
plot(c(1:ncol(posteriors[[j]]$post)),posteriors2[1,],type="l",col=1,ylim=range(0,1),ylab="Posterior Probabilities",xlab="Position",cex.main=1,font.main=4,main=c(paste(paste("Sequence",j),paste("Species:","Multiple"),paste("Probability:",species.class[,j][[2]]),sep="\n")))
      mtext("Note: If the plot above is noisy and uses all of the available positions,",side=1,line=4,cex=.8,font=3)
      mtext("this sequence may belong to a species that is not in the reference data set.",side=1,line=4.75,cex=.8,font=3)
      legend("topleft",c(rownames(testpriors)[which(posteriors[[j]]$post[,ncol(posteriors[[j]]$post)]==max(posteriors[[j]]$post[,ncol(posteriors[[j]]$post)]))]),title.col="red3",title="Warning: Multiple Species Identified",xjust=1,cex=.8,inset=.02,bty="n")

for(i in 2:nrow(posteriors[[j]]$post)){
lines(c(1:ncol(posteriors[[j]]$post)),posteriors2[i,],type="l",col=i)
}

      ifelse(plot.file=="pdf",dev.copy(pdf,paste("seq",j,".pdf")),ifelse(plot.file=="jpg",dev.copy(jpeg,paste("seq",j,".jpg")),ifelse(plot.file=="png",dev.copy(png,paste("seq",j,".png")),ifelse(plot.file=="wmf",dev.copy(win.metafile,paste("seq",j,".wmf")),ifelse(plot.file=="ps",dev.copy(postscript,paste("seq",j,".ps")),dev.copy(pdf,paste("seq",j,".pdf")))))))
      dev.off()
}

else {
plot(c(1:ncol(posteriors[[j]]$post)),posteriors[[j]]$post[1,],type="l",col=1,ylim=range(0,1),ylab="Posterior Probabilities",xlab="Position",cex.main=.8,font.main=4,main=c(paste(paste("Sequence",j),paste("Species:",species.class[,j][[1]]),paste("Probability:",species.class[,j][[2]]),sep="\n")))
      mtext("Note: If the plot above is noisy and uses all of the available positions,",side=1,line=4,cex=.8,font=3)
      mtext("this sequence may belong to a species that is not in the reference data set.",side=1,line=4.75,cex=.8,font=3)

for(i in 2:nrow(posteriors[[j]]$post)){
lines(c(1:ncol(posteriors[[j]]$post)),posteriors[[j]]$post[i,],type="l",col=i)
}

      ifelse(plot.file=="pdf",dev.copy(pdf,paste("seq",j,".pdf")),ifelse(plot.file=="jpg",dev.copy(jpeg,paste("seq",j,".jpg")),ifelse(plot.file=="png",dev.copy(png,paste("seq",j,".png")),ifelse(plot.file=="wmf",dev.copy(win.metafile,paste("seq",j,".wmf")),ifelse(plot.file=="ps",dev.copy(postscript,paste("seq",j,".ps")),dev.copy(pdf,paste("seq",j,".pdf")))))))
      dev.off()
}

}






    result <- list(k = k, totaltime = tottime,delta = delta, species.class = species.class, 
        priors = priors,posteriors=posteriors)

}

