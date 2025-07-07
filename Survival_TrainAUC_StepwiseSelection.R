library(nsROC)
library(caret)
library(survival)

forward_step <- function(dat, candid, fixvar, std.time, numSeed){
  forward_ls <- NULL
  for (s in seq(numSeed)){
    set.seed(s)
    trIdx<-createDataPartition(dat$OS,p=.7,list = FALSE,times = 1)
    trdat<-dat[trIdx,]
    tsdat<-dat[-trIdx,]
    for (g in setdiff(candid,fixvar)){
      f=as.formula(paste('Surv(OS.year,OS) ~ ',paste(fixvar,collapse = ' + '),' + ',g,collapse = ''))
      trdat1 <- trdat[complete.cases(trdat[,c('OS.year','OS',setdiff(c(fixvar,g),''))]),c('OS.year','OS',setdiff(c(fixvar,g),''))]
      tsdat1 <- tsdat[complete.cases(tsdat[,c('OS.year','OS',setdiff(c(fixvar,g),''))]),c('OS.year','OS',setdiff(c(fixvar,g),''))]
      CoxPHres<-coxph(f,data = trdat1)
      lptr <- predict(CoxPHres,trdat1)
      if (max(trdat1$OS.year)>=std.time){
        trauc<-cdROC(stime=trdat1$OS.year,status=trdat1$OS,marker = lptr,predict.time = std.time)$auc
      } else{
        trauc <- NA
      }
      lpts <- predict(CoxPHres,tsdat1)
      if (max(tsdat1$OS.year)>=std.time){
        tsauc<-cdROC(stime=tsdat1$OS.year,status=tsdat1$OS,marker = lpts,predict.time = std.time)$auc
      } else{
        tsauc <- NA
      }
      forward_ls <- rbind(forward_ls,c(s,paste(c(fixvar,g),collapse = ' + '),trauc,tsauc))
    }
  }
  forward_ls1 <- data.frame(forward_ls)
  forward_ls1$X3 <- as.numeric(forward_ls1$X3)
  forward_ls1$X4 <- as.numeric(forward_ls1$X4)
  AUCsumm <- NULL
  for (g in unique(forward_ls1$X2)){
    AUCsumm <- rbind(AUCsumm,c(g,colMeans(forward_ls1[which(forward_ls1$X2==g),-c(1,2)],na.rm = T)))
  }
  return(AUCsumm)
}

backward_step <- function(dat, backcandid, fixvar, std.time, numSeed){
  backward_ls <- NULL
  for (s in seq(numSeed)){
    set.seed(s)
    trIdx<-createDataPartition(dat$OS,p=.7,list = FALSE,times = 1)
    trdat<-dat[trIdx,]
    tsdat<-dat[-trIdx,]
    for (g in backcandid){
      f=as.formula(paste('Surv(OS.year,OS) ~ ',paste(setdiff(fixvar,g),collapse = ' + '),collapse = ''))
      trdat1 <- trdat[complete.cases(trdat[,c('OS.year','OS',setdiff(fixvar,g))]),]
      tsdat1 <- tsdat[complete.cases(tsdat[,c('OS.year','OS',setdiff(fixvar,g))]),]
      CoxPHres<-coxph(f,data = trdat1)
      lptr <- predict(CoxPHres,trdat1)
      if (max(trdat1$OS.year)>=std.time){
        trauc<-cdROC(stime=trdat1$OS.year,status=trdat1$OS,marker = lptr,predict.time = std.time)$auc
      } else{
        trauc <- NA
      }
      lpts <- predict(CoxPHres,tsdat1)
      if (max(tsdat1$OS.year)>=std.time){
        tsauc<-cdROC(stime=tsdat1$OS.year,status=tsdat1$OS,marker = lpts,predict.time = std.time)$auc
      } else{
        tsauc <- NA
      }
      backward_ls <- rbind(backward_ls,c(s,paste(setdiff(fixvar,g),collapse = ' + '),trauc,tsauc))
    }
  }
  backward_ls1 <- data.frame(backward_ls)
  backward_ls1$X3 <- as.numeric(backward_ls1$X3)
  backward_ls1$X4 <- as.numeric(backward_ls1$X4)
  AUCsumm <- NULL
  for (g in unique(backward_ls1$X2)){
    AUCsumm <- rbind(AUCsumm,c(g,colMeans(backward_ls1[which(backward_ls1$X2==g),-c(1,2)],na.rm = T)))
  }
  return(AUCsumm)
}

TrainAUCStepwise <- function(totvar,dat,fixvar,excvar,std.time,numSeed,imtres,addres1){
  while (length(setdiff(fixvar,excvar))<length(totvar)){
    candid <- setdiff(totvar,c(fixvar,excvar))
    forward_ls <- forward_step(dat, candid, fixvar, std.time, numSeed)
    forward.trauc1<-max(as.numeric(forward_ls[,2]))
    forward.var1 <- forward_ls[which.max(as.numeric(forward_ls[,2])),1];forward.tsauc1 <-forward_ls[which.max(as.numeric(forward_ls[,2])),3]
    forward.newstep <-matrix(c(forward.var1,forward.trauc1,forward.tsauc1),nrow=1)
    if (length(fixvar)==1){
      dir.create(addres1)
      imtres <- rbind(imtres, forward.newstep)
      colnames(forward.newstep)<-c('Variable','trainAUC','testAUC')
      eval(parse(text = paste("write.csv(forward.newstep,'./",addres1,"/Intermediate_Forward",nrow(imtres),".csv',row.names = F)",sep = "")))
      fixvar <- append(fixvar,forward.var1)
      forward.old <- forward.newstep
    } else{
      forward.old <- imtres[nrow(imtres),]
      if (as.numeric(forward.newstep[2]) > as.numeric(forward.old[2])){
        colnames(forward.newstep) <- colnames(imtres)
        imtres <- rbind(imtres, forward.newstep)
        colnames(forward.newstep)<-c('Variable','trainAUC','testAUC')
        eval(parse(text = paste("write.csv(forward.newstep,'./",addres1,"/Intermediate_Forward",nrow(imtres),".csv',row.names = F)",sep = "")))
        fixvar <- append(fixvar,gsub(' ','',strsplit(forward.var1,'\\+')[[1]])[length(gsub(' ','',strsplit(forward.var1,'\\+')[[1]]))])
        forward.old <- forward.newstep
        if (nrow(imtres)>2){
          backcandid<-fixvar[c(1:(length(fixvar)-2))]
          backward_ls<-backward_step(dat, backcandid, fixvar, std.time, numSeed)
          backward.trauc1<-max(as.numeric(backward_ls[,2]))
          backward.var1 <- backward_ls[which.max(as.numeric(backward_ls[,2])),1];backward.tsauc1 <-backward_ls[which.max(as.numeric(backward_ls[,2])),3]
          backward.newstep <-matrix(c(backward.var1,backward.trauc1,backward.tsauc1),nrow=1)
          if (backward.trauc1 > as.numeric(forward.old[2])){
            imtres <- rbind(imtres, backward.newstep)
            colnames(backward.newstep)<-c('Variable','trainAUC','testAUC')
            eval(parse(text = paste("write.csv(backward.newstep,'./",addres1,"/Intermediate_Backward",nrow(imtres),".csv',row.names = F)",sep = "")))
            fixvar <- gsub(' ','',strsplit(backward_ls[which.max(as.numeric(backward_ls[,2])),1],'\\+')[[1]])
            forward.old <- backward.trauc1
          }
        }
      } else{
        mat<-matrix(imtres[nrow(imtres),],nrow=1)
        colnames(mat)<-c('Variable','trainAUC','testAUC')
        colnames(imtres)<-c('Variable','trainAUC','testAUC')
        eval(parse(text = paste("write.csv(imtres,'./",addres1,"/Intermediate_Stepwise_Total.csv',row.names = F)",sep = "")))
        eval(parse(text = paste("write.csv(mat,'./",addres1,"/Final_Stepwise_Total.csv',row.names = F)",sep = "")))
        break
      }
    }
  }
  mat<-matrix(imtres[nrow(imtres),],nrow=1)
  colnames(mat)<-c('Variable','trainAUC','testAUC')
  colnames(imtres)<-c('Variable','trainAUC','testAUC')
  eval(parse(text = paste("write.csv(imtres,'./",addres1,"/Intermediate_Stepwise_Total.csv',row.names = F)",sep = "")))
  eval(parse(text = paste("write.csv(mat,'./",addres1,"/Final_Stepwise_Total.csv',row.names = F)",sep = "")))
}
