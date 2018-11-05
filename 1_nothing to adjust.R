k=k+1
var.rm <- names(sort(data.imp))[k]
ratio.num.den <- c(as.character(ratios$v1[which(var.rm == ratios$dfnames)[[1]]]),as.character(ratios$v2[which(var.rm == ratios$dfnames)[[1]]]))

if((ratio.num.den[1] %in% proteins.adj || ratio.num.den[2] %in% proteins.adj) == FALSE){
  output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=0,"var.rm.num"=ratio.num.den[1],"var.rm.den"=ratio.num.den[2], "corr.train.cor"=0, "p.train"=0, "corr.valid.cor"=0, "p.valid"=0, "RMSE.train"=0, "RMSE.valid"=0))
}
