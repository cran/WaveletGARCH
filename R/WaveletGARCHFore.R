
##############################################
#         Wavelet GARCH Forecast Function         #
##############################################


WaveletGARCHFore<-function(series,filtern,level,nofore)
{
  newenv<-new.env()

  MaxARParam      <- 5
  MaxMAParam      <- 5
  MaxAR_ARCH      <- 4
  MaxMA_ARCH      <- 4
  conditionDist   <- "norm"
  flag            <- 1
  ClusterExportData <- c("tsrss")
  series<-as.vector(series)
  tsrss<-NULL
  dflag<-0
  mod<-rep(c("no"),times=level+1)
  modell<-{}

  ##############################################
  #         Fitting Auto GARCH                 #
  ##############################################


  garchAutoTryFit = function(ll)
  {

    formula <- paste0("ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(",ll$order[3],",", ll$order[4],")),mean.model = list(armaOrder = c(",ll$order[1],",", ll$order[2],"), include.mean = TRUE),distribution.model = '",ll$dist,"')")
    fit = tryCatch(ugarchfit(spec=eval(parse(text=formula)),data=tsrss),error=function(err) TRUE,warning=function(warn) TRUE )
    if(is.logical(fit)){
      fit = NULL
    }
    return(fit)
  }



  ##############################################
  #         Auto GARCH Function                #
  ##############################################


  garchAuto = function(tsrss,min.order,max.order,cond.dists,arma.sum,ic,cores)
  {

    len = NROW(tsrss)

    models = list( )
    for( dist in cond.dists )
      for( p in min.order[1]:max.order[1] )
        for( q in min.order[2]:max.order[2] )
          for( r in min.order[3]:max.order[3] )
            for( s in min.order[4]:max.order[4] )
            {
              pq.sum = p + q
              if( pq.sum <= arma.sum[2] && pq.sum >= arma.sum[1] )
              {
                models[[length(models) + 1]] = list(order=c( p, q, r, s ),dist=dist)
              }
            }

    cl <- makeCluster(cores)

    c1<-clusterEvalQ(cl=cl, library(rugarch))

    c2<-clusterExport(cl, ClusterExportData, envir = environment())

    res = parSapply(cl, X=c(models),FUN=garchAutoTryFit)

    stopCluster(cl)

    GoodModel <- res[!sapply(res, is.null)]

    AICWithParam <- NULL
    if(length(GoodModel)==0)
    {

      dflag=-1
      return(NULL)
    }
    else
    {

    for(i in 1:length(GoodModel)){
      temp <- data.frame(infocriteria(GoodModel[[i]])[1,],GoodModel[[i]]@fit$LLH,length(GoodModel[[i]]@fit$coef),i,do.call(paste, c(as.list(names(GoodModel[[i]]@fit$coef)), sep="-")))
      AICWithParam <- rbind(AICWithParam,temp)
    }
    colnames(AICWithParam) <- c("AIC","Likelihood",'NoofParams','SeqenceNumb','Model-Specificaion')

    AICWithParamSorting <- AICWithParam[order(AICWithParam[,1],AICWithParam[,2],decreasing=FALSE),]

    best.fit <- GoodModel[[AICWithParamSorting[1,4]]]
    best.ic <- AICWithParamSorting[1,1]
    if( best.ic < arma.sum[2] )
    {
      return( best.fit )
    }
    return( NULL )
  }

  }
  dpoints<-length(series)
  w<-matrix(ncol=level, nrow= dpoints)
  v<-matrix(ncol=1,nrow = dpoints)
  nfor<-nofore


  dwt<-modwt(X=series, filter=filtern, n.levels=level)

  v[,1]<-dwt@V[[level]]

  for(i in (1:(level)))
  {

    w[,i]<-dwt@W[[i]]

  }

  total<-cbind(w,v)

  deco<-{}
  for(i in 1:(level+1))

  {

    ARCHLMOrigSeries <- as.numeric(ArchTest(total[,i])$p.value)


    assign("tsrss",as.list(total[,i]),envir = environment())

    ARCHLMCutoff<-0.05
    rrow<-NROW(tsrss)

    if(ARCHLMOrigSeries < ARCHLMCutoff)

    {


      garch.auto <-  garchAuto(tsrss,min.order=c(0,0,1,0),
                               max.order=c(MaxARParam=2,MaxMAParam=3,MaxAR_ARCH=2,MaxMA_ARCH=1),
                               cond.dists="snorm",#"sged", "snorm", "ged", "std", "sstd","snig","QMLE"
                               arma.sum=c(0,1e9),
                               cores=4,
                               ic="BIC")

      if(dflag==-1)
      {
        arima.auto<-auto.arima(total[,i])

        deco[flag]<-as.data.frame((forecast(arima.auto, h=nfor,level=level))$mean)
        coeff<-arima.auto$coef
        sig<-arima.auto$sigma2
        ai<-arima.auto$aic
        aii<-arima.auto$aicc
        bi<-arima.auto$bic
        orders<-arimaorder(arima.auto)
        arima.obj<-new("autoarima",coefficient=coeff,sigma2=sig,aic=ai,aicc=aii,bic=bi,order=orders)
        mod[flag]="Auto Arima"

        modell<-c(modell,arima.obj)
        flag=flag+1
        dflag=0

      }
      else
      {
        garchfore<-ugarchforecast(garch.auto, n.ahead=nfor)
        deco[flag]<-as.data.frame(garchfore@forecast$seriesFor)
        mod[flag]="auto garch"

        modell<-c(modell,garch.auto)
        flag=flag+1

      }





    }

    else

    {

      arima.auto<-auto.arima(total[,i])

      deco[flag]<-as.data.frame((forecast(arima.auto, h=nfor,level=level))$mean)
      coeff<-arima.auto$coef
      sig<-arima.auto$sigma2
      ai<-arima.auto$aic
      aii<-arima.auto$aicc
      bi<-arima.auto$bic
      orders<-arimaorder(arima.auto)
      arima.obj<-new("autoarima",coefficient=coeff,sigma2=sig,aic=ai,aicc=aii,bic=bi,order=orders)
      mod[flag]="Auto Arima"

      modell<-c(modell,arima.obj)
      flag=flag+1

    }


  }

  for(i in 1:level)
  {

    dwt@W[[i]]<-as.matrix(unlist(deco[i]))

  }

  dwt@V[[level]]<-as.matrix(unlist(deco[level+1]))

  forecastobject<-imodwt(dwt)

  forecastobject<-as.matrix(forecastobject)

  colnames(forecastobject)<-"Forecast Values"

  rownames(forecastobject)<-c(1:length(forecastobject))

  reslt<-list(r1=forecastobject,r2=mod,r3=modell,r4=level+1)
  class(reslt)<-"WaveletGARCHFore"

  return(reslt)

}



print.WaveletGARCHFore<-function(x,...)
{

  cat("                          **************************************\n")
  cat("                          *         Wavelet GARCH Results      *\n")
  cat("                          **************************************\n")
  cat("\n")

  cr1<-x$r1
  cr2<-x$r2
  cr3<-x$r3
  cr4<-x$r4

  cat("\n")

  cat("         Models used for differnt wavelet coefficients are \n\n")

  for (i in 1:cr4)
  {

    if(i==cr4)
    {
      cat("         For V",i-1,"model selected is: \n")
      cat("         -----------------------------------------------\n")
      cat("                          ",cr2[i])
      cat("\n")
      cat("         -----------------------------------------------\n")
      cat("\n")
      print(cr3[i])

    }

    else
    {
      cat("         For W",i,"model selected is: \n")
      cat("         -----------------------------------------------\n")
      cat("                          ",cr2[i])
      cat("\n")
      cat("         -----------------------------------------------\n")
      cat("\n")
      print(cr3[i])

    }
  }

  cat("         -----------------------------------------------\n")
  cat("                       List of Forecast Values")
  cat("\n")
  cat("         -----------------------------------------------\n")
  cat("\n")

  print(cr1)

}


