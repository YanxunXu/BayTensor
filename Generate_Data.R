#### Generate Simulation Data: Section 6.2######
Generate_data<-function(N,dimX,dimY,seed,B, t_and_s, SNR){
  set.seed(seed)
  nx<-length(dimX);ny<-length(dimY)
  
  ## Generate X
  X<- array(0.01 * rnorm(N*prod(dimX)),c(N,dimX)) # Change this such that sigma of error is 1
  E<-array(rnorm(N*prod(dimY)),c(N,dimY))
  
  if(!(class(B)=='array')){
    ## Assume core tensor dims (t1,t2,s1,s2) uniformly from (0,1)
    t<-rep(t_and_s[1],nx);s<-rep(t_and_s[2],ny)
    #t<-c(2,2);s<-c(2,2)
    B<-array(runif(prod(t)*prod(s)),c(t,s))
    #SNR<-5
    ## Factor matrices according to X dims
    Xnames<-paste("P",1:nx,sep="")
    for(i in 1:nx){
      assign(Xnames[i],matrix(rnorm(dimX[i]*t[i]),dimX[i],t[i]))
    }
    
    ## Factor matrices according to Y dims
    Ynames<-paste("Q",1:ny,sep="")
    for(i in 1:ny){
      assign(Ynames[i],matrix(rnorm(dimY[i]*s[i]),dimY[i],s[i]))
    }
    
    ## Calculate the true coefficient tensor B
    for(i in 1:nx){
      B<-tensor_matrix(B,get(Xnames[i]),i)
    }
    for(i in (nx+1):(nx+ny)){
      B<-tensor_matrix(B,get(Ynames[i-nx]),i)
    }
  }
  ## Calculate the true Y tensor
  Y_TEMP<-tensor_tensor(X,B,c(2:(nx+1)),c(1:nx))
  Y<-Y_TEMP+sqrt((sum(Y_TEMP^2)/sum(E^2))/SNR)*E
  ## Return a list of X, Y, and B
  data<-list(X=X,Y=Y,B=B)
  return(data)
}

########## Generate Correlated Simulation Data: Section 6.3 ###########
Generate_data_cor<-function(N,dimX,dimY,seed, B, t_and_s, SNR){
  set.seed(seed)
  nx<-length(dimX);ny<-length(dimY)
  XX <- rep(1:dimX[2], each=dimX[1])
  YY <- rep(1:dimX[1], dimX[2])
  coords <- data.frame(XX,YY)
  # distance matrix
  dist <- as.matrix(dist(coords))
  ## Correlation function
  str <- -0.1 # strength of autocorrelation, inv. proportional to str
  omega1 <- exp(str*dist)
  
  # calculate correlation weights, and invert weights matrix
  weights <- chol(solve(omega1))
  weights_inv <- solve(weights)
  
  # create an autocorrelated random field
  X<-array(0,c(N,dimX))
  for(i in 1:N){
    Xtemp<-weights_inv %*% (0.05 * rnorm(prod(dimX)))
    X[i,,]<-matrix(Xtemp,nrow=dimX[1])
  }
  #SNR<-50
  ## Generate X
  E<-array(rnorm(N*prod(dimY)),c(N,dimY))
  if(!(class(B)=='array')){
    ## Assume core tensor dims (t1,t2,s1,s2) uniformly from (0,1)
    #t<-rep(2,nx);s<-rep(2,ny)
    #t<-c(2,2);s<-c(2,2)
    t<-rep(t_and_s[1],nx);s<-rep(t_and_s[2],ny)
    B<-array(runif(prod(t)*prod(s)),c(t,s))
    
    ## Factor matrices according to X dims
    Xnames<-paste("P",1:nx,sep="")
    for(i in 1:nx){
      assign(Xnames[i],matrix(rnorm(dimX[i]*t[i]),dimX[i],t[i]))
    }
    
    ## Factor matrices according to Y dims
    Ynames<-paste("Q",1:ny,sep="")
    for(i in 1:ny){
      assign(Ynames[i],matrix(rnorm(dimY[i]*s[i]),dimY[i],s[i]))
    }
    
    ## Calculate the true coefficient tensor B
    for(i in 1:nx){
      B<-tensor_matrix(B,get(Xnames[i]),i)
    }
    for(i in (nx+1):(nx+ny)){
      B<-tensor_matrix(B,get(Ynames[i-nx]),i)
    }
  }
  ## Calculate the true Y tensor
  Y_TEMP<-tensor_tensor(X,B,c(2:(nx+1)),c(1:nx))
  Y<-Y_TEMP+sqrt((sum(Y_TEMP^2)/sum(E^2))/SNR)*E
  #print(sqrt((sum(Y_TEMP^2)/sum(E^2))/SNR))
  ## Return a list of X, Y, and B
  data<-list(X=X,Y=Y,B=B)
  return(data)
}
