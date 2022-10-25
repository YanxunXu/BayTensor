########## MAP Estimators ###########
MAP<-function(Y,X,core_dim,seed, total_iter = 100){
  ## core_dim = the chosen core tensor dimensions
  iter<-0;diff<-10;bbbb<-0 ## Initial the criteria
  tol<-1e-05;## Tolerance
  set.seed(seed)
  dimX<-dim(X)[-1];dimY<-dim(Y)[-1];N<-dim(X)[1]
  nx<-length(dimX);ny<-length(dimY)
  sigma<-1 ## Sigma^2 for error term
  
  ## Initialize core tensor
  G_est<-array(rnorm(prod(core_dim)),core_dim)
  
  ## Initialize P's and Q's
  Xnames<-paste("P",1:nx,sep="")
  for(i in 1:nx){
    assign(Xnames[i],matrix(runif(dimX[i]*core_dim[i]),dimX[i],core_dim[i]))
  }
  Ynames<-paste("Q",1:ny,sep="")
  Ynames_matrix<-paste("Y_matrix",1:ny,sep="")
  for(i in 1:ny){
    assign(Ynames[i],matrix(runif(dimY[i]*core_dim[i+nx]),dimY[i],core_dim[i+nx]))
    assign(Ynames_matrix[i],array_to_matrix(Y,(1:(ny+1))[-(i+1)],(i+1)))
  }
  Y_1<-array_to_matrix(Y,1,c(2,(ny+1)))
  
  ## While condition not satisfied
  while(diff > tol && iter<total_iter){
    ## Update P's
    for(i in 1:nx){
      temp<-G_est
      p_left<-(1:nx)[-i];x_left<-(2:(nx+1))[-i]
      for(j in p_left){
        temp<-tensor_matrix(temp,get(Xnames[j]),j)
      }
      for(k in 1:ny){
        temp<-tensor_matrix(temp,get(Ynames[k]),(k+nx))
      }
      Xtemp<-tensor_tensor(temp,X,p_left,x_left)
      ## Order of dimensions
      ord<-c((2+ny),(2:(ny+1)),(ny+3),1)
      ## Reorder Xtemp
      Xtemp<-aperm(Xtemp,ord)
      ## Convert to matrix
      Xtemp<-array_to_matrix(Xtemp,c(1:(ny+1)),c(ny+2,ny+3))
      DX<-dim(Xtemp)[2]
      ## New mu, and Sigma matrix
      mu<-solve(((t(Xtemp)%*%Xtemp/sigma)+diag(1,DX)),((t(Xtemp))%*%as.vector(Y)/sigma),tol=1e-25)
      assign(Xnames[i],matrix(mu,dimX[i],core_dim[i]))
    }
    
    ## Update Q's
    for(i in 1:ny){
      temp<-G_est
      q_left<-(1:ny)[-i];y_left<-(2:(ny+1))[-i]
      for(j in 1:nx){
        temp<-tensor_matrix(temp,get(Xnames[j]),j)
      }
      for(k in q_left){
        temp<-tensor_matrix(temp,get(Ynames[k]),(k+nx))
      }
      Xtemp<-tensor_tensor(temp,X,c(1:nx),c(2:(nx+1)))
      ord<-c((ny+1),q_left,i)
      Xtemp<-aperm(Xtemp,ord)
      Xtemp<-array_to_matrix(Xtemp,c(1:ny),(ny+1))
      DX<-dim(Xtemp)[2]
      ## New mu, and Sigma matrix
      mu<-solve(((t(Xtemp)%*%Xtemp/sigma)+diag(1,DX)),t(Xtemp)%*%get(Ynames_matrix[i])/sigma,tol=1e-25)
      mu <- mu + 1e-04 * matrix(runif(prod(dim(mu))), dim(mu)[1], dim(mu)[2]) # Small perturbation
      assign(Ynames[i],t(mu))
    }
    
    ## Update G
    ## Matrix outer product of P's
    Ax<-get(Xnames[nx])
    for(i in (nx-1):1){
      Ax<-matrix_outer(Ax,get(Xnames[i]))
    }
    ## Matrix outer product of Q's
    By<-get(Ynames[ny])
    for(i in (ny-1):1){
      By<-matrix_outer(By,get(Ynames[i]))
    }
    XX<-array_to_matrix(X,1,c(2:(nx+1)))
    DBy<-dim(By)[2]
    DXAx<-dim(XX%*%Ax)[2]
    ## X(1) multiply by U and make it block matrix
    NN<-matrix_outer(diag(DBy),XX%*%Ax)
    ## Matrix V(V^TV)^(-1)
    MM<-By%*%solve((t(By)%*%By),tol=1e-25)
    ## Vectorize of Y tilde
    Yv<-as.vector(Y_1%*%MM)
    
    Sig_temp<-matrix(0,N*DBy,N*DBy)
    for(i in 1:DBy){
      for(j in 1:DBy){
        row_s<-c(((i-1)*N+1):(i*N))
        col_s<-c(((j-1)*N+1):(j*N))
        tem<-t(MM[,i])%*%MM[,j]*sigma
        Sig_temp[row_s,col_s]<-diag(as.numeric(tem),N)
      }
    }
    
    Sig<-solve((t(NN)%*%solve(Sig_temp,tol=1e-25)%*%NN+diag(dim(NN)[2])),tol=1e-25)
    mu<-Sig%*%t(NN)%*%solve(Sig_temp,tol=1e-25)%*%Yv
    
    G_est<-array(mu,core_dim)
    
    
    #####################################################################################################
    
    ## Compare criteria
    iter<-iter+1
    B<-G_est
    for(i in 1:nx){
      B<-tensor_matrix(B,get(Xnames[i]),i)
    }
    for(i in (nx+1):(nx+ny)){
      B<-tensor_matrix(B,get(Ynames[i-nx]),i)
    }
    Y_est<-tensor_tensor(X,B,c(2:(nx+1)),c(1:nx))
    errY<-as.vector(Y_est)-as.vector(Y)
    
    ## Update sigma
    sigma<-1/rgamma(1,shape=N*(dim(By)[1])/2+1,scale=1/(1+sum(errY^2)/2))
    
    temp<-bbbb
    bbbb<-sum(errY^2)/sum(as.vector(Y)^2)
    diff<-abs(temp-bbbb)
  }
  
  
  
  BIC_N<-length(errY)
  par<-prod(core_dim)+sum(core_dim*c(dimX,dimY))
  BIC<-BIC_N*log(sigma)+sum(errY^2)/sigma+log(BIC_N)*par
  
  result<-list(B_est=B,Y_est=Y_est,BIC=BIC)
  return(result)
}

####### BayTensor Fast ##################
BayTensor_Fast<-function(Y, X, K0, seed, K = NULL){
  if(is.null(K)){
    K = K0 + 2
  }
  core_dim <- K0
  dim_record<-core_dim
  dim_total = length(core_dim)
  SA_p<-0
  SA_t<-100
  SA_gamma<-0.9
  sim<-MAP(Y, X, core_dim, seed, 20)
  dim_seen<-matrix(c(core_dim,sim$BIC),1)
  bic<-sim$BIC 
  SA_i <- 1
  count<-10
  while(SA_p<5 && count>1){
    SA_tp<-SA_gamma^SA_p*SA_t
    count<-0
    for(SA_k in c(1:10)){
      ## Each temperature for 10 iterations
      candidate<-Find_candidate(core_dim,K)
      set.seed(floor(1234*pi*(SA_k*(SA_p+1))))
      SA_choice<-sample(1:dim(candidate)[1],1)
      dim_temp<-core_dim+candidate[SA_choice,]
      
      cat('Cur dim:',core_dim,'\n')
      flag=0
      for(i in 1:dim(dim_seen)[1]){
        if(prod(dim_temp==dim_seen[i,1:dim_total])){
          BIC_temp<-dim_seen[i,dim_total+1]
          flag=1
          break
        }
      }
      if(flag==0){
        sim<-MAP(data$Y, data$X, dim_temp, seed, 20)
        BIC_temp<-sim$BIC
        dim_seen<-rbind(dim_seen,c(dim_temp,BIC_temp))
      }
      
      cat('BIC0:', bic[SA_i], 'BIC1:', BIC_temp, '\n')
      cat('Candidate Dim:', dim_temp,'\n')
      
      SA_u<-runif(1)
      SA_i<-SA_i+1
      if(BIC_temp<bic[SA_i-1] | SA_u<exp((bic[SA_i-1]-BIC_temp)/SA_tp)){
        bic<-c(bic, BIC_temp)
        core_dim<-dim_temp
        count<-count+1
      }else{
        bic<-c(bic, bic[SA_i-1])
      }
      dim_record<-rbind(dim_record,core_dim)
    }
    SA_p<-SA_p+1
  }
  cat('Core Dim:', core_dim, '\n')
  sim<-MAP(data$Y, data$X, core_dim, seed, 1000)
  return(sim)
}
