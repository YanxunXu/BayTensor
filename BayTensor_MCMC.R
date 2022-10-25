######### Bayesian Posterior ########
Estimate<-function(Y, X, core_dim, seed, b, total_iter = 200, collect = FALSE){
  ## core_dim = the chosen core tensor dimensions
  start_time = proc.time()
  iter<-0;#diff<-10;bbbb<-0 ## Initial the criteria
  n_burnin <- 100
  
  if(collect == TRUE){
    i_record <- 1
    result_post <- list()
    result_post[['B_estimate']] = list()
    result_post[['Y_estimate']] = list()
  }
  
  #tol<-1e-05;total_iter<-500 ## Tolerance
  set.seed(seed)
  dimX<-dim(X)[-1];dimY<-dim(Y)[-1];N<-dim(X)[1]
  nx<-length(dimX);ny<-length(dimY)
  
  gamma_alpha <- 1
  gamma_beta <- 1
  
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
  #while(diff > tol && iter<total_iter){
  while(iter < total_iter){
    ## Update P's
    for(i in 1:nx){
      #cat('Update P', i, '\n')
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
      
      mu<-solve(((t(Xtemp)%*%Xtemp*b/sigma)+diag(1,DX)),((t(Xtemp))%*%as.vector(Y)*b/sigma),tol=1e-25)
      Sig<-solve((t(Xtemp)%*%Xtemp*b/sigma)+diag(1,DX),tol=1e-25) 
      Sig <- Sig + diag(1e-08, dim(Sig)[1])
      
      ## Generate from posterior distribution
      Ptemp<-rmvn(1,mu,Sig)
      assign(Xnames[i], matrix(Ptemp, dimX[i], core_dim[i]))
    }
    
    ## Update Q's
    for(i in 1:ny){
      #cat('Update Q', i, '\n')
      
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
      
      mu<-solve(((t(Xtemp)%*%Xtemp*b/sigma)+diag(1,DX)),t(Xtemp)%*%get(Ynames_matrix[i])*b/sigma,tol=1e-25)
      Sig<-solve((t(Xtemp)%*%Xtemp*b/sigma)+diag(1,DX),tol=1e-25)
      Sig <- Sig + diag(1e-08, dim(Sig)[1])
      
      Vtemp<-mu
      ### Update each column of V^T
      for(ttt in 1:dim(mu)[2]){
        Vtemp[,ttt]<-rmvn(1,mu[,ttt],Sig)
      }
      Vtemp <- Vtemp + 1e-04 * matrix(runif(prod(dim(Vtemp))), dim(Vtemp)[1], dim(Vtemp)[2]) # Small perturbation
      assign(Ynames[i], t(Vtemp))
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
    Sig_temp <- Sig_temp + diag(1e-06, dim(Sig_temp)[1])
    Sig<-solve((b*t(NN)%*%solve(Sig_temp,tol=1e-25)%*%NN+diag(dim(NN)[2])),tol=1e-25)
    Sig<- Sig + diag(1e-06, dim(Sig)[1])
    
    mu<-b*Sig%*%t(NN)%*%solve(Sig_temp,tol=1e-25)%*%Yv
    GXY<-mvrnorm(1,mu,Sig)
    G_est<-array(GXY,core_dim)
    
    #####################################################################################################
    
    ## Compare criteria
    B<-G_est
    for(i in 1:nx){
      B<-tensor_matrix(B,get(Xnames[i]),i)
    }
    for(i in (nx+1):(nx+ny)){
      B<-tensor_matrix(B,get(Ynames[i-nx]),i)
    }
    Y_est<-tensor_tensor(X,B,c(2:(nx+1)),c(1:nx)) # Y Predictive Mean 
    errY<-as.vector(Y_est)-as.vector(Y)
    
    ## Update sigma
    sigma<-1/rgamma(1,shape=b*N*(dim(By)[1])/2+ gamma_alpha, scale=1/(gamma_beta+b*sum(errY^2)/2))
    
    if((collect == TRUE) & (iter >= n_burnin)){
      result_post[['B_estimate']][[i_record]] = B
      result_post[['Y_estimate']][[i_record]] = Y_est
      result_post[['Sigma2_estimate']][[i_record]] = sigma
      i_record = i_record + 1
    }
    
    iter<-iter+1
    if(iter %% 500 == 0){
      run_time = proc.time() - start_time
      cat('Iteration #', iter , '\n')
      #cat('Running Time is ', round(run_time[3] / 3600,2), 'hours', '\n')
    }
  }
  
  if(collect){
    return(result_post)
  }else{
    result<-list(B_est=B, Y_est = Y_est, Sigma=sigma)
    return(result)
  }
}

####### Calculate Log likelihood #########
Likelihood<-function(sim,Y){
  temp<-sum((as.vector(sim$Y_est)-as.vector(Y))^2)
  n<-length(as.vector(Y))
  L<- - n/2*log(sim$Sigma)-temp/sim$Sigma/2
  return(L)
}

####### BayTensor MCMC ##################
BayTensor_MCMC<-function(Y, X, K0, seed, b, n_iter, K = NULL, lambda = 1, I = 100){
  # K0 : Initial Guess
  # K : Upper bound for dimension
  # lambda : Poisson Prior with Lambda
  # I : Number of iteration to select core tensor dimension
  
  start_time = proc.time()
  if(is.null(K)){
    K = K0 + 2
  }
  set.seed(seed)
  seed_ls = sample(c(1:10000), I + 50, replace = FALSE)
  
  core_dim = K0
  dim_total = length(core_dim)
  seed_star = seed_ls[1]
  sim <- Estimate(Y, X, core_dim, seed_star, b, total_iter = 100, collect = FALSE)
  L0 <- Likelihood(sim, Y)
  
  dim_seen<-matrix(c(core_dim,L0),1)
  
  for(i in c(2:I)){
    cat('Iteration #', i, '\n')
    candidate<-Find_candidate(core_dim,K)
    choice<-sample(1:dim(candidate)[1],1)
    dim_temp<-core_dim+candidate[choice,]
    q_01<-1/dim(candidate)[1]
    q_10<-1/dim(Find_candidate(dim_temp,K))[1]
    
    ### Poission Prior with lambda
    ss<-which(candidate[choice,]!=0)
    if(sum(candidate[choice,])<0){
      ratio<-core_dim[ss]/lambda
    }else{
      ratio<-lambda/(core_dim[ss]+1)
    }
    
    flag=0
    ## See if dimension is seen before
    #for(j in 1:dim(dim_seen)[1]){
    # if(prod(dim_temp==dim_seen[j,1:dim_total])){
    #   L1<-dim_seen[j,dim_total+1]
    #   flag=1
    #   break
    # }
    #}
    ## Not seen
    if(flag==0){
      sim<-Estimate(Y, X, dim_temp, seed_ls[i], b, total_iter = 100, collect = FALSE)
      L1<-Likelihood(sim,data$Y)
      dim_seen<-rbind(dim_seen,c(dim_temp,L1))
    }
    u<-runif(1)
    R<-min(1,exp((1-b)*(L1-L0))*(q_01/q_10)*ratio)
    cat('Cur dim',core_dim,'\n')
    #cat('L0:', L0, 'L1:', L1, '\n')
    cat('Candidate dim',dim_temp,'\n')
    if(u<R){
      L0<-L1
      core_dim<-dim_temp
      unchange <- 0
      seed_star = seed_ls[i]
    }
    #if(i > I_burnin){
    #  sim_collect<-Estimate(data$Y, data$X, core_dim, seed, 1, total_iter = 200, collect = TRUE)
    #  result[[i-I_burnin]] = list(B_estimate = sim_collect$B_estimate, Sigma2_estimate = sim_collect$Sigma2_estimate)
    #}
    run_time = proc.time() - start_time
    cat('Running Time: ', round(run_time[3] / 3600,2), 'hours', '\n')
  }
  cat('Core Dimension:', core_dim, 'Start Collecting Posterior Samples \n')
  sim_collect<-Estimate(Y, X, core_dim, seed_star, 1, total_iter = n_iter, collect = TRUE)
  
  run_time = proc.time() - start_time
  cat('Total Running Time: ', round(run_time[3] / 3600,2), 'hours', '\n')
  return(sim_collect)
}
