########## Reshape Tensor to Matrix ###########
array_to_matrix<-function(Tensor,rdims,cdims){
  # rdims = which dimensions in tensor become row in new matrix
  # cdims = which dimensions in tensor become column in new matrix
  dimens<-dim(Tensor)
  r<-dimens[rdims] # Number of elements in this dimension
  c<-dimens[cdims] 
  if(prod(r)*prod(c)!=prod(dimens)){
    cat("Dimensions not compatible")
    return(Tensor)
  }else{
    new_tensor<-aperm(Tensor,c(rdims,cdims))
    new_matrix<-matrix(as.vector(new_tensor),prod(r),prod(c))
    return(new_matrix)
  }
}

########## Tensor and Tensor Multiplication (Contracted Product) ##########
tensor_tensor<-function(A,B,adims,bdims){
  ## adims = contracted dims in A, bdims = contracted dims in B
  ## the dimensions should match
  sizeA<-dim(A);modeA<-c(1:length(sizeA))
  sizeB<-dim(B);modeB<-c(1:length(sizeB))
  ## If adims=bdims=0, then outer product of A,B tensors are computed
  if(adims==0 && bdims==0){
    amatrix<-as.vector(A)
    bmatrix<-as.vector(B)
    newmatrix<-amatrix%*%t(bmatrix)
    newtensor<-array(as.vector(newmatrix),c(sizeA,sizeB))
    return(newtensor)
  }
  ## Check if the Adim, Bdim match
  if(prod(sizeA[adims]==sizeB[bdims])){
    amatrix<-array_to_matrix(A,modeA[-adims],modeA[adims])
    bmatrix<-array_to_matrix(B,modeB[bdims],modeB[-bdims])
    newmatrix<-amatrix%*%bmatrix
    newtensor<-array(as.vector(newmatrix),c(sizeA[modeA[-adims]],sizeB[modeB[-bdims]]))
    return(newtensor)
  }else{
    return(NULL)
  }
}

########## Tensor times Matrix ###########
tensor_matrix<-function(A,B,mod){
  ## mode 'mod' of tensor A match second dimension of matrix B
  sizeA<-dim(A);modeA<-c(1:length(sizeA))
  Amatrix<-array_to_matrix(A,mod,modeA[-mod])
  newmatrix<-B%*%Amatrix
  newmode<-c(dim(B)[1],sizeA[-mod]); permutation<-modeA
  permutation[1:mod]<-permutation[1:mod]+1; permutation[mod]<-1
  newtensor<-array(as.vector(newmatrix),newmode)
  tens<-aperm(newtensor,permutation)
  return(tens)
}

########## Matrix outer product ##########
matrix_outer<-function(A,B){
  dimA<-dim(A);dimB<-dim(B)
  temp<-A[1,1]*B
  for(k in 2:dimA[2]){
    temp<-cbind(temp,A[1,k]*B)  
  }
  result<-temp
  for(i in 2:dimA[1]){
    temp<-A[i,1]*B
    for(j in 2:dimA[2]){
      temp<-cbind(temp,A[i,j]*B)
    }
    result<-rbind(result,temp)
  }
  return(result)
}

########## Dimension Candidates from Neighbour ##########
Find_candidate<-function(core_dim,K){
  candidate<-matrix(0,1,length(core_dim))
  n<-1
  for(i in c(1:length(core_dim))){
    if(core_dim[i]==2){
      candidate[n,i]<-1
      n<-n+1
      candidate<-rbind(candidate,rep(0,length(core_dim)))
      next
    }
    if(core_dim[i]==K[i]){
      candidate[n,i]<- -1
      n<-n+1
      candidate<-rbind(candidate,rep(0,length(core_dim)))
      next
    }
    candidate[n,i]<-1
    n<-n+1
    candidate<-rbind(candidate,rep(0,length(core_dim)))
    candidate[n,i]<- -1
    n<-n+1
    candidate<-rbind(candidate,rep(0,length(core_dim)))
  }
  return(candidate[1:(n-1),])
}