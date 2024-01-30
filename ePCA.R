# The following code is adapted from the original ePCA algorithm as presented in MATLAB in the ePCA paper.

# Function to estimate the scaled true spike
standard_spiked_forward<- function(ell,gamma){
  k<- length(ell)
  lambda<- rep(0, k)
  cos_right<- rep(0, k)
  cos_left<- rep(0, k)
  v<- rep(0, k)
  gamma_minus<- (1-sqrt(gamma))^2; gamma_plus<- (1+sqrt(gamma))^2
  
  for (i in 1:k){
    if ((ell[i] < gamma^{1/2}) & (ell[i] > -gamma^{1/2})){
      lambda[i]<- (1+gamma^{1/2})^2
      cos_right[i]<- 0; cos_left[i]<- 0
    } else{
      lambda[i]<- (1+ell[i])*(1+gamma/ell[i])
      cos_right[i]<- (1-gamma/ell[i]^2)/(1+gamma/ell[i]);
      cos_left[i] = (1-gamma/ell[i]^2)/(1+1/ell[i]);
    }
    x<- lambda[i]
    im_mult = 1
    if ((x > gamma_minus) & (x < gamma_plus)){
      im_mult=1i
    }
    v[i]<- 1/(2*x)*(-(1+x-gamma)+im_mult*(abs((1+x-gamma)^2-4*x))^{1/2})
  }
  m = 1/gamma*v- (1-1/gamma)/lambda;
  return(v)
}

# Function to denoise the estimated S
denoise<- function(est_S, Y_bar1, D_n1, Y){
  X<- est_S %*% pinv(D_n1 + est_S) %*% t(Y) + D_n1%*% pinv(D_n1 + est_S) %*% matrix(Y_bar1) %*% rep(1, n)
  return(X)
}

# Main function for ePCA
epca<- function(obs_array){
  m1<- dim(obs_array)[1]; m2<- dim(obs_array)[2]; n<- dim(obs_array)[3]
  dim(obs_array)<- c(m1*m2, n)
  
  obs_array<- t(obs_array)
  Y<- obs_array
  Y_bar<- apply(Y, 2,mean)
  S<- (t(Y)-Y_bar)%*%t(t(Y)-Y_bar)/n
  D_n<- diag(Y_bar * (1-Y_bar)) ### binomial
  S_d<- S-D_n
  S_h<- sqrt(D_n^{1/2}) %*% S_d %*%sqrt(D_n^{1/2})
  svd_Sh<- svd(S_h); w<- svd_Sh$u; wt<- svd_Sh$v; lambda<- svd_Sh$d
  
  r=round(m1/2)*23
  gamma = m1*m2/n
  white_eval = lambda^2; E = lambda[1:r]^2;
  white_shr_eval<- vector("numeric", length(E))
  for (i in 1:length(E)){
    if (E[i] > ((1+sqrt(gamma))^2)){
      white_shr_eval[i]<- (E[i]+1-gamma+sqrt((E[i]+1-gamma)^2-4*E[i]))/2-1
    } else {
      white_shr_eval[i]<- 1+sqrt(gamma)-1
    }
  }
  
  S_h_eta<- w[,1:r] %*% diag(white_shr_eval) %*% wt[1:r,]
  S_he<- D_n^{1/2} %*% S_h_eta %*% D_n^{1/2} 
  
  
  lambda<- c(lambda[1:r], rep(0, m1*m2-r))
  S_h_eta<- w %*% diag(lambda) %*% wt
  S_he<- D_n^{1/2} %*% S_h_eta %*% D_n^{1/2} 
  recolor_eigend<- eigen(S_he); 
  recolor_eval<- sort(recolor_eigend$values, decreasing = TRUE); 
  recolor_v<- recolor_eigend$vectors[, order(recolor_eigend$values, decreasing = TRUE)]
  
  c2<- standard_spiked_forward(white_shr_eval,gamma)
  s2<- 1- c2
  tau<- (sum(D_n)*white_shr_eval)/(m1*m2*recolor_eval)
  alpha<- rep(0, r)
  for (i in 1:r){
    if (c2[i] > 0){
      alpha[i]<-  (1-s2[i]*tau[i])/c2[i]
    } else {
      alpha[i]<- 1
    }
  }
  
  white_shr_eval_scaled<- alpha*recolor_eval
  eigval<- diag(white_shr_eval_scaled[1:r])
  eigvec<- recolor_v[, 1:r]
  covar_est<- eigvec %*% eigval %*% t(eigvec)
  
  L_vec<- denoise(covar_est, Y_bar, D_n, Y)
  L_pred<- apply(L_vec, 1, mean)
  L_mat<- Re(matrix(L_pred, ncol=m1, nrow=m2, byrow=T))
  
  return(L_mat)
}