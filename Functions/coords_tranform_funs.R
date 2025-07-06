
hypersphere_coords_fun <- function(X){
  #assume that X is an n x m matrix 
  # radius 1, angles between 0 and 1 (pi/2)
  if(is.null(dim(X))) X <- matrix(X,nrow = 1)
  
  n_dim <- ncol(X) -1
  
  out <- matrix(nrow = nrow(X),
                ncol = n_dim)
  
  Xt <- sqrt(X)
  for(i in seq(n_dim - 1)){
    temp <- X[,-c(1:i)]
    temp <- sqrt(rowSums(temp))
    out[,i] <- atan2(temp,Xt[,i])
  } 
  
  out[,n_dim] <- 2*atan2(Xt[,n_dim + 1], 
                         Xt[,n_dim] + sqrt(rowSums(Xt[,-c(1:i)]^2)))
  return(out)
  
}

hsc_to_cart_fun <- function(x){
  cos_x <- cospi(x/pi)
  sin_x <- sinpi(x/pi)
  
  out <- rep(0, length(x) + 1)
  
  out[1] <- cos_x[1]
  
  for(i in 2:length(x)) out[i] <- cos_x[i] * prod(sin_x[1:(i-1)])
  
  out[length(x) + 1] = prod(sin_x)
  
  out^2
  
}

hsc_to_cart_mat_fun <- function(angles){
  
  cos_x <- cospi(angles/pi)
  sin_x <- sinpi(angles/pi)
  
  out <- cos_x[,1]
  psin <- 1
  
  for(i in 2:ncol(angles)) {
    psin <- psin * sin_x[,i-1]
    out <- cbind(out,
                 cos_x[,i] * psin)
  }
  
  out<- cbind(out, psin * sin_x[,i])
  
  return(out^2)
  
}