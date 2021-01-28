### this function will first create a generation based table with states
s0 = 7
i0 = 7
n = s0+i0
R =3.5
pFSfast = function(R,x,s0in,i0in)
  {
   #infection probability
    pinf =function(s){s*R/(s*R + n)};
    
   #create final size distribution for each s0 and i0 combination
   gen.mat.list = list(list())
   for(m in 1:length(s0in))
   {
     
     gen.mat <- matrix(rep(0, (s0+i0+1)*(s0+1)),nrow= (s0+1))
     gen.mat[s0+1,i0+1]<- 1;
    #fill first
    for(j in c((s0+1):1)){gen.mat[s0+1,j-1]<- gen.mat[s0+1,j]*(1-pinf(s0)) }
    gen.mat
  #next rows added by infection
  for(i in c(1:(s0+1))){
    for(j in c((s0+i0+1):1)){
      #enter state by infection (only for j > 2, j = 1 no infections, j = 2 is 1 infection)
      if(j>2){
        gen.mat[s0+1-i,j]<-gen.mat[s0+1-i,j] + gen.mat[s0+2-i,j-1]*pinf(s0+1-i)
       };
    if(j < s0+i0+1){
      #enter state by recovery
      gen.mat[s0+1-i,j]<- gen.mat[s0+1-i,j] + gen.mat[s0+1-i,j+1]*(1-pinf(s0-i));
      }
     
    }
   }
   }
}
sum(gen.mat[,1])

gen.tab <- function(R,s0,i0){
  gen.mat <- matrix(rep(0, s0+i0*s0),nrow= s0);
  gen.mat[i0]<- 1;
  p 
  for(j in 2:s0)
  {
    next.gen <- 
      gen.mat <- rbind(gen.mat,
                       next.gen)
  }
  
}
