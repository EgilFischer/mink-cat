### this function will first create a generation based table with states and then gives the final size distribution as 
#output
distFSfast <- function(R,s0in,i0in ,r0in = NULL)
  {
    #set size of the population 
    nin = s0in+i0in;
    if(!is.null(r0in)) 
           nin = s0in+i0in+r0in;
    
    #infection probability
    pinf =function(s,n){s*R/(s*R + n)};
    
     #create final size distribution for each s0 and i0 combination
     max.cases <- max(s0in);
     gen.mat.out = matrix(0,ncol = (max.cases+1),nrow = length(s0in));
     for(m in 1:length(s0in))
       {
       s0= s0in[m];
       i0= i0in[m];
       n = nin[m];
       #set the state matrix
       gen.mat <- matrix(rep(0, (s0+i0+1)*(s0+1)),nrow= (s0+1));
       gen.mat[s0+1,i0+1]<- 1;
       #fill first row 
       for(j in c((i0+1):1)){
         gen.mat[s0+1,j-1]<- gen.mat[s0+1,j]*(1-pinf(s0,n)) 
         }
          
        #next rows added by infection
        for(i in c(1:(s0+1))){
          for(j in c((s0+i0+1):1)){
            #enter state by infection (only for j > 2, j = 1 no infections, j = 2 is 1 infection)
            if(j>2){
              gen.mat[s0+1-i,j]<-gen.mat[s0+1-i,j] + gen.mat[s0+2-i,j-1]*pinf(s0+1-i,n)
             };
          if(j < s0+i0+1){
            #enter state by recovery
            gen.mat[s0+1-i,j]<- gen.mat[s0+1-i,j] + gen.mat[s0+1-i,j+1]*(1-pinf(s0-i,n));
            }
           
          }
        }
    #store the final size distribution
    gen.mat.out[m,1:(s0+1)] = gen.mat[,1];
     }
     return(gen.mat.out)
}

system.time(tst<-distFSfast(3.5,c(49,19,5,4),c(1,1,1,1)))
tst

###Function to determine the probability of more extreme values given R####
# r = R, x = final number of cases, s0 = initial susceptibles, i0 = initial infectious, comp = the type of extreme 
rm(pExtremes)
pExtremes<-  function(r,x,s0in,i0in,comp = `<`){
  #create all possible outcomes of these transmission experiments
  #this means all possibilities between 0 and s0 contact infection (hence s0 + 1 options per trial)
  out <- matrix(ncol = length(s0),nrow = prod(s0+1))
  #repeat the outcome as many times as the previous trial possibilities
  for(k in c(1:length(s0)))
  {
    #check how often to repeat the same number given previous trials
    repetition <- ifelse(k > 1,prod(s0[1:k-1]+1),1)
    #put it in the matrix
    out[,k]<- matrix(sapply(c(0:s0[k]),FUN = function(x){rep(x,repetition)}), ncol = 1, nrow =prod(s0+1))[,1]
  }
  #produce final size distribution for value r
  final.size.dist <- distFSfast(r,s0in,i0in);
  #define function for this distribution for the probability of a certain number of cases x in each of the trials
  pFS<- function(v,m){return(prod(prod(mapply(function(i,j)final.size.dist[i,j],c(1:length(s0in)),v))))}
  #select the extremes by selecting those for which the total number of cases
  #and the sum of x are given by the comparison "comp"  thus either <,>,<= or >= 
  #calculate for each extreme the probability of the Final Size under the hypothesis R = r
  #and sum all probabilities of these extreme outcomes
  #If only one extreme outcome possible: 
  if(is.null(dim(out[comp(apply(out,1,sum),sum(x)),]))){
    return(pFS(out[comp(apply(out,1,sum),sum(x)),],m))}
  #else
  return(c(sum(apply(out[comp(apply(out,1,sum),sum(x)),],1,function(ext){pFS(ext,m)}))))
}
