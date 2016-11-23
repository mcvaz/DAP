AddBurnin = function(sims.array,burnin,n.thin){
  nchains=dim(sims.array)[2]
  end=dim(sims.array)[1]
  start=burnin/n.thin+1
  
  Results.sims.array=sims.array[start:end,,]
  
  sims.matrix=Results.sims.array[,1,]
  if(nchains>1){
    for(k in 2:nchains){
      sims.matrix=rbind(sims.matrix,Results.sims.array[,k,])
    }
    
  }
  
  Output.Matrix=mat.or.vec(dim(sims.matrix)[2],5)
  
  for(q in 1:dim(sims.matrix)[2]){
    Output.Matrix[q,1]=mean(sims.matrix[,q])
    Output.Matrix[q,2]=sd(sims.matrix[,q])
    Output.Matrix[q,3]=quantile(sims.matrix[,q],.025)
    Output.Matrix[q,4]=quantile(sims.matrix[,q],.975)	
    Output.Matrix[q,5]=length(sims.matrix[,q][sims.matrix[,q]>0])/length(sims.matrix[,q])
  }
  
  rownames(Output.Matrix)<-colnames(sims.matrix)
  colnames(Output.Matrix)<-c("mu.vect","sd.vect","2.5%","97.5%", "P>0")
  
  JAGS.Result=list(Burnin.sims.array=Results.sims.array,Burnin.sims.matrix=sims.matrix,Burnin.Summary=Output.Matrix)
  return(JAGS.Result)
}