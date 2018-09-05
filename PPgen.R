PPgen=function(lambda){
  X=0
  Sum=0
  flag=0
  while (flag==0){
    E=-log(runif(1))
    Sum=Sum+E
    if (Sum < lambda) { X=X+1} else { flag=1}
  }
  return(X)
}