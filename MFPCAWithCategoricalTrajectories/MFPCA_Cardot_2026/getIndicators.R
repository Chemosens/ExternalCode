getIndicators <- function(state, times, times.out=seq(0,1,length=100))
{
  npoints <- length(times.out)
  if(any(diff(times)<0)){stop("times should be ordered")}
  indicatrices <- matrix(0,length(unique(state)),npoints)
  rownames(indicatrices) <- unique(state)
  for(i in 1:(length(state)-1))
  {     
    current_state <- state[i]
    current_time <- times[i]
    ind.times <- times.out >= current_time & times.out < times[i+1]
    indicatrices[state[i],ind.times] <- 1  
  }
  indicatrices[state[length(state)],times.out>=times[length(state)]] <- 1   
  return(indicatrices)
}
