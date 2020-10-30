extract_segment = function(epath, tpath, tmin, tmax) {
  # Extract the portion of the complete trajectory that covers the time range 
  # [tmin, tmax].
  # 
  # Parameters:
  #  epath - edges visited in complete trajectory
  #  tpath - transition times in complete trajectory
  
  # indices of transitions that immediately precede and follow [tmin, tmax]
  istart = max(which(tpath <= tmin))
  iend = min(which(tmax <= tpath))
  
  # correct end index if last path transition occurs before tmax
  if(is.infinite(iend)) {
    iend = length(tpath)
  }
  
  
  # extract and return trajectory segment that covers [tmin, tmax]
  inds = istart:iend
  res = list(epath = epath[inds], tpath = tpath[inds])
  
  if(tail(res$tpath, 1) < tmax) {
    res$tpath = c(res$tpath, tmax)
  }
    
  res
}