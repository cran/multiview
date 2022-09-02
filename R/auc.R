auc=function(y,prob,w){
  if(missing(w))
    survival::concordance(y~prob)$concordance
  else
    survival::concordance(y~prob,weights=w)$concordance
}
