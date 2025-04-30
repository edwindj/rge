flip_edges <- function(x, flip){
  for (n in colnames(x)[-c(1:2)]){
    x[[n]][flip] <- -x[[n]][flip]
  }
  from_flip <- x$to[flip]
  x$to[flip] <- x$from[flip]
  x$from[flip] <- from_flip
  x
}

normalize_edges <-  function(x){
  flip_edges(x, unclass(x$from) > unclass(x$to))
}

denormalize_edges <-  function(x){
  flip_edges(x, x$weight < 0)
}
