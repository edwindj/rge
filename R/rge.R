restricted_gradient_extraction <- function( x
                                     , tol = .Machine$double.eps
                                     , method=c("svd")
                                     , max_N = 10
                                     , div_scale = 1
                                     , verbose = FALSE
){
  method <- match.arg(method)
  # node id to factor...
  if (is.factor(x$from)){
    x$from <- as.character(x$from)
  }
  if (is.factor(x$to)){
    x$to <- as.character(x$to)
  }
  # above code is needed to make sure, concatenating is going well.
  f <- factor(c(x$from, x$to))
  l <- levels(f)

  # recoding node ids
  x$from <- factor(x$from, levels = l)
  x$to <- factor(x$to, levels = l)

  # remove self links
  x <- x[x$from != x$to,]

  # normalize el: always from < to (by negating weight).
  xn <- flip_edges(x, unclass(x$from) > unclass(x$to))

  # browser()
  # take the net flow + range
  xn <- xn[
    , .( weight = sum(weight)
         , w_max = max(c(weight, 0))
         , w_min = min(c(weight, 0))
    )
    , by = .(from, to)
  ]
  # x

  v <- data.table(id = l, key = "id")

  # get divergence
  v$div <- calc_div(xn)
  v$pot <- calc_pot(xn)
  # v$div <- tapply(x$weight, x$from, sum, default = 0) -
  #          tapply(x$weight, x$to,   sum, default = 0)
  #
  # v$div <- as.numeric(v$div)

  # L <- laplacian( xn
  #               , n = length(l)
  #               , sym = TRUE
  #               )

  # x <- xn[,]
  # # start, we copy the weight into w_r (remaining weight)
  # x[, w_r := weight]
  # x[, w_g := 0] # we start with an empty gradient
  #
  # # subset the network, so that the edges that are
  # # fully gradient are removed
  # # xn <- x[w_r > 0, .(from, to, weight = w_r, w_max, w_min)]
  # pot <- calc_pot(x, div = div)
  #
  # v$pot <- pot
  #
  # # constraint 1: only allow edges in positive direction or with
  # # a negative weight that is smaller than current collected gradient
  # # that is a correction
  # x[, g := (pot[from] - pot[to]) |> pmax(-w_g)]
  #
  # # constraint 2: gradient must be smaller than weight
  # x[, w_g := pmin(w_g + g, weight)]
  #
  # # the current non-gradient component
  # x[, w_r := weight - w_g]
  #
  # x[, w_g := pmin(weight)]
  #
  #
  # # potential is translation invariant, so make lowest potential 0.
  # # v$pot <- v$pot - min(v$pot)
  #
  # # the gradient flow is the difference in potential
  # e_gradient <- within(xn, {
  #   w_g <- v$pot[from] - v$pot[to]
  # })

  x_r <- xn[,]
  l <- list()
  div_r <- calc_div(x_r)
  div_total <- sum(abs(div_r))
  ds <- div_total
  for (i in seq_len(max_N)){
    div_sum <- sum(abs(div_r))

    if (min(ds) < div_sum){
      message("divergence increasing")
      break
    }


    ds <- c(ds, div_sum)
    message("loop ", i, ", div left: ", div_scale * div_sum / div_total)

    if (all(div_r*div_r < tol)){
      message("no divergence left")
      break
    }

    pot <- calc_pot(x_r, div = div_r)
    x_r[, g := pot[from] - pot[to]]

    # e_test[, g_r := ifelse(g<min, min, ifelse(g>max, max, g))]
    x_r[, remove := FALSE]
    x_r[g < 0, remove := g - w_min < tol]
    x_r[g > 0, remove := w_max - g < tol]

    x_r[g < 0, g_r := pmax(w_min, g)]
    x_r[g_r < -tol, w_min := w_min - g_r]

    x_r[g >= 0, g_r := pmin(w_max, g)]
    x_r[g_r >= tol, w_max := w_max - g_r]
    # x_r[g > 0, w_min := 0] # block the opposite path

    x_r[, w_r := w_min + w_max]

    if (verbose){
      print(list(div_r = div_r, x_r = x_r))
    }

    l <- c(l, list(
      x_r[, .(from,to, weight, g_r, w_min, w_max)]
    ))

    if (sum(x_r$remove) == 0){
      message("no removal")
      break
    }

    # r <- x_r[remove == TRUE, ]
    # if (nrow(r) == 0){
    #   break
    # }

    # subtract rest component
    div_r <- div_r - calc_div(x_r[, .(from,to,weight = g_r)])
    x_r[g_r < 0, w_max := 0] # block the opposite path
    x_r[g_r > 0, w_min := 0] # block the opposite path

    x_r <- x_r[remove == FALSE,]
    if (nrow(x_r) == 0){
      break
    }
    x_r[, weight := w_r]
    # calculate the residual div
  }

  plot(ds, type = "b")

  # denormalize gradient
  g_grad <- rbindlist(l)[, .(weight = sum(g_r)), by = .(from, to)]
  g_grad <- flip_edges(g_grad, g_grad$weight < 0)
  g_grad <- g_grad[weight > tol,]

  g_comp <- rbindlist(
    list( x[, .(from,to,weight)]
          , g_grad[, .(from, to, weight = -weight)]
    )
  )

  g_comp <- g_comp[, .(weight = sum(weight)), by = .(from,to)]
  g_comp <- g_comp[ weight > tol,]

  list( v          = v
        , g_grad     = g_grad
        , g_comp     = g_comp
        , div_scale  = div_scale * min(ds)/div_total
  )
}
