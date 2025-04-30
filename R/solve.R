solve_svd <- function(A, b, tol = .Machine$double.eps){
  # could also use sparsesvd!

  A_inv <- with(svd(A), {
    is_pos <- d > max(tol * d[1L], 0)
    # so that 1/d is zero
    d[!is_pos] <- Inf
    # alternatively, may be faster, drop all !is_pos columns
    v %*% ((1/d) * t(u))
  })
  as.numeric(A_inv %*% b)
}

solve_qr <- function(A, b, tol = .Machine$double.eps){
  as.numeric(
    solve( qr(A, tol = tol)
           , b
    )
  )
}
