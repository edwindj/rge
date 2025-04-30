library(data.table)
source("example/journal.R")
n_cg <- 20
rank_cg <- 5

# number of enterprises
n_v <- 500*n_cg
rank_v <- rank_cg + 1
max_outdegree <- 3

set.seed(1)
cg <- data.table(
  id = "",
  volume = 1, # runif(n_cg, 10, 1e4),
  rank = sample(rank_cg, size = n_cg, replace = TRUE),
  key = c("rank,volume")
)[, id := paste0("CG_", rank, ".", seq_len(.N)), by = rank]

suppliers <- data.table(
  id = paste0("E_", seq_len(n_v)),
  rank = sample(rank_v, size = n_v, replace = TRUE),
  turnover = rlnorm(n_v),
  key = c("rank,turnover")
)[, id := paste0("E_",rank, ".", seq_len(.N)), by = "rank"]

users <- suppliers[,.(to=id, rank = rank - 1, tt=turnover)][rank > 0,]
#setkey(users, rank)

edges <- suppliers[users,.(from=id, to, share = tt * turnover), on="rank",by = .EACHI]
edges <- edges[, .SD[sample(nrow(.SD), max_outdegree, replace = TRUE)], by = .(rank, from)]

edges <- edges[cg, .(from, to, cg = id, share, volume), on="rank", by =.EACHI]

edges <- edges[, .SD[sample(nrow(.SD), max_outdegree, replace = TRUE)], by = .(cg, from)]
edges[, share := share/sum(share), by = .(cg)]
edges[, volume := share * volume]

edges <- edges[, .(from,to,cg,volume)]
edges

library(ggplot2)

edges |>
  ggplot(aes(x = volume)) +
  geom_histogram(binwidth = 0.00001) +
  # coord_cartesian(xlim = c(0)) +
  labs(x = "", y = "", title = "volume of edges" ) +
  theme_journal()

edges |>
  dplyr::count(to) |>
  ggplot(aes(x=n)) +
  geom_histogram(binwidth = 1) +
  labs(x = "indegree", y = "count", title="") +
  theme_journal()

save_journal("fig/synth1_indegree_ent.pdf")

edges[, weight := volume]

# project

project_cg <- function(edges){
  use <- edges[,.(u = sum(weight)), by = .(id=to, use=cg)]
  supply <- edges[,.(s = sum(weight)), by = .(id=from, supply=cg)]
  supply <- supply[, s_ij := s/sum(s), by = id] # share of supplies

  su <- use[supply,,on="id", nomatch = NULL, allow.cartesian = TRUE] # only for supply-use connections

  # su[is.na(use), `:=`(u = 0, use = supply)]

  # both a supplier and user, so make it a net user/supplier
  self_supply <-  su[use == supply, .(u = u-s, s = s-u)][u < 0, u:=0][s < 0, s:=0]
  su[use == supply,]$s <- self_supply$s
  su[use == supply,]$u <- self_supply$u


  su <- su[s > 0,] # remove division by zero because of self supply
  su[, r_ijk := u * s_ij] # share the use of a cg

  cg <- su[, .(m = sum(r_ijk, na.rm = TRUE)), by = .(from=use,to=supply)]
  cg[, weight := m]
  cg[weight > 0, ]
}

edges_cg <- project_cg(edges)

edges_cg[, .(weight = sum(weight)), by = .(cg = to)] |>
  ggplot() +
  geom_histogram(aes(x = weight, fill = cg)
                , binwidth = 1)

g <- gradient_graph_edgelist(edges_cg, tol = 0.001)

s <- g$g_grad[, .(weight = sum(weight)), by = from]
s[, .(m = mean(weight), sd = sd(weight))]

l <- list()
for (f in seq(0, 1, by = 0.05)){
  if (f == 0) next
  message("adding in ", f, "random links")
  n_cyc <- (f*nrow(edges)) |> as.integer()

  cyc <- data.table( from = sample(suppliers$id, size = n_cyc, replace = TRUE)
                   , to = sample(suppliers$id, size = n_cyc, replace = TRUE)
                   , cg = sample(cg$id, size = n_cyc, replace=TRUE)
                   , volume = sample(edges$weight, size = n_cyc, replace = TRUE)
                   # , volume = 1e-4  #sample(edges$weight, size = n_cyc, replace = TRUE)
  )[from != to, ]

  n_cyc <- nrow(cyc)



  # E_opt <- sapply(seq_len(rank_v), function(i){
  #   sample(suppliers[rank==i, id], replace=TRUE, n_cyc)
  # })
  #
  # cyc$from <- E_opt[cbind(seq_len(n_cyc), cyc$rfrom)]
  # cyc$to <- E_opt[cbind(rev(seq_len(n_cyc)), cyc$rto)]
  # cyc$volume <- mean(edges$volume)
  # # cyc$volume <- sample(edges$volume, size = n_cyc, replace = TRUE)
  # exclude self loops
  #cyc <- cyc[from != to,]

  # revcyc <- cyc[, .(to=from, from=to, cg, volume)]
  # revcyc$from <- E_opt[cbind(rev(seq_len(n_cyc)), revcyc$rfrom)]
  # revcyc$to <- E_opt[cbind(seq_len(n_cyc), revcyc$rto)]
  #revcyc <- revcyc[from != to,]

  # cyc2 <- rbindlist(list(cyc[, .(from,to, cg, volume)], revcyc[, .(from,to, cg, volume)]))

  m <- rbindlist(list( edges[, .(from,to,cg, volume)]
                     , cyc
                     )
                )
  m <- m[, .(weight = sum(volume)), by = .(from,to,cg)]

  edges_cg <- m |>
    project_cg()

  W <- edges_cg[, sum(weight, na.rm = TRUE)]
  g <- gradient_graph_edgelist(edges_cg, tol = 0.001)
  gl <- g$g_grad[, .(w=sum(weight)), by = from]
  gl2 <- gl[, .(w = mean(w), sd = sd(w), f = f, t = sum(w)/W, edge_s = nrow(edges_cg))]
  l <- c(l, list(gl2))
}

b <- rbindlist(l)
b
# b |>
#   ggplot(aes(x = f, y = )) +
#   geom_point() +
#   theme_journal()

b |>
  ggplot(aes(x = f, y = w)) +
  geom_ribbon(aes(ymin = w - 2*sd, ymax = w + 2*sd), alpha=0.1) +
  geom_line() +
  theme_journal()


# generate a cg network
N_c <- 501

d_in <- c(
  rnorm(2*N_c/3, mean = 200, sd = 25),
  rnorm(N_c/3, mean = 100, sd = 25)
)

d_in <- c(0,cumsum(d_in/sum(d_in)))

M <- round(0.3 * N_c * (N_c - 1))

d_out <- rlnorm(N_c, 1)
d_out <- d_out/sum(d_out)
p = (N_c - 1)/M

for (i in 1:10){
  d_out <- pmin(d_out, p)
  d_out <- d_out/sum(d_out)
}

d_out <- c(0, cumsum(d_out))

from <- findInterval(runif(1.5*M), d_in)
to <- findInterval(runif(1.5*M), d_out)

d <- data.table(from, to)
d2 <- d[from != to,] |> head(M)
d2

d2[, .N, by = from] |>
  ggplot(aes(x = N)) + geom_histogram(binwidth = 5) +
  coord_cartesian(xlim = c(0, N_c)) +
  labs(y = "number of commodities", x = "outdegree (synthetic)", title="") +
  theme_journal()

save_journal("fig/C_synth_outdegree.pdf")

d2[, .N, by = to] |>
  ggplot(aes(x = N)) + geom_histogram(binwidth = 5) +
  labs(y = "number of commodities", x = "indegree (synthetic)", title="") +
  theme_journal()

save_journal("fig/C_synth_indegree.pdf")


d2[, weight := rlnorm(M, meanlog = log(5e3))]
d2 <- d2[, .(from, to, weight)]

d2 |>
  ggplot(aes(x = weight)) +
  scale_x_log10(guide="axis_logticks", labels=label_log()) +
  geom_histogram(binwidth = .05) +
  theme_journal()

g <- d2 |>
  gradient_graph_edgelist()

g1 <- d2
g1[, n := "C[0[]]"]

g2 <- g$g_grad
g2[, n := "C[0[nabla]]"]

g3 <- rbindlist(list(g1, g2))

library(scales)
g3 |>
  ggplot(aes(x = weight, linetype=n)) +
  # geom_freqpoly(binwidth = 0.05) +
  geom_density_journal() +
  # stat_density(
  #   aes(y = after_stat(count)),
  #   geom="line",
  #   position="identity"
  #   , alpha = 0.85
  # ) +
  # scale_y_continuous(expand = c(0,0)) +
  scale_x_log10(labels=label_log(), breaks=10^(1:7), limits=10^c(2,7), guide="axis_logticks") +
  labs(col = "", linetype = "") +
  scale_color_discrete(labels = parse_format()) +
  scale_linetype(labels = parse_format()) +
  labs(y = "number of transaction edges", x = "weight (synthetic)") +
  geom_vline(xintercept = 10^3, linetype=2, alpha=0.5) +
  theme_journal()

save_journal("fig/C_synth_weight_dist.pdf")


g3[, .(outdegree = .N), by = .(commodity=from, n)] |>
  ggplot(aes(x = outdegree, linetype = n)) +
  geom_density_journal(bw=2.5) +
  # stat_density(
  #   aes(y = after_stat(count)),
  #   geom="line",
  #   position="identity"
  #   , bw = 2
  # ) +
  # scale_y_continuous(expand = c(0,0)) +
  # geom_density(bw = 2) +
  labs(
    col = "",
    linetype = "",
    y = "number of commodities",
    x = "outdegree (synthetic)"
  ) +
  scale_color_discrete(labels = parse_format()) +
  scale_linetype(labels = parse_format()) +
  lims(x = c(0, 550)) +
  theme_journal()

save_journal("fig/C_synth_outdegree2.pdf")

g3[, .(indegree = .N), by = .(commodity=to, n)] |>
  ggplot(aes(x = indegree, linetype = n)) +
  geom_density_journal(bw = 2.5) +
  labs(
    col = "",
    linetype = "",
    y = "number of commodities",
    x = "indegree (synthetic)"
  ) +
  scale_color_discrete(labels = parse_format()) +
  scale_linetype(labels = parse_format()) +
  theme_journal()

save_journal("fig/C_synth_indegree2.pdf")

g3[, .(indegree = .N), by = .(commodity=to, n)] |>
  ggplot(aes(x = indegree, linetype = n)) +
  geom_density_journal(bw = 2.5) +
  labs(
    col = "",
    linetype = "",
    fill = "",
    y = "number of commodities",
    x = "indegree (synthetic)"
  ) +
  scale_linetype_manual(labels = parse_format(), values = c(1,2)) +
  # scale_fill_manual(labels = parse_format(), values = c("#cccccc", "#cccccc00")) +
  # scale_color_manual(labels = parse_format(), values = c("#cccccc", "#333333")) +
  theme_journal()


save_journal("fig/C_synth_indegree2.pdf")
