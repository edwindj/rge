source("example/gradient_graph.R")
library(data.table)
library(ggplot2)
library(scales)
source("example/journal.R")

ggg <- fread("data/gg_projection.csv")
gg_desc <- fread("data/gg_agt_iot.csv")

labels <- sub("^0+", "", gg_desc$gg_agt)
lookup <- gg_desc$desc_agt |> setNames(labels)

library(igraph)
g <- graph_from_data_frame(ggg)

ggg |>
  ggplot(aes(x = waarde)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10(labels=label_log(), guide = "axis_logticks") +
  theme_journal()

use_cg <- ggg[, .N, by = .(gg= use_gg)][order(N, decreasing = TRUE), ]
use_cg[, cg := lookup[gg |> as.character()]]
use_cg

#library(ggcbs)
library(ggplot2)

ggplot(use_cg) +
  geom_histogram(aes(x=N), binwidth = 5) +
  labs(title="", y = "number of commodities", x = "indegree C (DPN2018)") +
  coord_cartesian(xlim=c(0, NA), expand = FALSE) +
  theme_journal()

save_journal("fig/C_DPN_indegree.pdf")

supply_cg <- ggg[, .N, by = .(gg= supply_gg)][order(N, decreasing = TRUE), ]
supply_cg[, cg := lookup[gg |> as.character()]]
supply_cg

ggplot(supply_cg) +
  geom_histogram(aes(x=N), binwidth = 5) +
  labs(y = "number of commodities", x = "outdegree C (DPN2018)") +
  theme_journal()

save_journal("fig/C_DPN_outdegree.pdf")

ggg[order(waarde, decreasing = TRUE),]
ggg[, use_desc := lookup[use_gg |> as.character()]]
ggg[, supply_desc := lookup[supply_gg |> as.character()]]

ggg2 <- ggg[, .(to=supply_gg, from = use_gg, weight=waarde, N)]


source("example/gradient_graph.R")
hodge <- hodge_graph_edgelist(ggg2, tol=.1)
d <- calc_div(hodge$g_grad)
sum(abs(d))

hd_grad <- gradient_graph_edgelist(ggg2, tol=.1, max_N = 1e2)
grad <- hd_grad$g_grad

library(igraph)
g <- graph_from_data_frame(grad)
message("is_dag: ", is_dag(g))

# <<<<<<< HEAD
# grad[, from_desc := lookup[from]]
# grad[, to_desc := lookup[to]]
#
# =======
# grad <- grad[g_r > 0, .(from, to, w = g_r)]
grad[, from_desc := lookup[from |> as.character()]]
grad[, to_desc := lookup[to |> as.character()]]

library(scales)

# >>>>>>> 0a45992f0072dd4bd5f0d6b7c092f2375ff110f7
grad |>
  ggplot(aes(x = weight)) +
  geom_histogram(binwidth = .1, color="NA") +
  geom_vline(xintercept = 10^3, linetype="dashed") +
  scale_x_log10(labels=label_log(), breaks=10^(0:7)) +
  labs(x = "weight", y = "", title="weight distribution") +
  theme_minimal()

ggsave("fig/grad_value_dist.pdf", width = 7, height = 5, scale = 0.5)

grad <- grad[ weight >= 1e3, ]

#View(grad)

G <- graph_from_data_frame(grad)
d_out <- distances(G, mode="out", weights = NA)
is.na(d_out) <- !is.finite(d_out)

dl <- data.table(
  from   = rownames(d_out),
  to =  colnames(d_out) |> rep(each=nrow(d_out)),
  distance = as.numeric(d_out)
)

dl <- dl[distance>0, ]
dl <- dl[order(distance, decreasing = TRUE), ]

dl[, from_desc := ..lookup[from |> as.character()]]
dl[, to_desc := ..lookup[to |> as.character()]]

grad2 <- data.table(hd_grad$g_grad)[, .(from, to, w = weight)]
setkey(grad2, from,to)
dl$w <- grad2[.(dl$from, dl$to),w]
dl

dl[w > 0,][order(distance, decreasing = T),]

dl[, .N, by = .(distance)] |>
ggplot(aes(x = distance, y = N)) +
  geom_line() +
  geom_point() +
  labs(y = "number of supply chains", x = "path length (DPN2018)") +
  scale_x_continuous(breaks = 1:5, limits = c(1,5)) +
  scale_y_log10(labels = label_log(), limits = 10^(c(1,5)), guide="axis_logticks") +
  theme_journal()

save_journal("fig/C_DPN_grad_path_length.pdf")


fwrite(dl, "data/distances.csv")
# View(dl)

pot <- as.data.table(hd_grad$v)

pot[, desc := ..lookup[id |> as.character()]]
pot <- pot[order(pot, decreasing = TRUE),]

lookup_p <- pot$pot |> setNames(pot$id)

dl[, value := lookup_p[from |> as.character()] - lookup_p[to |> as.character()]]
# View(dl)

# gg2[, to_desc := ..lookup[as.character(to)]]
# gg2[, from_desc := ..lookup[as.character(from)]]
#
# View(gg2)
#View(pot)
f <- normalize_edges(ggg2)

U <- f[, .(w2 = min(abs(weight))), by = .(from, to)]

sum(U$w2)

flow <- fread("data/flow.csv")
flow

e <- gradient_graph_edgelist(ggg2)
#View(e$v)

pot <- e$v
# pot$naam <- gg$desc_agt[match(pot$id, gg$gg_agt)]
# View(pot)

library(scales)
# vergelijken grad en org
c_dpn <- rbindlist(list(
  ggg[, .(value = waarde, from = supply_gg, to = use_gg, what="C[]")],
  grad[, .(value = weight,from, to, what = "C[nabla]")]
))

c_dpn |>
  ggplot(aes(x = value, linetype = what)) +
  geom_density_journal(bw = .02) +
#  geom_freqpoly(binwidth = 0.05) +
  scale_x_log10(
    labels = label_log(), breaks=10^(1:7), limits=c(100,NA), guide = "axis_logticks"
               ) +
  labs(col = "", linetype = "") +
  scale_color_discrete(labels = parse_format()) +
  scale_linetype(labels = parse_format()) +
  labs(y = "number of transaction edges", x = "weight (DPN2018)") +
  geom_vline(xintercept = 10^3, linetype=2, alpha=0.5) +
  # annotate("text", x = 10^4, y = 10, "test") +
  theme_journal()

save_journal("fig/C_DPN_weight_dist.pdf")

c_dpn

c_dpn[, .(indegree = .N), by = .(commodity=to, what)] |>
  ggplot(aes(x = indegree, linetype = what)) +
  geom_density_journal(bw = 2.5) +
  labs(
    col = "",
    linetype = "",
    fill = "",
    y = "number of commodities",
    x = "indegree (DPN2018)"
  ) +
  scale_linetype_manual(labels = parse_format(), values = c(1,2)) +
  # scale_fill_manual(labels = parse_format(), values = c("#cccccc", "#cccccc00")) +
  # scale_color_manual(labels = parse_format(), values = c("#cccccc", "#333333")) +
  theme_journal()

save_journal("fig/C_DPN_indegree2.pdf")

c_dpn[, .(outdegree = .N), by = .(commodity=from, n = what)] |>
  ggplot(aes(x = outdegree, linetype = n)) +
  geom_density_journal(bw = 2.5) +
  labs(
    col = "",
    linetype = "",
    fill = "",
    y = "number of commodities",
    x = "outdegree (DPN2018)"
  ) +
  scale_linetype_manual(labels = parse_format(), values = c(1,2)) +
  # scale_fill_manual(labels = parse_format(), values = c("#cccccc", "#cccccc00")) +
  # scale_color_manual(labels = parse_format(), values = c("#cccccc", "#333333")) +
  theme_journal()

save_journal("fig/C_DPN_outdegree2.pdf")
