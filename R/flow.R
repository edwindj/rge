library(data.table)
d <- fread("./data/flow.csv")


library(ggplot2)

d |>
  ggplot(aes(y = from, x = to, fill = value)) +
  geom_tile()


library(ggraph)

id <- unique(c(d$from, d$to))

nodes <- data.frame(id = id)
edges <- d


g <- igraph::graph_from_data_frame(d)
g |>
  create_layout("linear") |>
  ggraph() +
  geom_edge_arc(aes(width = value, color = value), arrow = arrow()) +
  scale_fill_viridis() +
  geom_node_label(aes(label=name))
