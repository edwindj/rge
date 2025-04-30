library(data.table)

gg_gg <- fread("./data/gg_projection.csv")
gg <- gg_gg[, .(to = supply_gg, from = use_gg, weight = waarde, n=N)]

source("example/gradient_graph.R")
grad <- gradient_graph_edgelist(gg)
setDT(grad$v)

goederen <- fread("data/goederengroep.csv", encoding="Latin-1")
setkey(goederen, gg)
grad$v[nchar(id) == 6, id := paste0("0", id)]
grad$v$description <- goederen[(grad$v$id), description ]

View(grad$v)


grad$v
setDT(grad$e_gradient)
e_gradient <- grad$e_gradient
e_gradient[, from := as.character(from)]
e_gradient[, to := as.character(to)]
e_gradient[nchar(from) == 6, from := paste0("0", from)]
e_gradient[nchar(to) == 6, to := paste0("0", to)]
e_gradient[, from_desc := goederen[(from), description]]
e_gradient[, to_desc := goederen[(to), description]]
View(e_gradient)


