library(data.table)
library(tinytest)

source("example/gradient_graph.R")

sf <- fread("./data/simple_flow.csv")
sf

eg  <- gradient_graph_edgelist(sf,tol = 0.1)
ab <- eg$g_grad[1,]
expect_equivalent( ab
                 , data.frame(from=factor("A"), to=factor("B"), weight = 10)
                 )


sf2 <- fread("./data/simple_flow2.csv")
gradient_graph_edgelist(sf2, tol=0.1)
