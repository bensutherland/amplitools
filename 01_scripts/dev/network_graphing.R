# Visualize network of relatedness 
# Sutherland Bioinformatics
# 2023-08-02

# Install.packages("igraph")
library("igraph")

# Set working directory
setwd("~/Documents/00_sbio/GBMF_UBC_Pacific_oyster/amplicon_panel_USDA_pilot/amplitools/")

# Read in input
sib.df <- read.delim2(file = "03_results/parent_fs_broodstock_pw_logl_5.txt"
                      , header = T, sep = "\t"
                      )

nrow(sib.df) # how many interactions? 

# Reduce length of names
sib.df$D2_indiv <- gsub(pattern = "G00", replacement = "", x = sib.df$D2_indiv)
sib.df$D1_indiv <- gsub(pattern = "G00", replacement = "", x = sib.df$D1_indiv)

# # Create a graph
# # for example
# g <- make_graph(edges = c(1,2, 1,5), n=10, directed = FALSE)
# summary(g)
# plot(g)
# 
# g <- make_graph('Zachary')
# plot(g)
# 
# # #todo: need to assign name attribute as well
# 
# 
# test.df
# 
# g <-        make_graph(edges = c("G005498", "G005499",  "G001945", "G001951")
#                    #, n = 4
#                    , directed = F
#                    )
# 
# g <- make_graph(edges = c("G005498", "G005499"))
# plot(g)
# 
# g <- graph_from_literal(G005498--G005499)
# plot(g)
# 
# g <- graph_from_data_frame(test.df, directed = FALSE, vertices = NULL)
# plot(g)
# 
# 
# 
# # need a vector of sequential interactions
# interactions <- c(paste0(test.df$D2_indiv, ",", test.df$D1_indiv, collapse = ","))
# str(interactions)
# make_graph(edges = interactions)
# 
# 
# interactions <- paste0(sib.df$D2_indiv, ",", sib.df$D1_indiv, collapse = ",")
# interactions <- gsub(pattern = " ", replacement = "_", x = interactions)
# g <- make_graph(edges = interactions)


sibs.plot <- graph_from_data_frame(sib.df, directed = F)
print(sibs.plot, e=TRUE, v=TRUE)

plot.igraph(sibs.plot, label.cex = 0.2, vertex.color = "grey"
            , vertex.size = 10
            , vertex.label.color = "black"
            , vertex.frame.color="gray"
            , vertex.label.cex=0.8
            , vertex.label.dist=2
            , edge.curved=0.2
            )

# What about one with higher required logl cutoff? 
head(sib.df)
sib.df <- sib.df[sib.df$logl_ratio > 20, ]
nrow(sib.df)
sibs.plot <- graph_from_data_frame(sib.df, directed = F)

pdf(file = "03_results/relatedness_network_parent_sibs_2023-08-02.pdf", width = 5, height = 5)
plot.igraph(sibs.plot
            , vertex.color = "grey"
            , vertex.label.color = "black"
            , vertex.frame.color="gray"
            , vertex.size = 5
            , label.cex = 1
            , vertex.label.cex=0.7
            , vertex.label.dist=1
            , edge.curved=0.2
)
dev.off()
