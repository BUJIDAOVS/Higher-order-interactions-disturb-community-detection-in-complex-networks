# Higher-order Interactions Disturb Community Detection in Complex Networks

This repository contains code and data processing steps used in the study "Higher-order interactions disturb community detection in complex networks." The study focuses on the impact of higher-order interactions on community detection within complex networks.

## Network Data Preprocessing

1. `bip2one.m` - Converts bipartite networks to unilayer networks.
2. `high_order_nw.m` - Transforms bipartite networks into unilayer and multilayer networks.

## Model Construction

1. `build_modle.m` - Constructs the cross-community cooperation model. In this model, higher-order interactions are more likely to occur across communities.
2. `build_modle_rand.m` - Establishes a random model where interactions occur freely, independent of the orders.

## Analysis of Imbalance in Number of Interactions and Edges

1. `statistic.m` - Demonstrates the imbalance in the number of interactions (i.e., papers) and edges per order.

## High-order Interference in Community Partitioning in the Model

1. `high_order_nosie.m` - Script for plotting Figure 3, revealing the interference of high-order interactions in community detection within the model.

## Higher Probability of Cross-Community Cooperation in Higher Orders

1. `order_weight.m` - Compares empirical networks and models to validate cross-community behaviors in high-order interactions in real systems.

## Removing High-order Information Benefits Community Partitioning

1. `pacs_deal_data.m` - Processes the APS dataset, correlating authors with PACS codes of papers.
2. `pacs_sim_big_author.m` - Validates the consistency and effectiveness of field identification before and after removing high-order interactions.

---

For details on the methods used in these scripts, as well as the analysis of results and conclusions, please refer to the corresponding sections of the research paper.
