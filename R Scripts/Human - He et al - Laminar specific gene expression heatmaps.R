#Load required libraries
library(reshape2)

#Create lists of layer-specific genes

Zeng_layer1_gene_list <- Zeng_Layer1$gene_symbol
Zeng_layer2_gene_list <- Zeng_Layer2$gene_symbol
Zeng_layer3_gene_list <- Zeng_Layer3$gene_symbol
Zeng_layer4_gene_list <- Zeng_Layer4$gene_symbol
Zeng_layer5_gene_list <- Zeng_Layer5$gene_symbol
Zeng_layer6_gene_list <- Zeng_Layer6$gene_symbol

#Subset values with above lists

He_Layer1_scaled <- He_Layer1_scaled %>%
  select(Zeng_layer1_gene_list)
He_Layer2_scaled <- He_Layer2_scaled %>%
  select(Zeng_layer2_gene_list) 
He_Layer3_scaled <- He_Layer3_scaled %>%
  select(Zeng_layer3_gene_list)   
He_Layer4_scaled <- He_Layer4_scaled %>%
  select(Zeng_layer4_gene_list)
He_Layer5_scaled <- He_Layer5_scaled %>%
  select(Zeng_layer5_gene_list) 
He_Layer6_scaled <- He_Layer6_scaled %>%
  select(Zeng_layer6_gene_list) 

#Heatmaps

transform_tibble <- function(x) {
  t(x) %>%
    melt() %>%
    as_tibble()
}

#Heatmap of layer 1 specific genes
t_vHL1_test <- transform_tibble(values_He_Layer1_test)
ggplot(data = t_vHL1_test, mapping = aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

#Heatmap of layer 2 specific genes
t_vHL2_test <- transform_tibble(values_He_Layer2_test)
ggplot(data = t_vHL2_test, mapping = aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

#Heatmap of layer 3 specific genes
t_vHL3_test <- transform_tibble(values_He_Layer3_test)
ggplot(data = t_vHL3_test, mapping = aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

#Heatmaps of layer 4 specific genes
t_vHL4_test <- transform_tibble(values_He_Layer4_test)
ggplot(data = t_vHL4_test, mapping = aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

#Heatmap of layer 5 specific genes
t_vHL5_test <- transform_tibble(values_He_Layer5_test)
ggplot(data = t_vHL5_test, mapping = aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

#Heatmap of layer 6 specific genes
t_vHL6_test <- transform_tibble(values_He_Layer6_test)
ggplot(data = t_vHL6, mapping = aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

