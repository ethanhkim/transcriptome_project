# Refer to 2. Human - Preprocessing - Zeng et al for preprocessing steps

library(magrittr)
library(dplyr)

# Compare V1 and V2 expression levels

Zeng_V1_V2_comparison <- Zeng_dataset_updated %>%
  mutate(Expression_comparison = ifelse(
    test = (V1_expression_level == V2_expression_level),
    yes = "SAME",
    no = "DIFFERENT"
  )) %>%
  mutate(Pattern_comparison = ifelse(
    test = (V1_expression_pattern == V2_expression_pattern),
    yes = "SAME",
    no = "DIFFERENT"
  ))

# Compare number of different vs. same expression levels

Zeng_V1_V2_comparison %>%
  group_by(Expression_comparison) %>%
  summarize(n())

## 73 genes show different expression levels - 2 show NA ##
## NA genes: MAOA, PHGDH

# Compare number of different. same expression patterns 

Zeng_V1_V2_comparison %>%
  group_by(Pattern_comparison) %>%
  summarize(n())

## 277 genes show different expression patterns ##

# Compare genes which had donor differences for expression level and pattern

test <- Zeng_V1_V2_comparison %>%
  filter(h_donor_difference == "diff", .keep_all = TRUE)

Zeng_V1_V2_comparison %>%
  filter(h_donor_difference == "diff", .keep_all = TRUE) %>%
  group_by(Pattern_comparison) %>%
  summarize(n())

Zeng_V1_V2_comparison %>%
  filter(h_donor_difference == "diff", .keep_all = TRUE) %>%
  group_by(Expression_comparison) %>%
  summarize(n())
