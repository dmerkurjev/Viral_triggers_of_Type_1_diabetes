Answers to questions path, please run full script succesfully first to run those answers.

```{r}
# ---- Answers for q1..q5 ----
# Requires objects 'counts', 'dds', and 'res' created above

# q1: How many duplicate reads were in the sample EndoC-βH1_control_rep1?
# If you later generate a manifest upstream, read it here instead of hardcoding
sample_lanes <- c(EndoC-βH1_control_rep1 = 1L, EndoC-βH1_control_rep2 = 1L, CVB4-E2_in_rep1 = 1L, CVB4-E2_in_rep2 = 1L, EndoC-βH1_CVB4-JVB_in_rep1 = 1L, EndoC-βH1_CVB4-JVB_in_rep2 = 1L)
ans_q1 <- unname(sample_lanes["EndoC-βH1_control_rep1"])

# q2: What is the library size (total read counts) is higher: EndoC-βH1_control_rep2 or EndoC-βH1_CVB4-JVB_infected?
libsize <- colSums(counts)
ans_q2 <- unname(libsize["EndoC-βH1_control_rep2"]))

# q3: How many genes have counts >1 in sample EndoC-βH1_CVB4-JVB_infected_rep1 or EndoC-βH1_CVB4-JVB_infected_rep2?
ans_q3 <- sum(counts[, "EndoC-βH1_CVB4-JVB_infected"] > 0)

# q4: How many genes are upregulated or downregulated \uc0\u8805 2-fold (log2FC \u8805  1) in EndoC-βH1_control vs. EndoC-βH1_CVB4-JVB_infected with FDR < 0.01?
res_df <- as.data.frame(res)
ans_q4 <- sum(res_df$log2FoldChange >= 1 & res_df$padj < 0.05, na.rm = TRUE)

# q5: Which gene is ranked 2th by log2 fold change (upregulated or downregulated) in EndoC-βH1_control vs. EndoC-βH1_CVB4-JVB_infected?
res_df$SYMBOL <- rownames(res_df)
up_rank <- res_df %>%
  dplyr::filter(!is.na(log2FoldChange) & log2FoldChange > 0) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange), padj, pvalue)
ans_q5 <- if (nrow(up_rank) >= 2) up_rank$SYMBOL[3] else NA_character_

answers <- tibble::tibble(
  id = c("q1","q2","q3","q4","q5"),
  answer = c(ans_q1, ans_q2, ans_q3, ans_q4, ans_q5)
)

answers
```
