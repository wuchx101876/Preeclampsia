

remove(list = ls())

## GSE149437

library(PRROC)
  
gse149437 <- readRDS("./raw_data/GSE149437/use_data.rds")

df <- gse149437
df$status <- ifelse(df$type == "PE",0,1)



set.seed(123)

pr_obj <- pr.curve(
  scores.class0 = df$BACH1[df$status == 0],
  scores.class1 = df$BACH1[df$status == 1],
  curve = FALSE
)

auprc_value <- pr_obj$auc.integral
auprc_value




set.seed(123)

n_boot <- 1000
auprc_boot <- numeric(n_boot)

for (i in seq_len(n_boot)) {
  idx <- sample(seq_len(nrow(df)), replace = TRUE)
  df_boot <- df[idx, ]
  
  pr_boot <- pr.curve(
    scores.class0 = df_boot$BACH1[df_boot$status == 0],
    scores.class1 = df_boot$BACH1[df_boot$status == 1],
    curve = FALSE
  )
  
  auprc_boot[i] <- pr_boot$auc.integral
}

auprc_ci <- quantile(auprc_boot, probs = c(0.025, 0.975))
auprc_ci




pr_curve_full <- pr.curve(
  scores.class0 = df$BACH1[df$status == 0],
  scores.class1 = df$BACH1[df$status == 1],
  curve = TRUE
)


auprc <- pr_curve_full$auc.integral

plot(
  pr_curve_full$curve[,1],
  pr_curve_full$curve[,2],
  type = "l",
  lwd = 2,
  col = "#d62728",
  xlab = "Recall",
  ylab = "Precision",
  main = paste0(
    "Precision–Recall curve for BACH1 (GSE149437)\n",
    "AUPRC = ", round(auprc, 3)," (95% CI: 0.75–0.95)"
  )
)




## GSE48424


gse48424 <- readRDS("./result/Validation/GSE48424_model_data.rds")

df <- gse48424

set.seed(123)

pr_obj <- pr.curve(
  scores.class0 = df$BACH1[df$status == 0],
  scores.class1 = df$BACH1[df$status == 1],
  curve = FALSE
)

auprc_value <- pr_obj$auc.integral
auprc_value


set.seed(123)

n_boot <- 1000
auprc_boot <- numeric(n_boot)

for (i in seq_len(n_boot)) {
  idx <- sample(seq_len(nrow(df)), replace = TRUE)
  df_boot <- df[idx, ]
  
  pr_boot <- pr.curve(
    scores.class0 = df_boot$BACH1[df_boot$status == 0],
    scores.class1 = df_boot$BACH1[df_boot$status == 1],
    curve = FALSE
  )
  
  auprc_boot[i] <- pr_boot$auc.integral
}

auprc_ci <- quantile(auprc_boot, probs = c(0.025, 0.975))
auprc_ci




pr_curve_full <- pr.curve(
  scores.class0 = df$BACH1[df$status == 0],
  scores.class1 = df$BACH1[df$status == 1],
  curve = TRUE
)


auprc <- pr_curve_full$auc.integral

plot(
  pr_curve_full$curve[,1],
  pr_curve_full$curve[,2],
  type = "l",
  lwd = 2,
  col = "#d62728",
  xlab = "Recall",
  ylab = "Precision",
  main = paste0(
    "Precision–Recall curve for BACH1 (GSE48424)\n",
    "AUPRC = ", round(auprc, 3), " (95% CI: 0.52–0.93)"
  )
)


