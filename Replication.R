# add regularization & 50 times
# Load required packages
library(farff)
library(e1071)
library(pls)
library(pROC)
library(LiblineaR)
library(caret)

# Load Yeast ARFF Dataset
df <- as.data.frame(readARFF("/Users/yanfeng/Downloads/file2754771351f4.arff"))
X_full <- scale(df[, 1:103])
Y_full <- as.matrix(data.frame(lapply(df[, 104:ncol(df)], function(x) as.integer(as.logical(x)))))

# Metrics
compute_metrics <- function(Y_true, Y_pred) {
  L <- ncol(Y_true)
  f1s <- numeric(L); aucs <- numeric(L)
  tp <- fp <- fn <- 0
  for (l in 1:L) {
    true <- Y_true[, l]; pred <- Y_pred[, l]
    tp_l <- sum(true == 1 & pred == 1)
    fp_l <- sum(true == 0 & pred == 1)
    fn_l <- sum(true == 1 & pred == 0)
    prec <- ifelse(tp_l + fp_l == 0, 0, tp_l / (tp_l + fp_l))
    rec  <- ifelse(tp_l + fn_l == 0, 0, tp_l / (tp_l + fn_l))
    f1s[l] <- ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
    aucs[l] <- tryCatch(roc(true, pred)$auc, error = function(e) NA)
    tp <- tp + tp_l; fp <- fp + fp_l; fn <- fn + fn_l
  }
  micro_p <- tp / (tp + fp); micro_r <- tp / (tp + fn)
  list(
       macro_f1 = mean(f1s, na.rm = TRUE),
       micro_f1 = 2 * micro_p * micro_r / (micro_p + micro_r),
       auc = mean(aucs, na.rm = TRUE))
}

# SVM Evaluation Wrapper
evaluate <- function(Z_train, Y_train, Z_test, Y_test) {
  L <- ncol(Y_train)
  Y_pred <- matrix(0, nrow = nrow(Z_test), ncol = L)
  for (l in 1:L) {
    if (length(unique(Y_train[, l])) < 2) {
      Y_pred[, l] <- 0
      next
    }
    model <- LiblineaR(Z_train, factor(Y_train[, l], levels = c(0, 1)), type = 1, cost = 100)
    Y_pred[, l] <- as.integer(as.character(predict(model, Z_test)$predictions))
  }
  compute_metrics(Y_test, Y_pred)
}

# Safe matrix inversion via eigenvalue thresholding
safe_inverse <- function(mat, epsilon = 1e-3) {
  eig <- eigen(mat, symmetric = TRUE)
  eig$values <- pmax(eig$values, epsilon)
  eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors)
}

# SPPCA / S2PPCA EM with regularization
sppca_em <- function(X_labeled, Y_labeled, X_unlabeled = NULL, K = 10, max_iter = 1000, tol = 1e-5) {
  X_labeled <- as.matrix(X_labeled)
  Y_labeled <- as.matrix(Y_labeled)
  if (!is.null(X_unlabeled)) X_unlabeled <- as.matrix(X_unlabeled)
  
  X_all <- if (is.null(X_unlabeled)) X_labeled else rbind(X_labeled, X_unlabeled)
  mu_x <- colMeans(X_all)
  X_all <- scale(X_all, center = mu_x, scale = FALSE)
  Y_labeled <- scale(Y_labeled, center = TRUE, scale = FALSE)
  
  N_l <- nrow(X_labeled); N <- nrow(X_all); M <- ncol(X_all); L <- ncol(Y_labeled)
  Wx <- matrix(rnorm(M * K, 0, 0.1), M, K)
  Wy <- matrix(rnorm(L * K, 0, 0.1), L, K)
  sigma_x <- 1e-5; sigma_y <- 1e-5
  Ez <- matrix(0, N, K); Wx_old <- Wx
  
  for (iter in 1:max_iter) {
    Wy_old <- Wy; sigma_x_old <- sigma_x; sigma_y_old <- sigma_y
    Ezz <- array(0, c(K, K, N))
    
    A_l <- safe_inverse(t(Wx) %*% Wx / sigma_x^2 + t(Wy) %*% Wy / sigma_y^2 + diag(K))
    for (n in 1:N_l) {
      bx <- t(Wx) %*% X_labeled[n, ] / sigma_x^2
      by <- t(Wy) %*% Y_labeled[n, ] / sigma_y^2
      Ez[n, ] <- A_l %*% (bx + by)
      Ezz[,,n] <- A_l + Ez[n,] %*% t(Ez[n,])
    }
    
    if (!is.null(X_unlabeled)) {
      A_u <- safe_inverse(t(Wx) %*% Wx / sigma_x^2 + diag(K))
      for (n in (N_l + 1):N) {
        bx <- t(Wx) %*% X_all[n, ] / sigma_x^2
        Ez[n, ] <- A_u %*% bx
        Ezz[,,n] <- A_u + Ez[n,] %*% t(Ez[n, ])
      }
    }
    
    C <- apply(Ezz, c(1, 2), sum)
    Wx <- t(X_all) %*% Ez %*% solve(C)
    Wy <- t(Y_labeled) %*% Ez[1:N_l, , drop = FALSE] %*% solve(C)
    sigma_x <- sqrt(mean((X_all - Ez %*% t(Wx))^2))
    sigma_y <- sqrt(mean((Y_labeled - Ez[1:N_l, ] %*% t(Wy))^2))
    
    if (
      max(abs(Wx - Wx_old)) < tol &&
      max(abs(Wy - Wy_old)) < tol &&
      abs(sigma_x - sigma_x_old) < tol &&
      abs(sigma_y - sigma_y_old) < tol
    ) break
    
    Wx_old <- Wx
  }
  
  return(list(Wx = Wx, Wy = Wy, Ez = Ez, sigma_x = sigma_x, mu_x = mu_x))
}

posterior_z <- function(X_new, Wx, sigma_x, mu_x) {
  A_inv <- safe_inverse(t(Wx) %*% Wx + sigma_x^2 * diag(ncol(Wx)))
  t(A_inv %*% t(Wx) %*% t(scale(X_new, center = mu_x, scale = FALSE)))
}
# experiment with 50 repetitions
set.seed(694)
num_reps <- 50
k_values <- c(5, 10, 20)
methods <- c("FULL", "PCA", "PLS", "SPPCA", "S2PPCA")
metrics <- c("macro_f1", "micro_f1", "auc")

# Initialize results storage
results_array <- array(0, dim = c(length(methods), length(metrics), length(k_values), num_reps),
                       dimnames = list(methods, metrics, paste0("K", k_values), 1:num_reps))

for (rep in 1:num_reps) {
  cat("\n=== Repetition", rep, "/", num_reps, "===\n")
  set.seed(rep) 
  
  # Create train/test split for each repetition
  label_ids <- 1:ncol(Y_full)
  train_indices <- c()
  
  # Select 5 positive samples per class with safety checks
  for (l in label_ids) {
    pos_idx <- which(Y_full[, l] == 1)
    if (length(pos_idx) >= 5) {
      selected <- sample(pos_idx, 5)
    } else {
      selected <- pos_idx  # Use all available if <5
    }
    train_indices <- union(train_indices, selected)
  }
  
  X_labeled <- X_full[train_indices, ]
  Y_labeled <- Y_full[train_indices, ]
  X_unlabeled <- X_full[-train_indices, ]
  Y_unlabeled <- Y_full[-train_indices, ]
  
  for (k_idx in seq_along(k_values)) {
    K <- k_values[k_idx]
    cat("  Running K =", K, "\n")
    
    # PCA
    tryCatch({
      pca_model <- prcomp(X_labeled, rank. = K)
      Z_train_pca <- predict(pca_model, X_labeled)
      Z_test_pca <- predict(pca_model, X_unlabeled)
      res_pca <- evaluate(Z_train_pca, Y_labeled, Z_test_pca, Y_unlabeled)
      results_array["PCA", , k_idx, rep] <- unlist(res_pca)
    }, error = function(e) cat("PCA failed:", e$message, "\n"), silent=TRUE)
    
    # PLS
    tryCatch({
      pls_model <- plsr(Y_labeled ~ X_labeled, ncomp = K)
      Z_train_pls <- scores(pls_model)
      Z_test_pls <- predict(pls_model, newdata = X_unlabeled, type = "scores")[, 1:K]
      res_pls <- evaluate(Z_train_pls, Y_labeled, Z_test_pls, Y_unlabeled)
      results_array["PLS", , k_idx, rep] <- unlist(res_pls)
    }, error = function(e) cat("PLS failed:", e$message, "\n"), silent=TRUE)
    
    # SPPCA
    tryCatch({
      sppca <- sppca_em(X_labeled, Y_labeled, K = K)
      Z_train_sppca <- sppca$Ez[1:nrow(X_labeled), ]
      Z_test_sppca <- posterior_z(X_unlabeled, sppca$Wx, sppca$sigma_x, sppca$mu_x)
      res_sppca <- evaluate(Z_train_sppca, Y_labeled, Z_test_sppca, Y_unlabeled)
      results_array["SPPCA", , k_idx, rep] <- unlist(res_sppca)
    }, error = function(e) cat("SPPCA failed:", e$message, "\n"), silent=TRUE)
    
    # S2PPCA
    tryCatch({
      s2ppca <- sppca_em(X_labeled, Y_labeled, X_unlabeled, K = K)
      Z_train_s2 <- s2ppca$Ez[1:nrow(X_labeled), ]
      Z_test_s2 <- posterior_z(X_unlabeled, s2ppca$Wx, s2ppca$sigma_x, s2ppca$mu_x)
      res_s2ppca <- evaluate(Z_train_s2, Y_labeled, Z_test_s2, Y_unlabeled)
      results_array["S2PPCA", , k_idx, rep] <- unlist(res_s2ppca)
    }, error = function(e) cat("S2PPCA failed:", e$message, "\n"), silent=TRUE)
    
    # FULL
    tryCatch({
      res_full <- evaluate(X_labeled, Y_labeled, X_unlabeled, Y_unlabeled)
      results_array["FULL", , k_idx, rep] <- unlist(res_full)
    }, error = function(e) cat("FULL failed:", e$message, "\n"), silent=TRUE)
  }
}

# Calculate mean and standard deviation across repetitions
final_results <- list(
  mean = apply(results_array, 1:3, mean, na.rm = TRUE),
  sd = apply(results_array, 1:3, sd, na.rm = TRUE)
)

# Print formatted results
for (k in dimnames(final_results$mean)[[3]]) {
  cat("\nResults for", k, "components:\n")
  cat("Method       Macro-F1 (SD)    Micro-F1 (SD)    AUC (SD)\n")
  for (method in dimnames(final_results$mean)[[1]]) {
    m <- round(final_results$mean[method, , k], 4)
    s <- round(final_results$sd[method, , k], 4)
    cat(sprintf("%-8s  %5.4f (±%.4f)  %5.4f (±%.4f)  %5.4f (±%.4f)\n",
                method,
                m["macro_f1"], s["macro_f1"],
                m["micro_f1"], s["micro_f1"],
                m["auc"], s["auc"]))
  }
}













