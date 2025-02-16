# Load necessary libraries
library(foreach)
library(doParallel)
library(class)
library(e1071)

# Load PCA data and feature importance data
pca_data <- read.table("../data/intermediate/features.PCA.txt")
importance_data <- read.table("../data/intermediate/feature.combine.accuracy.txt", sep = " ")
best_feature_combination <- importance_data$V1[importance_data$V2 == max(importance_data$V2)]
best_feature_combination <- strsplit(best_feature_combination, " ")[[1]][1]
best_features <- strsplit(best_feature_combination, "_")[[1]]

# Register parallel backend
registerDoParallel(50)

# Subset PCA data with the best features
selected_data <- pca_data[, best_features]
times <- 500

# Get protein-coding and lncRNA sample IDs
pc_samples <- rownames(selected_data)[grep("protein_coding", rownames(selected_data))]
total_pc <- length(pc_samples)
lnc_samples <- rownames(selected_data)[grep("lncRNA", rownames(selected_data))]

# Perform parallel computation
model_results <- foreach(i = 1:times, .combine = rbind, .packages = c("e1071", "class")) %dopar% {
  # Print progress every 10 iterations
  if (i %% 10 == 0) {
    print(paste("Completed", i, "iterations"))
  }
  
  # Randomly sample training data
  random_train_samples <- c(sample(pc_samples, total_pc), sample(lnc_samples, total_pc))
  train_data <- selected_data[random_train_samples, ]
  
  # Create labels
  labels <- ifelse(grepl("protein_coding", rownames(train_data)), "PC", "lncRNA")
  train_data <- data.frame(label = as.factor(labels), train_data)
  
  # Test data (samples not in training data)
  test_data <- selected_data[!(rownames(selected_data) %in% rownames(train_data)), ]
  
  # KNN
  classifier_knn <- knn(train = train_data[, -1], test = test_data, cl = train_data$label, k = 10)
  knn_result <- data.frame(time = i, method = "KNN", sample = rownames(test_data), res = as.vector(classifier_knn))
  
  # Naive Bayes
  nb_model <- naiveBayes(label ~ ., data = train_data)
  nb_pred <- predict(nb_model, newdata = test_data)
  nb_result <- data.frame(time = i, method = "naiveBayes", sample = rownames(test_data), res = as.vector(nb_pred))
  
  # Logistic Regression
  train_data$label <- as.numeric(train_data$label == "PC")
  glm_model <- glm(label ~ ., data = train_data, family = binomial)
  glm_pred <- predict(glm_model, newdata = as.data.frame(test_data), type = "response")
  glm_result <- data.frame(time = i, method = "glm", sample = rownames(test_data), res = glm_pred)
  
  # Combine results
  return(rbind(knn_result, nb_result, glm_result))
}
stopImplicitCluster()

knn_res <- model_results[model_results$method == "KNN", ]
nb_res <- model_results[model_results$method == "naiveBayes", ]
glm_res <- model_results[model_results$method == "glm", ]
  
knn_table <- table(knn_res$sample, knn_res$res)
knn_pre <- rownames(knn_table)[knn_table[, "PC"] > knn_table[, "lncRNA"]]

nb_table <- table(nb_res$sample, nb_res$res)
nb_pre <- rownames(nb_table)[nb_table[, "PC"] > nb_table[, "lncRNA"]]

glm_avg <- tapply(as.numeric(glm_res$res), glm_res$sample, mean)
glm_pre <- names(glm_avg)[glm_avg > 0.7]
  
common_predictions <- table(c(knn_pre, nb_pre, glm_pre))
common_predictions <- names(common_predictions[common_predictions > 1])

write.table(common_predictions, "../data/output/CodingLncRNAsORF.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

