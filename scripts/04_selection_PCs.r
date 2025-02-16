library(ggfortify)
library(class)
library(e1071)
library(foreach)
library(doParallel)

# Load PCA data and feature accuracy data
pca_data <- read.table("../data/intermediate/features.PCA.txt")
feature_accuracy <- read.table("../data/intermediate/feature.accuracy.txt", sep = " ")
feature_accuracy <- feature_accuracy[order(feature_accuracy$V2, decreasing = TRUE), ]
ordered_features <- feature_accuracy$V1

# Function to perform cross-validation and evaluate models
perform_cv <- function(pca_data, ordered_features, cv = 10, times = 500, output_file = "../data/intermediate/feature.combine.accuracy.txt") {
  unlink(output_file)  # Clear the output file
  
  for (num_features in 2:length(ordered_features)) {
    print(paste("Evaluating feature combination:", num_features))
    
    # Select the top N features
    selected_features <- ordered_features[1:num_features]
    feature_data <- pca_data[, selected_features, drop = FALSE]
    
    # Split data into protein-coding (PC) and lncRNA samples
    pc_samples <- rownames(feature_data)[grep("protein_coding", rownames(feature_data))]
    pc_samples <- sample(pc_samples)  # Shuffle PC samples
    lnc_samples <- rownames(feature_data)[grep("lncRNA", rownames(feature_data))]
    
    # Split PC samples into CV chunks
    chunks <- split(pc_samples, cut(seq_along(pc_samples), cv, labels = FALSE))
    cv_accuracy <- c()
    
    for (chunk_name in names(chunks)) {
      test_samples <- chunks[[chunk_name]]
      train_samples <- setdiff(pc_samples, test_samples)
      
      # Parallel computation
      registerDoParallel(30)
      model_results <- foreach(i = 1:times, .combine = rbind, .packages = c("e1071", "class")) %dopar% {
        temp_results <- c()
        random_train_samples <- c(train_samples, sample(lnc_samples, length(train_samples)))
        train_data <- data.frame(feature_data[random_train_samples, , drop = FALSE])
        
        # Create labels
        labels <- ifelse(grepl("protein_coding", rownames(train_data)), "PC", "lncRNA")
        train_data <- data.frame(label = as.factor(labels), train_data)
        
        test_data <- feature_data[test_samples, , drop = FALSE]
        
        # KNN
        knn_pred <- knn(train = train_data[, -1, drop = FALSE], test = test_data, cl = train_data$label, k = 10)
        temp_results <- rbind(temp_results, data.frame(time = rep(i, nrow(test_data)), method = "KNN", sample = rownames(test_data), res = as.vector(knn_pred)))
        
        # Naive Bayes
        nb_model <- naiveBayes(label ~ ., data = train_data)
        nb_pred <- predict(nb_model, newdata = test_data)
        temp_results <- rbind(temp_results, data.frame(time = rep(i, nrow(test_data)), method = "naiveBayes", sample = rownames(test_data), res = as.vector(nb_pred)))
        
        # Logistic Regression
        train_data$label <- ifelse(train_data$label == "PC", 1, 0)
        glm_model <- glm(label ~ ., data = train_data, family = binomial)
        glm_pred <- predict(glm_model, newdata = as.data.frame(test_data), type = "response")
        temp_results <- rbind(temp_results, data.frame(time = rep(i, nrow(test_data)), method = "glm", sample = rownames(test_data), res = glm_pred))
        return(temp_results)
      }
      stopImplicitCluster()
      
      # Process results
      knn_results <- model_results[model_results$method == "KNN", ]
      nb_results <- model_results[model_results$method == "naiveBayes", ]
      glm_results <- model_results[model_results$method == "glm", ]
      
      # Extract predictions
      knn_predictions <- rownames(table(knn_results$sample, knn_results$res))[table(knn_results$sample, knn_results$res)[, "PC"] > table(knn_results$sample, knn_results$res)[, "lncRNA"]]
      nb_predictions <- rownames(table(nb_results$sample, nb_results$res))[table(nb_results$sample, nb_results$res)[, "PC"] > table(nb_results$sample, nb_results$res)[, "lncRNA"]]
      glm_predictions <- names(which(tapply(as.numeric(glm_results$res), glm_results$sample, mean) > 0.7))
      
      # Combine predictions
      combined_predictions <- table(c(glm_predictions, knn_predictions, nb_predictions))
      combined_predictions <- names(combined_predictions)[combined_predictions > 1]
      accuracy <- length(combined_predictions) / length(test_samples)
      cv_accuracy <- c(cv_accuracy, accuracy)
    }
    
    # Write results to file
    write.table(t(c(paste(selected_features, collapse = "_"), mean(cv_accuracy),paste(cv_accuracy,collapse=";"))), output_file, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
}

# Perform cross-validation
perform_cv(pca_data, ordered_features)

