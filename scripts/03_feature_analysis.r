library(ggfortify)
library(class)
library(e1071)
library(foreach)
library(doParallel)
library(corrplot)
library(Hmisc)
library(RColorBrewer)

# Load custom color palette
source("../data/input/COL2.txt")

pc_peptides_features <- read.table("../data/intermediate/protein_coding_peptides_features.txt")
pc_peptides_features <- pc_peptides_features[nchar(pc_peptides_features$V6) <= 100, ] 
rownames(pc_peptides_features) <- pc_peptides_features$V1
pc_peptides_features <- pc_peptides_features[, -1] 
colnames(pc_peptides_features) <- c("aromaticity", "instability_index", "isoelectric_point", "gravy", "peptide")

lncRNA_peptides_features <- read.table("../data/intermediate/lncrna_peptides_features.txt")
rownames(lncRNA_peptides_features) <- lncRNA_peptides_features$V1
lncRNA_peptides_features <- lncRNA_peptides_features[, -1] 
colnames(lncRNA_peptides_features) <- c("aromaticity", "instability_index", "isoelectric_point", "gravy", "peptide")



merge_datasets <- function(data1, data2) {
  common_ids <- intersect(rownames(data1), rownames(data2))
  data1 <- data1[common_ids, ]
  data2 <- data2[common_ids, ]
  combined_data <- cbind(data1, data2)
  return(combined_data)
}


pc_nucleotide_features <- read.table("../data/intermediate/protein_coding_nucleotide_features.txt")
lncrna_nucleotide_features <- read.table("../data/intermediate/lncrna_nucleotide_features.txt")

pc_combined_data <- merge_datasets(pc_nucleotide_features, pc_peptides_features)
lnc_combined_data <- merge_datasets(lncrna_nucleotide_features, lncRNA_peptides_features)

combined_data <- rbind(pc_combined_data, lnc_combined_data)
combined_data=combined_data[,-17]

gene_biotype <- ifelse(grepl("protein_coding", rownames(combined_data)), "protein_coding", "lncRNA")


# Function to perform statistical analysis and generate plots
perform_statistical_analysis <- function(data, gene_biotype, output_dir = ".") {
  feature_pvalues <- c()
  for (i in 1:ncol(data)) {
    feature_name <- colnames(data)[i]
    output_file <- file.path(output_dir, paste0(feature_name, ".pdf"))
    pdf(output_file)
    
    if (!(feature_name %in% c("ARE", "PAS", "Kozak"))) {
      boxplot(data[, i] ~ gene_biotype, ylab = feature_name, main = feature_name)
      res <- wilcox.test(data[, i] ~ gene_biotype)
      p <- res$p.value
      text(1.5, 0.9 * max(data[, i]), paste("P value =", p))
    } else {
      counts <- table(data[, i], gene_biotype)
      res <- chisq.test(data[, i], gene_biotype)
      p <- res$p.value
      mosaicplot(counts, xlab = paste("P value =", p), 
                 ylab = '', main = feature_name, col = 'orange')
    }
    dev.off()
    feature_pvalues <- c(feature_pvalues, p)
  }
  return(feature_pvalues)
}



# Function to perform correlation analysis
perform_correlation_analysis <- function(data, peptide_length, output_dir = ".") {
  feature_cor_pvalues <- c()
  feature_cor_r <- c()
  for (feature_name in colnames(data)) {
    cor_res <- cor.test(peptide_length, data[, feature_name])
    output_file <- file.path(output_dir, paste0(feature_name, ".cor.pdf"))
    pdf(output_file)
    plot(peptide_length, data[, feature_name], xlab = "Peptide length", ylab = feature_name, main = feature_name)
    text(median(peptide_length), 0.9 * max(data[, feature_name]), paste("p value =", cor_res$p.value))
    text(median(peptide_length), 0.8 * max(data[, feature_name]), paste("r =", cor_res$estimate))
    dev.off()
    feature_cor_pvalues <- c(feature_cor_pvalues, cor_res$p.value)
    feature_cor_r <- c(feature_cor_r, cor_res$estimate)
  }
  return(list(pvalues = feature_cor_pvalues, r_values = feature_cor_r))
}

# Function to perform PCA and save results
perform_pca <- function(data, output_dir = ".") {
  pca_result <- prcomp(data, center = TRUE, scale = TRUE)
  pca_data <- pca_result$x
  
  write.table(pca_data, file.path(output_dir, "features.PCA.txt"), sep = "\t", quote = FALSE)
  write.table(pca_result$rotation, file.path(output_dir, "features.PCA.rotation.txt"), sep = "\t", quote = FALSE)
  write.table(summary(pca_result)$importance, file.path(output_dir, "features.PCA.importance.txt"), sep = "\t", quote = FALSE)
  
  return(pca_data)
}

# Function to evaluate PCs
evaluate_PCs <- function(pca_data, pc_data, lnc_data, total_pc, output_dir = ".") {
  unlink(file.path(output_dir, "feature.accuracy.txt"))
  set.seed(123)
  all_pcs <- colnames(pca_data)
  
  for (each_feature in all_pcs) {
    print(each_feature)
    feature_data <- data.frame(pca_data[, each_feature])
    
    times <- 500
    registerDoParallel(30)
    
    result1 <- foreach(i = 1:times, .combine = rbind, .packages = c("e1071", "class")) %dopar% {
      random_train_sample <- c(sample(rownames(pc_data), total_pc), sample(rownames(lnc_data), total_pc))
      train <- data.frame(feature_data[random_train_sample, ])
      colnames(train) <- each_feature
      rownames(train) <- random_train_sample
      
      label <- ifelse(grepl("protein_coding", rownames(train)), "PC", "lncRNA")
      train <- data.frame(label = as.factor(label), train)
      
      test <- data.frame(feature_data[rownames(pc_data), ])
      rownames(test) <- rownames(pc_data)
      colnames(test) <- each_feature
      
      # KNN
      classifier_knn <- knn(train = train[, -1], test = test, cl = train$label, k = 10)
      knn_result <- data.frame(time = i, method = "KNN", sample = rownames(test), res = as.vector(classifier_knn))
      
      # Naive Bayes
      nb <- naiveBayes(label ~ ., data = train)
      nb_pred <- predict(nb, newdata = test)
      nb_result <- data.frame(time = i, method = "naiveBayes", sample = rownames(test), res = as.vector(nb_pred))
      
      # Logistic Regression
      train$label <- as.numeric(train$label == "PC")
      glm.fit <- glm(label ~ ., data = train, family = binomial)
      glm_pred <- predict(glm.fit, newdata = test, type = "response")
      glm_result <- data.frame(time = i, method = "glm", sample = rownames(test), res = glm_pred)
      
      return(rbind(knn_result, nb_result, glm_result))
    }
    stopImplicitCluster()
    
    # Calculate accuracy and save results
    accuracy <- calculate_accuracy(result1, total_pc)
    write.table(t(c(each_feature, accuracy)), file.path(output_dir, "feature.accuracy.txt"), append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
}

# Function to calculate accuracy
calculate_accuracy <- function(results, total_pc) {
  knn_res <- results[results$method == "KNN", ]
  nb_res <- results[results$method == "naiveBayes", ]
  glm_res <- results[results$method == "glm", ]
  
  knn_table <- table(knn_res$sample, knn_res$res)
  knn_pre <- rownames(knn_table)[knn_table[, "PC"] > knn_table[, "lncRNA"]]

  nb_table <- table(nb_res$sample, nb_res$res)
  nb_pre <- rownames(nb_table)[nb_table[, "PC"] > nb_table[, "lncRNA"]]

  glm_avg <- tapply(as.numeric(glm_res$res), glm_res$sample, mean)
  glm_pre <- names(glm_avg)[glm_avg > 0.7]
  
  common_predictions <- table(c(knn_pre, nb_pre, glm_pre))
  common_predictions <- names(common_predictions[common_predictions > 1])
  accuracy <- length(common_predictions) / total_pc
  return(accuracy)
}

# Perform statistical analysis
feature_pvalues <- perform_statistical_analysis(combined_data, gene_biotype, output_dir = "../data/intermediate")

# Extract significant features
significant_features <- colnames(combined_data)[feature_pvalues < 0.000001]

# Perform correlation analysis
correlation_results <- perform_correlation_analysis(pc_combined_data[, significant_features], nchar(pc_combined_data$peptide), output_dir = "../data/intermediate")

# Extract correlated features
correlated_features <- significant_features[correlation_results$pvalues < 0.05 & abs(correlation_results$r_values) > 0.2]
significant_features <- setdiff(significant_features, correlated_features)

# Perform PCA
pca_data <- perform_pca(combined_data[, significant_features], output_dir = "../data/intermediate")

# Evaluate PCs
evaluate_PCs(pca_data, pc_combined_data, lnc_combined_data, total_pc = nrow(pc_combined_data), output_dir = "../data/intermediate")
