require(tidyverse)
require(ggplot2)
require(pbapply)
require(ape)
require(familiar)
# Function to read in and get per-species gene copy number
# for each orthogroup
get_per_spp_og_counts <-
  function(results_dir = NULL,
           out_dir = NULL) {
    orthogroup_dir <-
      list.files(
        path = paste0(results_dir, "orthofinder/complete_dataset/"),
        pattern = "Results_Inflation",
        full.names = TRUE
      )
    # read in per-species gene copy number per gene family
    og_counts <-
      read.table(
        paste0(orthogroup_dir, "/Orthogroups/Orthogroups.GeneCount.tsv"),
        header = TRUE,
        check.names = FALSE
      )
    colnames(og_counts) <- gsub("\\..*", "", colnames(og_counts))
    # And calculate the number of species in each gene family
    og_counts$NumSpecies <-
      rowSums(og_counts[, -c(1, ncol(og_counts))] > 0)
    # Write out to file if an output directory is provided
    if (!is.null(out_dir)) {
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      # And write out to file
      write.table(
        og_counts,
        file = paste0(out_dir, "og_copy_num_per_spp.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
    }
    return(og_counts)
  }

# Define a function to read in the event counts per-node (terminal and internal)
get_og_events_per_node <-
  function(i, per_spp_og_events = per_spp_og_events, tree) {
    # Read in table of per-species event counts for this gene family
    tmp <- read.table(per_spp_og_events[i], check.names = FALSE)
    colnames(tmp) <- c("node", "s", "sl", "d", "t")
    tmp <- tmp[which(tmp$node %in% tree$tip.label), ]
    return(tmp)
  }

# function to fit a Lasso regression model, fitting based on transformed data:
fit_lasso_counts <-
  function(count_data,
           tree,
           response,
           response_name,
           predictor_name,
           cluster_method = "hclust") {
    # Conduct a phylogenetic transformation of count data
    # Data Preparation
    phyvcv <- ape::vcv(tree)
    l <- chol(phyvcv)
    # Ensure species are ordered the same as in the tree
    data <-
      count_data[, match(tree$tip.label, colnames(count_data))]
    # Transform the observed count variables
    transf_data <- as.matrix(data) %*% t(l)
    # Transpose the data so each row is a species
    transf_data <- as.data.frame(t(transf_data))
    # Remove any gene families that are invariant
    keep <-
      which(unlist(lapply(apply(
        transf_data, 2, unique
      ), length)) > 0)
    transf_data <- transf_data[, keep]
    # Add response status
    transf_data$response <-
      as.factor(as.numeric(as.factor(response)) - 1)
    # Add a sample id column:
    transf_data$sample_id <- rownames(transf_data)
    # Conduct lasso binomial regression as implemented in the
    # familiar package which optimized hyperparameters and
    # feature selection
    dir.create(
      paste0("trait_prediction/", response_name, "_lasso/"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    mod_out_dir <-
      file.path(paste0(
        "trait_prediction/",
        response_name,
        "_lasso/",
        predictor_name
      ))
    # Fit the model:
    summon_familiar(
      data = transf_data,
      experiment_dir = mod_out_dir,
      outcome_type = "binomial",
      outcome_column = "response",
      sample_id_column = "sample_id",
      experimental_design = "cv(fs+mb, 10, 1)",
      fs_method = "lasso_binomial",
      learner = "lasso",
      normalisation_method = "quantile",
      cluster_method = cluster_method,
      cluster_similarity_threshold = "0.5",
      parallel = TRUE,
      parallel_nr_cores = 10
    )
    print("Finished fitting model!")
    # Set the path to the output
    model_directory_path <-
      file.path(mod_out_dir, "trained_models", "lasso", "lasso_binomial")
    # And set a path directly to the model itself
    model_paths <-
      file.path(model_directory_path,
                list.files(model_directory_path, pattern = "model"))
    # Read in fitted models
    print("Reading in ensemble and model collections")
    ensemble <- as_familiar_ensemble(model_paths)
    collection_path <-
      list.files(paste0(mod_out_dir, "/familiar_collections"),
                 full.names = TRUE)
    collection <- readRDS(collection_path)
    # Pull out variable importance data:
    vimps <- get_vimp_table(model_paths)
    ncol_vimps <- c()
    for (i in 1:10) {
      if (!is.null(vimps[[i]])) {
        vimps[[i]]$model_id <- i
        ncol_vimps <- c(ncol_vimps, ncol(vimps[[i]]))
      } else {
        ncol_vimps <- c(ncol_vimps, 0)
      }
    }
    vimps <-
      as.data.frame(do.call(rbind, vimps[which(ncol_vimps == 4)]))
    # Now, read in the models so we can get the coefficients of each
    # feature for each model
    models <- list()
    for (i in 1:10) {
      models[[i]] <- readRDS(model_paths[[i]])
      tmp <- coef(models[[i]])
      if (!is.null(tmp)) {
        if (i == 1) {
          if (is.numeric(tmp)) {
            tmp <-
              data.frame(
                set = names(rev(tmp)),
                s1 = rev(tmp),
                row.names = NULL
              )
          } else {
            tmp <- data.frame(set = rownames(tmp), as.matrix(tmp))
          }
          coefs <- as.data.frame(as.matrix(tmp))
          coefs$model_id <- i
        } else {
          if (is.numeric(tmp)) {
            tmp <-
              data.frame(
                set = names(rev(tmp)),
                s1 = rev(tmp),
                row.names = NULL
              )
          } else {
            tmp <- data.frame(set = rownames(tmp), as.matrix(tmp))
          }
          tmp$model_id <- i
          coefs <- rbind(coefs, tmp)
        }
      }
    }
    # Write variable importance and coefs out to file:
    write.table(
      vimps,
      paste0(mod_out_dir, "/feature_importance_cv.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    write.table(
      coefs,
      paste0(mod_out_dir, "/feature_coefficients_cv.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    # Print out a bunch of diagnostic plots
    print("Plotting model performance")
    mod_perform <- plot_model_performance(
      object = ensemble,
      draw = TRUE,
      facet_by = "metric",
      detail_level = "ensemble",
      data = transf_data,
      metric = c("auc", "accuracy", "brier", "f1_score")
    )
    ggsave(
      mod_perform[[1]],
      file = paste0(mod_out_dir, "/model_performance.pdf"),
      height = 7,
      width = 7
    )
    print("Plotting variable importance")
    mod_var_import <-
      plot_variable_importance(
        object = collection,
        type = "feature_selection",
        rotate_x_tick_labels = TRUE,
        draw = TRUE,
        facet_by = "metric",
        data = transf_data
      )
    ggsave(
      mod_var_import[[1]],
      file = paste0(mod_out_dir, "/model_variable_importance.pdf"),
      height = 7,
      width = 7
    )
    print("Plotting univariate feature importance")
    univar_import <-
      plot_univariate_importance(collection,
                                 draw = TRUE)
    ggsave(
      univar_import[[1]],
      file = paste0(mod_out_dir, "/univariate_variable_importance.pdf"),
      height = 7,
      width = 7
    )
    print("Plotting permutation variable importance")
    perm_var_importance <-
      plot_permutation_variable_importance(
        object = ensemble,
        detail_level = "ensemble",
        data = transf_data,
        estimation_type = "bias_correction",
        rotate_x_tick_labels = TRUE,
        draw = TRUE
      )
    ggsave(
      perm_var_importance[[1]],
      file = paste0(mod_out_dir, "/permutation_variable_importance.pdf"),
      height = 7,
      width = 7
    )
    print("Plotting confusion matrix")
    conf_mat <-
      plot_confusion_matrix(object = collection,
                            data = transf_data,
                            draw = TRUE)
    ggsave(
      conf_mat[[1]],
      file = paste0(mod_out_dir, "/confusion_matrix.pdf"),
      height = 7,
      width = 7
    )
    print("Plotting feature similarity")
    feat_simil <-
      plot_feature_similarity(
        object = collection,
        data = transf_data,
        feature_similarity_metric = "spearman",
        rotate_x_tick_labels = TRUE,
        draw = TRUE
      )
    pdf(
      paste0(mod_out_dir, "/feature_similarity.pdf"),
      height = 7,
      width = 7
    )
    grid::grid.draw(feat_simil[[1]])
    dev.off()
    print("Plotting sample clustering")
    pdf(
      paste0(mod_out_dir, "/sample_feature_clustering.pdf"),
      height = 7,
      width = 7
    )
    grid::grid.draw(plot_sample_clustering(collection,
                                           rotate_x_tick_labels = TRUE)[[1]])
    dev.off()

    return(
      list(
        models = models,
        collection = collection,
        ensemble = ensemble,
        var_importance = vimps,
        coefs = coefs
      )
    )
  }
