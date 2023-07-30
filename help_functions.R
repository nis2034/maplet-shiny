
# Extract the object names from result SE "D"------------
get_obj_name <- function(D){
  obj_list <- data.frame()
  for (i in seq_along(metadata(D)$results)) {
    for (j in seq_along(metadata(D)$results[[i]]$fun)) {
      obj_list[i, j] <- metadata(D)$results[[i]]$fun[j]
    }
  }
  
  #############################
  
  ###################################
  
  
  # Extract the stat_name
  stat_name <- data.frame(stat_name=NA)
  for (i in seq_along(metadata(D)$results)) {
    stat_name[i, 1] <- if ('stat_name' %in% names(metadata(D)$results[[i]]$args)) {
      metadata(D)$results[[i]]$args$stat_name
    } else if ('stat_list' %in% names(metadata(D)$results[[i]]$args) & class(metadata(D)$results[[i]]$args$stat_list)!="call") {
      metadata(D)$results[[i]]$args$stat_list
    } else {
      NA
    }
  }
  
  # Merge object names and stat_name
  order_id <- 1:nrow(obj_list)
  obj_name <- cbind(order_id, obj_list, stat_name) 
  obj_name$stat_name <- 
    ifelse(is.na(obj_name$stat_name), 
           "(no stat_name)", 
           obj_name$stat_name)
  
  # count number of objects in the list of output of each metadata
  count_list <- c()
  for (i in seq_along(metadata(D)$results)){
    count_list[i] <- length(metadata(D)$results[[i]]$output)
  }
  obj_name$cnt <- count_list
  return(obj_name)
}

# get pathway annotations
get_pathway_annotations <- function(D, pwvar) {
  rd <- rowData(D) %>% as.data.frame %>%
    dplyr::mutate(name = rownames(rowData(D))) %>%
    dplyr::select(name,!!sym(pwvar), BIOCHEMICAL)
  if (class(rd[[pwvar]][1]) == "AsIs") {
    # remove rows with missing pathway annotations and unnest pathway column
    miss_idx <-
      apply(rd, 1, function(x) {
        x[[pwvar]][[1]] %>% is.null() %>% unname()
      })
    rd <- rd[!miss_idx, ] %>% tidyr::unnest(!!sym(pwvar))
    # extract pathway_name column form pathways data frame
    rd %<>% dplyr::left_join(metadata(D)$pathways[[pwvar]], by = setNames("ID", pwvar)) %>%
      dplyr::select(name, BIOCHEMICAL, pathway_name)
    # replace value for pathway column variable
    pwvar <- "pathway_name"
  }
  
  rd
}

# Get the data set -----------------------
get_data_by_name <- function(D, args.name, plot.nmae, stat.name) {
  plots <- mtm_res_get_entries(D, c("plots", plot.nmae))
  for (i in seq_along(plots)) {
    if (plots[[i]]$args[args.name] == stat.name) {
      data <- plots[[i]]$output[[1]]$data
    }
  }
  data
}


# create reverselog_trans for log10-SCALE in volcano plot-----------------
reverselog_trans <- function (base = exp(1)){
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf))
}

# get the threshold for significance (extracted from corresponding stat_bar plot)
get_threshold_for_p_adj <- function(D, stat_name) {
  # define threshold for significance (extracted from corresponding stat_bar plot)
  stats_plots <- mtm_res_get_entries(D, c("plots", "stats"))
  for (plot in stats_plots) {
    if (plot$args$stat_list == stat_name) {
      alpha <- plot$args$feat_filter[3][[1]]
    }
  }
  alpha
}


# define PCA output function for mod3 referring 'mt_plots_pca'------------
mod3_plots_pca <- function(D, title = "PCA", 
                           ## scale argument
                           scale_data,
                           ## color argument
                           color,
                           categorizing,
                           pc1 = 1, pc2 = 2, 
                           data_type,
                           hover,
                           ...){
  X = t(assay(D))
  if (any(is.na(X))) 
    stop("Data matrix for PCA cannot contain NAs")
  if (!(data_type %in% c("scores", "loadings"))) 
    stop("Show must be either 'scores' or 'loadings'")
  # scale scores if scale checkbox=T
  if (scale_data) 
    X <- scale(X)
  pca <- stats::prcomp(x = as.matrix(X), center = F, scale = F)
  expvar <- (pca$sdev)^2/sum(pca$sdev^2)
  ## reactivate axis labels
  pc1name <- sprintf("PC%d (%.1f%%)", pc1, expvar[pc1] * 
                       100)
  pc2name <- sprintf("PC%d (%.1f%%)", pc2, expvar[pc2] * 
                       100)
  
  if (data_type == "scores") {
    df <- data.frame(x = pca$x[, pc1], 
                     y = pca$x[, pc2], 
                     colData(D)
    )
    ## reactivate plot title
    plot_title <- paste0(ifelse(scale_data,
                                "Scaled ", 
                                "Non-scaled "),
                         title, " - ", color)
    ## categorize coloring if color checkbox=T
    if(categorizing){
      df[, color] <- factor(df[, color])
    }
    # customize hover text
    # hover_text <- paste0(names(data.frame(colData(D)))[as.numeric(hover)], ": ", 
    #                     data.frame(colData(D))[[as.numeric(hover)]])
    # draw ggplot
    p <- ggplot(data = df, 
                do.call(aes_string, as.list(structure(c("x","y",color,hover), names = c("x","y","colour",hover))))
    ) + 
      geom_point() + 
      xlab(pc1name) + 
      ylab(pc2name) + 
      ggtitle(plot_title) +
      ## reactive legend title
      labs(color = color)
  } else {
    df = data.frame(x = pca$rotation[, pc1], 
                    y = pca$rotation[, pc2], 
                    rowData(D)
    )
    ## reactivate plot title
    plot_title <- paste0(ifelse(scale_data, 
                                "Scaled ", 
                                "Non-scaled "),
                         title)
    ## categorize coloring if color checkbox=T
    if(categorizing){
      df[, color] <- factor(df[, color])
    }
    # customize hover text
    # hover_text <- paste0(names(data.frame(rowData(D)))[as.numeric(hover)], ": ", 
    #                      data.frame(rowData(D))[[as.numeric(hover)]])
    # draw ggplot
    p <- ggplot(data = df, 
                do.call(aes_string, as.list(structure(c("x","y",color,hover), names = c("x","y","colour",hover))))
    ) + 
      geom_point() + 
      xlab(pc1name) + 
      ylab(pc2name) + 
      ggtitle(plot_title)
  }
  
  # draw plotly
  ggplotly(p, tooltip = c(color, hover)) %>%
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = 0.5, ## set position of legend
                         y = -0.2,
                         tracegroupgap = 5),
           autosize = TRUE
    )
}


# define UMAP output function for mod3 referring 'mt_plots_umap'---------------
mod3_plots_umap <- function (D, title = "UMAP", 
                             ## scale argument
                             scale_data,  
                             ## color argument
                             color,
                             categorizing,
                             n_neighbors, 
                             hover,
                             ...) {
  X <- t(assay(D))
  if (any(is.na(X))) 
    stop("Data matrix for UMAP cannot contain NAs")
  ## dependent on scale checkbox
  if (scale_data) 
    X <- scale(X)
  um <- umap::umap(d = as.matrix(X), n_neighbors = as.numeric(n_neighbors))
  df <- data.frame(x = um$layout[, 1], y = um$layout[, 2], colData(D))
  ## reactivate plot title
  plot_title <- paste0(ifelse(scale_data,
                              "Scaled ", 
                              "Non-scaled "), title, 
                       " with ", n_neighbors, 
                       " neighbors colored by ", color)
  ## categorize coloring if color checkbox=T
  if(categorizing){
    df[, color] <- factor(df[, color])
  }
  # customize hover text
  # hover_text <- paste0(names(data.frame(colData(D)))[as.numeric(hover)], ": ", 
  #                      data.frame(colData(D))[[as.numeric(hover)]])
  # draw ggplot
  p <- ggplot(data = df,
              do.call(aes_string, as.list(structure(c("x","y",color,hover), names = c("x","y","color","hover"))))
  ) + 
    geom_point() + 
    xlab("comp 1") + 
    ylab("comp 2") + 
    ggtitle(plot_title) +
    ## reactive legend title
    labs(color = color)
  
  # draw plotly
  ggplotly(p, tooltip = c(color, hover)) %>% 
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = .5,
                         y = -.2,
                         tracegroupgap = 5),
           autosize = TRUE
    )
  
}

# define PLS output function for mod3 referring 'mt_plots_pls'---------------
# mod3_plots_pls <- function (D, title = "PLS",
#                             subgroupvar,hover, ...) {
# 
# 
# 
# 
# 
#   p <- NULL
# 
#   Dplot <- D
#   # subgroupvar <- "GROUP_ID"
# 
# 
#   pls_D <- autonomics:::pls(Dplot, subgroupvar,color =!!ensym(subgroupvar))
# 
#   p <- autonomics:::biplot(pls_D ,pls1, pls2,color =!!ensym(subgroupvar) )
# 
#   # draw plotly
#   ggplotly(p, tooltip = c(subgroupvar, hover)) %>%
#     layout(legend = list(orientation = "h",   # show entries horizontally
#                          xanchor = "center",  # use center of legend as anchor
#                          x = .5,
#                          y = -.2,
#                          tracegroupgap = 5),
#            autosize = TRUE
#     )
# 
# }


library(ggplot2)
library(plotly)

mod3_plots_pls <- function(D, title = "PLS", subgroupvar, hover, ...) {
  
  p <- NULL
  Dplot <- D
  
  # Perform partial least squares analysis
  pls_D <- autonomics:::pls(Dplot, subgroupvar)
  
  print("done")
  
  if (identical(D, pls_D))
  {print("no change")}
  else{
    print("pls done")
    print(str(pls_D))
  }
  
  # Extract the scores for the x and y axes from the PLS result
  p <- autonomics:::biplot(pls_D ,pls1, pls2,color =!!ensym(subgroupvar) )

    # draw plotly
    ggplotly(p, tooltip = c(subgroupvar, hover)) %>%
      layout(legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = .5,
                           y = -.2,
                           tracegroupgap = 5),
             autosize = TRUE
      )
  
}






# define boxplot function in module 5---------------
mod5_boxplot <- function(D, x, x_cate, y, y_cate, fill, hover, ...){
  df <- data.frame(colData(D))
  ## categorize variable if user think it's not continuous
  if(x_cate==FALSE){
    df[, x] <- factor(df[, x])
  }
  if(y_cate==FALSE){
    df[, y] <- factor(df[, y])
  }
  # to reproduce jitter plot
  set.seed(4017)
  p <- ggplot(df,
              do.call(aes_string, as.list(structure(c(x, y, fill, hover), names = c("x","y", "fill", "hover"))))
  ) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .2, alpha = 0.5) +
    theme(legend.title = element_blank()) +
    ggtitle("Boxplot with ignored outliers and Jitter with Set Seed")
  
  # get number of boxplots from plot object
  nbox <- p$data %>% dplyr::pull(p$mapping$x[[2]]) %>% unique %>% length
  
  # draw plotly
  fig <- ggplotly(p) %>%
    layout(legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3))
  
  # remove outlier dots from boxplot layer (but not from jitter layer)
  lapply(1:nbox, function(i){
    fig$x$data[i] <<- lapply(fig$x$data[i], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
  })
  # return fig
  fig
  
}


# define scatterplot function in mod5------------------
mod5_scatter <- function(D, x, y, hover, ...){
  df <- data.frame(colData(D))
  p <- ggplot(df,
              do.call(aes_string, as.list(structure(c(x, y, hover), names = c("x","y", "hover"))))
  ) + geom_point()
  
  ggplotly(p)
  

}


# define barplot function in mod5--------------------
mod5_barplot <- function(D, x, fill, hover, ...){
  df <- data.frame(colData(D))
  p <- ggplot(df,
              do.call(aes_string, as.list(structure(c(x, fill, hover), names = c("x","fill", "hover"))))
  ) + geom_bar()
  
  ggplotly(p) %>%
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = 0.5, ## set position of legend
                         y = -0.2,
                         tracegroupgap = 5),
           autosize = TRUE
    )
  
}

# SE generating function for Mod6
diff_analysis_func <- function(D,
                               var,
                               binary=F,
                               analysis_type="lm",
                               mult_test_method="BH",
                               alpha=0.05,
                               group_col_barplot,
                               color_col_barplot=NULL){
  D %<>%
    mt_reporting_heading(heading = sprintf("%s Differential Analysis",var), lvl = 2) %>%
    {.}
  
  if(analysis_type=="lm"){
    D %<>%
      mt_stats_univ_lm(formula = as.formula(sprintf("~  %s",var)), stat_name = sprintf("%s analysis",var)) %>%
      {.}
  }else{
    D %<>%
      mt_stats_univ_cor(in_col = var, stat_name = sprintf("%s analysis",var),method = analysis_type) %>%
      {.}
  }
  
  if(binary){
    D %<>%
      mt_post_fold_change(stat_name = sprintf("%s analysis",var))
  }
  D %<>%
    mt_post_multtest(stat_name = sprintf("%s analysis",var), method = mult_test_method) %>%
    mt_reporting_stats(stat_name = sprintf("%s analysis",var), stat_filter = p.adj < alpha) %>%
    mt_plots_volcano(stat_name = sprintf("%s analysis",var),
                     x = !!sym(ifelse(binary,"fc","statistic")),
                     feat_filter = p.adj < alpha,
                     color = p.adj < alpha) %>%
    mt_plots_box_scatter(stat_name = sprintf("%s analysis",var),
                         x = !!sym(var),
                         plot_type = ifelse(binary,"box","scatter"),
                         feat_filter = p.adj < alpha,
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
    mt_plots_stats_pathway_bar(stat_list = sprintf("%s analysis",var),
                               y_scale = "count",
                               feat_filter = p.adj < alpha,
                               group_col = group_col_barplot,
                               color_col = color_col_barplot) %>%
    {.}
  
  return(D)
}


# SE generating function for Mod7
diff_analysis_func_tab <- function(D,
                                   var,
                                   binary=F,
                                   sample_filter = NULL,
                                   filter_val =NULL,
                                   covar_col_select = NULL,
                                   covar_row_select = NULL,
                                   analysis_type="lm",
                                   mult_test_method="BH",
                                   alpha=0.05,
                                   group_col_barplot,
                                   color_col_barplot=NULL){
  
  
  
  print("samp fileter")
  print(sample_filter)
  
  print(filter_val)
  if(is.null(covar_col_select) && is.null(covar_row_select))
  {
    covar = NULL
  } else
  {
    covar <- paste("+", paste(c(covar_row_select,covar_col_select), collapse = "+"), sep = "")
  }
  D %<>%
    mt_reporting_heading(heading = sprintf("%s Differential Analysis",var), lvl = 2) %>%
    {.}
  
  if(analysis_type=="lm" ){ 
    if (!is.null(sample_filter))
    {
      
      D %<>%
        mt_stats_univ_lm(formula = as.formula(sprintf("~  %s%s",var, replace(covar, is.null(covar),""))), stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                         samp_filter = (!!sym(sample_filter) %in% filter_val)) %>%
        
        {.}
    }else{
      D %<>%
        mt_stats_univ_lm(formula = as.formula(sprintf("~  %s%s",var, replace(covar, is.null(covar),""))), stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),""))) %>%
        
        {.}
    }
  }
  
  else{
    D %<>%
      mt_stats_univ_cor(in_col = var, stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),method = analysis_type) %>%
      {.}
  }
  
  if(binary){
    print("binary analysis")
    D %<>%
      mt_post_fold_change(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")))
  }
  D %<>%
    mt_post_multtest(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")), method = mult_test_method) %>%
    mt_reporting_stats(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")), stat_filter = p.adj < alpha) %>%
    mt_plots_volcano(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                     x = !!sym(ifelse(binary,"fc","statistic")),
                     feat_filter = p.adj < alpha,
                     color = p.adj < alpha) %>%
    mt_plots_box_scatter(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                         x = !!sym(var),
                         plot_type = ifelse(binary,"box","scatter"),
                         feat_filter = p.adj < alpha,
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
    mt_plots_stats_pathway_bar(stat_list = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                               y_scale = "count",
                               feat_filter = p.adj < alpha,
                               group_col = group_col_barplot,
                               color_col = color_col_barplot) %>%
    {.}

  return(D)
}



# SE generating function for for DE part of mod7
diff_analysis_func_tab_new <- function(D,
                                       var,
                                       binary=F,
                                       sample_filter = NULL,
                                       filter_val =NULL,
                                       covar_col_select = NULL,
                                       covar_row_select = NULL,
                                       analysis_type="lm",
                                       mult_test_method="BH",
                                       alpha=0.05){
  
  
  
  print("samp fileter")
  print(sample_filter)
  
  print(filter_val)
  if(is.null(covar_col_select) && is.null(covar_row_select))
  {
    covar = NULL
  } else
  {
    covar <- paste("+", paste(c(covar_row_select,covar_col_select), collapse = "+"), sep = "")
  }
  D %<>%
    mt_reporting_heading(heading = sprintf("%s Differential Analysis",var), lvl = 2) %>%
    {.}
  
  if(analysis_type=="lm" ){ 
    if (!is.null(sample_filter))
    {
      
      D %<>%
        mt_stats_univ_lm(formula = as.formula(sprintf("~  %s%s",var, replace(covar, is.null(covar),""))), stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                         samp_filter = (!!sym(sample_filter) %in% filter_val)) %>%
        
        {.}
    }else{
      D %<>%
        mt_stats_univ_lm(formula = as.formula(sprintf("~  %s%s",var, replace(covar, is.null(covar),""))), stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),""))) %>%
        
        {.}
    }
  }
  
  else{
    D %<>%
      mt_stats_univ_cor(in_col = var, stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),method = analysis_type) %>%
      {.}
  }
  
  if(binary){
    D %<>%
      mt_post_fold_change(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")))
  }
  D %<>%
    mt_post_multtest(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")), method = mult_test_method) %>%
    mt_reporting_stats(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")), stat_filter = p.adj < alpha) %>%
    mt_plots_volcano(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                     x = !!sym(ifelse(binary,"fc","statistic")),
                     feat_filter = p.adj < alpha,
                     color = p.adj < alpha) %>%
    mt_plots_box_scatter(stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                         x = !!sym(var),
                         plot_type = ifelse(binary,"box","scatter"),
                         feat_filter = p.adj < alpha,
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
    {.}
  
  return(D)
}




# get the volcano plot in mod2--------------------
mod2_plot_vol <- function(D, inputs, legend_name, d, pwvar, alpha) {
  # Get volcano data set if the plot1 is null
  data_vol <- get_data_by_name(D, "stat_name", "volcano", inputs[1])
  
  # Get volcano data set if the plot1 is bar
  if (inputs[2] == "bar") {
    if (!is.null(d)) {
      # get the click information for the bar plot 
      data_bar <- get_data_by_name(D, "stat_list", "stats", inputs[1])
      lvls <- rev(levels(data_bar$label))
      label <- lvls[round(as.numeric(d$y))]
      sub_pathway_name <- data_bar[data_bar$label == label, ]$name
      # get the variable name by sub_pathway_name
      row_data <- rowData(D) %>% data.frame()
      names <- unlist(row_data[row_data[pwvar] == sub_pathway_name,]$name)
      data_vol <- data_vol[data_vol$name %in% names, ]
    }
  }
  # Set the legend color column
  data_vol[, legend_name] <- ifelse(data_vol$p.adj < alpha, TRUE, FALSE)
  
  # draw the plot2
  plot <- data_vol %>%
    ggplot(aes(x = statistic, y = p.value, color = !!sym(legend_name), label = name)) +
    geom_point() +
    scale_y_continuous(trans = reverselog_trans(10),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(y = "p-value (10^(-y))") +
    ggtitle(paste0(inputs[1])) +
    scale_color_brewer(palette="Dark2")
  
  if (inputs[2] == "bar") {
    plot <- plot +
      geom_point(size = 3) +
      ggtitle(paste0(sub_pathway_name, "-", inputs[1]))
  }
  
  plot
}



# get the equalizer plot in mod2--------------------
mod2_plot_eq <- function(D, inputs, rd, alpha, pwvar, path_name, d) {
  # add pathway annotations to results
  res <- maplet::mtm_get_stat_by_name(D, inputs[1]) %>%
    dplyr::left_join(rd, by=c("var"="name")) %>%
    dplyr::mutate(x = sign(statistic)*log10(p.adj)) %>%
    dplyr::filter(!is.na(BIOCHEMICAL))
  # compute multiple testing correction line
  sel <- res %>% 
    dplyr::filter(p.adj < alpha)
  if(nrow(sel)>0){
    xfine <- sel %>% .$p.adj %>% max(., na.rm = T)
  } else {
    xfine <- Inf
  }
  
  # If plot1 is not bar
  if (inputs[2] == "null") {
    df <- res %>%
      dplyr::filter(!!sym(pwvar)==path_name)
  } else {
    if (!is.null(d)) {
      # get the click information for the bar plot
      data_bar <- get_data_by_name(D, "stat_list", "stats", inputs[1])
      lvls <- rev(levels(data_bar$label))
      label <- lvls[round(as.numeric(d$y))]
      sub_pathway_name <- data_bar[data_bar$label == label,]$name
      
      # create equalizer plots
      df <- res %>%
        dplyr::filter(!!dplyr::sym(pwvar) == sub_pathway_name)
    }
  }
  
  # colors
  clrs <- c("#9494FF","red")
  # x axis limits
  a = max(abs(df$x)) + 0.3
  plot <-
    ggplot(df, aes(x = x, y = BIOCHEMICAL)) +
    geom_vline(xintercept = 0, color = "gray") +
    (if (!is.infinite(xfine)) {
      geom_vline(
        xintercept = c(-log10(xfine), log10(xfine)),
        color = "red",
        alpha = 0.4
      )
    }) +
    geom_point(pch = 22,
               fill = clrs[1],
               size = 3) +
    facet_grid(as.formula(sprintf("%s~.", pwvar)),
               scales = "free_y",
               space = "free_y") +
    theme(
      strip.background = element_rect(fill = NA),
      strip.text = element_text(colour = 'black', face = "bold"),
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_rect(fill = NA, color =
                                        "black")
    ) +
    ylab("") +
    xlab("sign(statistic)*log10(p.adj)") +
    scale_x_continuous(limits = c(-a, a))
  
  if (inputs[2] == "bar") {
    plot <- plot +
      (if (!is.infinite(xfine)) {
        ggtitle(paste0(sub_pathway_name, "-", inputs[1],"-",
                       sprintf("Differential Metabolites at alpha %.2f", alpha)))
      } else{
        ggtitle(sprintf("No significant results at alpha %.2f", alpha))
      })
  }
  
  plot
}


# get the box/scatter plot in mod2--------------------
mod2_plot_box_scatter <- function(D,
                                  inputs,
                                  d.bar,
                                  d.eq,
                                  d.vol,
                                  rd,
                                  pwvar,
                                  path_name,
                                  alpha,
                                  is_categorical,
                                  data) {
  
  # get the bar name if plot1 is selected "bar"
  if (inputs[2] == "bar") {
    if (!is.null(d.bar)) {
      data_bar <- get_data_by_name(D, "stat_list", "stats", inputs[1])
      lvls <- rev(levels(data_bar$label))
      label <- lvls[round(as.numeric(d.bar$y))]
      sub_pathway_name <- data_bar[data_bar$label == label,]$name
    }
  }
  
  # Get the metabolite name by click information
  if (inputs[3] == "equalizer") {
    if (!is.null(d.eq)) {
      res <- maplet::mtm_get_stat_by_name(D, inputs[1]) %>%
        dplyr::left_join(rd, by=c("var"="name")) %>%
        dplyr::mutate(x = sign(statistic)*log10(p.adj)) %>%
        dplyr::filter(!is.na(BIOCHEMICAL))
      
      # get the data set of the equalizer plot
      if (inputs[2] == "bar") {
        # Filter the data by selected pwvar
        df <- res %>%
          dplyr::filter(!!dplyr::sym(pwvar) == sub_pathway_name)
      } else {
        df <- res %>%
          dplyr::filter(!!dplyr::sym(pwvar) == path_name)
      }
      # get the metabolite name by click info for equalizer plot
      metabolite <- df[as.numeric(d.eq$pointNumber) + 1,]$var[1]
      term <- df$term[1]
    }
  } else { # get the click info in volcano plot
    if (!is.null(d.vol)) {
      data_vol <- get_data_by_name(D, "stat_name", "volcano", inputs[1])
      if (inputs[2] == "bar") {
        # get the variable name by sub_pathway_name
        row_data <- rowData(D) %>% data.frame()
        names <- unlist(row_data[row_data[pwvar] == sub_pathway_name,]$name)
        data_vol <- data_vol[data_vol$name %in% names, ]
      }
      # set the column curveNumber by color legend
      p.adj.significant <- alpha
      data_vol[, "curveNumber"] <- ifelse(data_vol$p.adj < as.numeric(p.adj.significant), 1, 0)
      data_vol_true <- data_vol[data_vol$curveNumber==1, ]
      data_vol_false <- data_vol[data_vol$curveNumber==0, ]
      # By using click info (curveNumber & ponitNumber) to get the metabolite name
      metabolite <- ifelse(d.vol$curveNumber == 1,
                           data_vol_true[d.vol$pointNumber + 1, ]$var[1],
                           data_vol_false[d.vol$pointNumber + 1, ]$var[1])
      term <- data_vol$term[1]
      
    }
  }
  
  # Filter the data by metabolite name
  data <- data[data$var == metabolite, ]
  
  # Treat as categorical or not?
  if (is_categorical) {
    data[, term] <- factor(data[, term])
  } else {
    data[, term] <- as.numeric(data[, term])
  }
  
  # Draw the plot3
  if (inputs[4] == "scatter") {
    plot <- data %>%
      ggplot(aes(x = !!sym(term), y = value)) +
      geom_point(size = 3) +
      geom_smooth(method = "lm", se = T, color = "black") + 
      ggtitle(metabolite)
  } else {
    plot <- data %>%
      ggplot(aes(x = !!sym(term), y = value)) +
      geom_boxplot() +
      geom_jitter(size = 3, width = 0.2) +
      ggtitle(metabolite)
  }
  
  plot
}


# get log text from a SE object
get_log_text <- function(D){
  dt <- metadata(D)$results
  text <- lapply(names(dt), function(x){
    sprintf("%s\n",dt[[x]]$logtxt)
  }) %>% {do.call(cat,.)}
  # return text
  text
}


# get box plots from a SE object
get_plots_SE <- function(D) {
  dt <- metadata(D)$results
  plots <- list()
  n <- length(names(dt))
  j <- 1
  
  for (i in 1:n) {
    if (dt[[i]]$fun[1] == "plots") { #and fun[2] = pathway for pathway enrichment function
      plots[[j]] <- dt[[i]]$output
      j <- j + 1
    }
  }
  
  plots
}





get_plots_SE_differ <- function(D){
  dt <- metadata(D)$results
  plots <- list()
  n= length(names(dt))
  j=1
  for(i in 1:n){
    if(dt[[i]]$fun[1] =="plots")
    {
      
      if((replace(dt[[i]]$fun[2], is.na(dt[[i]]$fun[2]),"") =="volcano") || (replace(dt[[i]]$fun[3],is.na(dt[[i]]$fun[3]),"") =="scatter")||(replace(dt[[i]]$fun[3], is.na(dt[[i]]$fun[3]),"") =="pathway") || (replace(dt[[i]]$fun[2], is.na(dt[[i]]$fun[2]),"") =="net"))
      {
        print("nets")
        plots[[j]] <- dt[[i]]$output
        
        j = j+1
        
      }
    } 
  } 
  plots
}

get_plots_SE_preprocess <- function(D, title){
  dt <- metadata(D)$results
  plots <- list()
  n= length(names(dt))
  j=1
  for(i in 1:n){
    if(dt[[i]]$fun[1] =="plots")
    { 
      
      if (replace(dt[[i]]$args$title, is.null(dt[[i]]$args$title),"")  %in% title)
        
        
      {
        plots[[j]] <- dt[[i]]$output
        
        j = j+1
        
      }
    } 
  } 
  plots
}

# dt[[1]]$args$title
# 
# dt <- metadata(D1)$results
# 
# dt[[12]]$fun[1]
# title = c("Original","After imputation","After filtering 50% missing features")
# 
#  if (replace(dt[[12]]$args$title, is.null(dt[[i]]$args$title),"") %in% title)
#   { print("yes")
#  }else {
#   print("no")}



mt_anno_pathways_hmdb_new <- function(D,
                                      in_col,
                                      out_col,
                                      pwdb_name = "SMP",
                                      db_dir = system.file("extdata", "precalc/hmdb/", package = "maplet") ,
                                      db_file,
                                      raw_db_outfile) {
  
  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  if (!dir.exists(db_dir))
    stop(glue::glue("{db_dir} does not exist. input a valid db_dir."))
  
  if (missing(db_file)) {
    # NOTE: taking the tail work only for HMDB versions below 9, since
    # list.files sorts the names alphabetically. Current HMDB version
    # is 4.0
    db_file <- list.files(db_dir,pattern='*.rds') %>% utils::tail(1)
    if (length(db_file) == 0)
      stop(glue::glue("no pathway files were found in {db_dir}"))
  }
  
  if (!pwdb_name %in% c("SMP", "KEGG"))
    stop(glue::glue("pwdb_name is invalid. please use either SMP or KEGG for pwdb_name."))
  
  if(!in_col %in% names(rowData(D)))
    stop(glue::glue("in_col is invalid. please enter an existing column name."))
  
  # HMDB IDs here pertain to secondary accessions
  #mti_logstatus(glue::glue("reading the {pwdb_name} database from {db_file}"))
  pwdb <-
    readr::read_rds(file.path(db_dir, db_file)) %>%
    dplyr::select(HMDB_id, ID = pwdb_name, pathway_name, accession)
  
  # create a dataframe that enables the mapping between pathway
  # names and IDs. Included also are num_total, num_measured,
  # num_pw_total, and num_pw_measured (see below for further details)
  
  # num_total - overall number of metabolites in entire database (this will
  # be a redundant, repeating number, identical in every row… but it’s the
  # easiest way to store it right now)
  # pwdb$accession is used here, since there is a many-to-one mapping
  # between HMDB_id (secondary accessions) and accession
  num_total <-
    pwdb$accession %>%
    unique() %>%
    length()
  
  # num_measured - overall number of measured metabolites found in the DB
  num_measured =
    pwdb$HMDB_id %>%
    # for this intersect, similar to the one in pwdb_summary (below),
    # it is assumed that the HMDB IDs in the dataset D is non-redundant.
    # else this number may not be accurate.
    intersect(rowData(D)[[in_col]]) %>%
    length()
  
  #mti_logstatus(glue::glue("summarizing {pwdb_name} pathway information"))
  pwdb_summary <- pwdb
  # using methods from data.table reduces runtime by almost 10x as compared
  # to dplyr
  pwdb_summary <-
    data.table::setDT(pwdb_summary)[, `:=`(
      
      num_total = num_total,
      num_measured = num_measured,
      
      # num_pw_total - the total number of metabolites in that pathway
      # here, accession is used.
      # (overall DB background)
      num_pw_total = data.table::uniqueN(accession),
      
      # num_pw_measured - the number of measured metabolites in that pathway
      # (for the M type of analysis on the actual measured background).
      num_pw_measured =
        sum(HMDB_id %in% rowData(D)[[in_col]], na.rm = TRUE)
    ),
    by = ID] %>%
    # some IDs might have more than two names, however, these will be discarded
    # for now
    unique(by = c("ID")) %>%
    subset(!is.na(ID),
           select = -c(HMDB_id, accession))
  
  
  #mti_logstatus(glue::glue("nesting {pwdb_name} pathway IDs"))
  # nest all the pathway IDs given our lieblings input IDs (HMDB IDs for now)
  pwdb_reduced <-
    pwdb %>%
    dplyr::group_by(HMDB_id) %>%
    dplyr::filter(!is.na(ID)) %>%
    dplyr::distinct(HMDB_id, ID) %>%
    #tidyr::nest(ID, .key = ID) %>%
    tidyr::nest(.key = "ID") %>%
    dplyr::mutate(ID =
                    ID %>%
                    unlist(recursive = FALSE) %>%
                    as.list())
  
  
  
  # match the nested pathways to our lieblings IDs
  match_idx <-
    match(rowData(D)[[in_col]],
          pwdb_reduced$HMDB_id)
  
  pw_col <- pwdb_reduced$ID[match_idx]
  
  # add the pathway IDs into D
  rowData(D)[[out_col]] <- pw_col
  
  # add pathway map to the metadata of D
  metadata(D)$pathways[[out_col]] <-
    pwdb_summary
  
  
  if(!missing(raw_db_outfile)) {
    openxlsx::write.xlsx(pwdb, raw_db_outfile)
  }
  
  
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('added pathway annotations using the %s pathway database', pwdb_name)
    )
  
  D
}


mt_stats_pathway_enrichment_new <- function(D, stat_name, pw_col, cutoff = 0.05) {
  
  stopifnot("SummarizedExperiment" %in% class(D))
  
  meta_D <- metadata(D)
  
  if(!"pathways" %in% names(meta_D)) stop("'pathways' does not exist in current SummarizedExperiment input")
  
  # Check if given pathway column actually exists
  if (!pw_col %in% names(meta_D$pathways)) stop(sprintf("'%s' not found in metabolite annotations.", pw_col))
  
  # have a check for wheter stat_name exists in D?
  
  pw_id_map <-
    meta_D$pathways[[pw_col]] %>%
    dplyr::distinct(ID, pathway_name)
  
  
  geneset <-
    rowData(D) %>%
    as.data.frame() %>%
    dplyr::select(COMP_IDstr, !!rlang::sym(pw_col)) %>%
    dplyr::mutate(met_ID=rownames(.)) %>%
    dplyr::filter(!!rlang::sym(pw_col) != "NULL") %>%
    tidyr::unnest_longer(!!rlang::sym(pw_col)) %>%
    dplyr::filter(!!rlang::sym(pw_col) != "NA") %>%
    dplyr::distinct()
  
  
  # Algorithm summary:
  # - calculate metabolite p-values per group
  # - adjust p-values using FDR
  # - assign significance to adjusted p-value at 0.05 level
  # - perform Fisher's exact test
  # - adjust Fisher's exact test p-value
  # - calculate mean fold change based on mean log value of cases and ctrls
  
  enrichment_results <-
    mtm_get_stat_by_name(D, stat_name) %>%
    
    # assign significance
    dplyr::mutate(significant = dplyr::if_else(p.adj < cutoff, TRUE, FALSE),
                  n_total = dplyr::n(),
                  n_total_sig = sum(significant)) %>%
    dplyr::inner_join(geneset, by = c("var" = "met_ID")) %>%
    dplyr::group_by(!!rlang::sym(pw_col)) %>%
    
    # calculate summary numbers for Fisher's test
    dplyr::summarise(n_total = unique(n_total),
                     n_total_sig = unique(n_total_sig),
                     n_pw = dplyr::n(),
                     n_pw_sig = sum(significant),
                     mean_fc = mean(fc)) %>%
    dplyr::filter(n_pw >= 5) %>%
    
    # calculate contingency table entries
    dplyr::mutate(s_p = n_pw_sig,
                  ns_p = n_pw - n_pw_sig,
                  s_np = n_total_sig - n_pw_sig,
                  ns_np = n_total - (s_p + ns_p + s_np)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p_value =
                    stats::fisher.test(matrix(c(s_p, s_np, ns_p, ns_np), nrow = 2)) %>%
                    .$p.value) %>%
    dplyr::ungroup() %>%
    dplyr::rename(ID = !!rlang::sym(pw_col)) %>%
    dplyr::left_join(pw_id_map, by = "ID") %>%
    dplyr::transmute(pathway_name,
                     pathway_ID = ID,
                     p_value,
                     p_value_adjusted = stats::p.adjust(p_value, method = "fdr"),
                     mean_foldchange = mean_fc) %>%
  dplyr::arrange(p_value)
  
  metadata(D)$pathways$enrichment_results <-
    dplyr::as_tibble(enrichment_results)
  
  
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("performed pathway enrichment on %s pathways using Fihser's exact test",
                       nrow(enrichment_results))
    )
  
  D
}



mt_plots_sample_boxplot_new <- function(D,
                                        title = "Sample boxplot",
                                        show_legend = TRUE,
                                        ylabel = "Feature concentrations",
                                        plot_logged = FALSE,
                                        ggadd = NULL,
                                        ...) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("plottitle" %in% dot_args) stop("You used the old MT naming convention plottitle. Should be: title")
  if ("legend" %in% dot_args) stop("You used the old MT naming convention legend. Should be: show_legend")
  if ("manual_ylab" %in% dot_args) stop("You used the old MT naming convention manual_ylab. Should be: ylabel")
  if ("manual_ylabel" %in% dot_args) stop("You used the old MT naming convention manual_ylabel. Should be: ylabel")
  if ("logged" %in% dot_args) stop("You used the old MT naming convention logged. Should be: plot_logged")
  
  # plot_logged?
  Dplot = D
  if (plot_logged) {
    assay(Dplot) <- log2(assay(Dplot))
    ylabel = sprintf("%s [log2]", ylabel)
  }
  
  # merge with sample annotations, only keep the ones that were actually used
  cd <- Dplot %>% colData() %>% as.data.frame() %>% tibble::rownames_to_column("merge.primary")
  keep <- c(mti_extract_variables(quos(...)), "merge.primary")
  cd <- cd[, colnames(cd) %in% keep, drop = FALSE]
  df <- cbind(cd, t(assay(Dplot)))
  
  # generate ggplot
  p <- df %>% tidyr::gather(metab, value, dplyr::one_of(rownames(Dplot))) %>%
    ggplot(aes(x = merge.primary, y = value, ...)) +
    geom_boxplot() +
    ylab(ylabel) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # remove legend?
  if (!show_legend) p = p + theme(legend.position = "none")
  
  # add custom elements?
  if (!is.null(ggadd)) p <- p + ggadd
  
  # Convert ggplot object to plotly object
  interactive_plot <- plotly::ggplotly(p)
  
  # Return the interactive plot
  interactive_plot
}


mt_plots_volcano_new <- function(D,
                                 x = fc,
                                 stat_name,
                                 feat_filter = p.value < 0.05,
                                 xlabel=gsub("~","",as.character(x)),
                                 vline=NA,
                                 hline,
                                 ggadd=NULL,
                                 ...) {
  
  x <- dplyr::enquo(x)
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_name))
    stop("stat_name must be given for volcano plot")
  
  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("statname" %in% dot_args) stop("You used the old MT naming convention statname Should be: stat_name")
  if ("metab_filter" %in% dot_args) stop("You used the old MT naming convention metab_filter Should be: feat_filter")
  
  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))
  # remove rows not needed for plotting
  vars <- c(mti_extract_variables(c(dplyr::enquo(x), dplyr::enquo(feat_filter), quos(...))),"var","name")
  rd <- rd[,colnames(rd) %in% vars,drop=F]
  
  
  ## stat
  data_plot <- mtm_get_stat_by_name(D, stat_name)
  if(quo_name(x) %in% colnames(data_plot)==F) stop(glue::glue("Column {quo_name(x)} not found in stat table."))
  data_plot %<>% dplyr::inner_join(rd, by = "var") %>%
    dplyr::mutate(xxx = !!x)
  
  ## SCALE -log10
  reverselog_trans <- function (base = exp(1)){
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
  }
  
  ## determine if and where to draw hline
  if (!missing(hline)) {
    hliney <- mtm_get_stat_by_name(D, stat_name) %>%
      dplyr::inner_join(rd, by = "var") %>%
      dplyr::mutate(xxx = !!x) %>% dplyr::filter(!!dplyr::enquo(hline)) %>% .$p.value %>% max()
  } else {
    hliney <- NA
  }
  
  ## sanity check that there is something to plot
  if (all(is.na(data_plot$p.value))) stop("All p-values for Volcano plot are NA")
  
  ## CREATE PLOT
  p <- data_plot %>%
    ## do plot
    ggplot(aes(x = xxx, y = p.value)) +
    # vline?
    (if(!is.na(vline)){geom_vline(xintercept = c(-vline, vline), linetype='dashed', color='#F8766D')}else{NULL}) +
    # hline?
    (if(!is.na(hliney)){geom_hline(yintercept = hliney, linetype='dashed', color='#F8766D')}else{NULL}) +
    # points
    geom_point(aes(...)) +
    scale_y_continuous(trans = reverselog_trans(10),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x = xlabel, y = "p-value") +
    ggtitle(stat_name)
  
  ## ADD FEATURE LABELS
  if(!missing(feat_filter)){
    #mti_logstatus("add label")
    feat_filter_q <- dplyr::enquo(feat_filter)
    data_annotate <- data_plot %>%
      dplyr::filter(!!feat_filter_q)
    p <- p + ggrepel::geom_text_repel(data = data_annotate,
                                      aes(label = name), max.overlaps = Inf)
  }
  
  ## ADD AXIS GROUPS
  d <- mtm_get_stat_by_name(D, stat_name, fullstruct=T)
  if ("groups" %in% names(d) && length(d$groups)==2) {
    p <- mti_add_leftright_gg(p, paste0(d$groups[1],' high'), paste0(d$groups[2],' high'))
  }
  
  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd
  
  ## add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("volcano plot, aes: %s", mti_dots_to_str(...)),
      output = p
    )
  
  # Convert ggplot object to plotly object
  interactive_plot <- plotly::ggplotly(p)
  
  # Return the interactive plot
  return(interactive_plot)
  
  
}

mti_get_setting <-  function(D, sname) {
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  # check if settings missing
  if (!("settings" %in% names(metadata(D)))) {
    # no, use default
    # build list
    fulllist <- mti_settings_list()
    lst <- names(fulllist) %>% lapply(function(s){fulllist[[s]]$default})
    names(lst) <- names(fulllist)
  } else {
    # yes, take from metadata
    lst <- metadata(D)$settings
  }
  # error control
  if (!(sname %in% names(lst))){stop(glue::glue("Invalid pipeline setting name: '{sname}'"))}
  # return
  lst[[sname]]
}

mti_add_leftright_gg <- function(gg, left, right) {
  
  # with backward compatibility
  if (utils::compareVersion(as.character(utils::packageVersion("ggplot2")),"3.3.0")>=0) { # at least 3.3.0
    # ggplot2 version >= 3.3.0
    ggbld <- ggplot2::ggplot_build(gg)
    xticks =  ggbld$layout$panel_params[[1]]$x$get_breaks() # needed for >=3.3.0
    xticks.minor = ggbld$layout$panel_params[[1]]$x$get_breaks_minor() # needed for >=3.3.0
    xlims = ggbld$layout$panel_params[[1]]$x.range
    
    # add positions of annotation labels
    xticks = c(xlims[1], xticks, xlims[2])
    # get breaks labels
    xtlabs = ggbld$layout$panel_params[[1]]$x$get_labels() # needed for >=3.3.0
    
  } else {
    # ggplot2 version < 3.3.0
    ggbld <- ggplot2::ggplot_build(gg)
    xticks = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.major_source # needed for <3.3.0
    xticks.minor = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.minor_source # needed for <3.3.0
    xlims = ggbld$layout$panel_params[[1]]$x.range
    
    # add positions of annotation labels
    xticks = c(xlims[1], xticks, xlims[2])
    # get breaks labels
    xtlabs = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.labels # needed for <3.3.0
  }
  
  # align with \n
  txt <- c(left, right)
  txt = paste0("\n", txt)
  xtlabs = paste0(xtlabs, "\n")
  xtlabs = c(txt[1], xtlabs, txt[2])
  
  # return
  gg + ggplot2::scale_x_continuous(breaks = xticks, labels = xtlabs, minor_breaks = xticks.minor)
  
}



mti_extract_variables <- function(lst) {
  # filter down only to the variables needed for plotting
  # need to parse x and ... list
  # browser()
  # q <- quos(...)
  vars = c()
  if (length(q) > 0) {
    vars <- lst %>% unlist() %>% lapply(function(x){x %>% all.vars()}) %>% unlist() %>% as.vector()
  }
  vars <- unique(vars)
  
  # return
  vars
}




mt_plots_net_new <- function(D,
                             stat_name,
                             cor_filter = p.value < 0.05,
                             node_coloring,
                             html_outfile,
                             height = 500) {
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_name))
    stop("stat_name must be given to plot the network")
  
  ## rowData
  rd1 <- subset(rowData(D), select=which(names(rowData(D))=="name")) %>%
    as.data.frame() %>%
    dplyr::mutate(var1 = rownames(D))
  colnames(rd1)[colnames(rd1)=="name"] <- "name1"
  rd1$name1 %<>% make.names()
  rd2 <- subset(rowData(D), select=which(names(rowData(D))=="name")) %>%
    as.data.frame() %>%
    dplyr::mutate(var2 = rownames(D))
  colnames(rd2)[colnames(rd2)=="name"] <- "name2"
  rd2$name2 %<>% make.names()
  
  ## stat
  data_plot <- mtm_get_stat_by_name(D, stat_name) %>%
    dplyr::inner_join(rd1, by = "var1") %>%
    dplyr::inner_join(rd2, by = "var2")
  
  # define node attributes
  nodes <- as.data.frame(unique(rbind(cbind(ids=data_plot$var1,label=data_plot$name1),cbind(ids=data_plot$var2,label=data_plot$name2))))
  
  if(!(missing(node_coloring))) {
    test <- mtm_get_stat_by_name(D, node_coloring);
    test <- test[match(nodes$ids, test$var),];
    nodes$map <- sign(test$statistic)*log10(test$p.value);
    nodes$node_color <- grDevices::colorRampPalette(c('blue', 'white', 'red'))(length(nodes$map))[rank(nodes$map)]
  } else nodes$node_color <- rep("lightblue",times=length(nodes$ids))
  
  ## apply filter on correlations
  if(!missing(cor_filter)){
    mti_logstatus("filter correlations")
    cor_filter_q <- dplyr::enquo(cor_filter)
    data_plot <- data_plot %>%
      dplyr::filter(!!cor_filter_q)
  }
  
  ## define edge attributes
  # rescale correlation to [0 10]
  cor.scaled <- (abs(data_plot[,which(colnames(data_plot)=="statistic")])-min(abs(data_plot[,which(colnames(data_plot)=="statistic")])))/(max(abs(data_plot[,which(colnames(data_plot)=="statistic")]))-min(abs(data_plot[,which(colnames(data_plot)=="statistic")])))*10
  
  edges <- data.frame(from = data_plot[,which(colnames(data_plot)=="name1")], to = data_plot[,which(colnames(data_plot)=="name2")], value = cor.scaled)
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$color <- "black"
  edges$color[edges$value<0] <- "red"
  
  ## plot
  e <- edges
  n <- data.frame(id=nodes$label, label= nodes$label, color=nodes$node_color)
  p_vis <- visNetwork::visNetwork(n,e, height = height, width = "100%") %>%
    # disable physics
    visNetwork::visIgraphLayout(layout = "layout_nicely", physics = F, smooth = F)
  
  # create ggnetwork object
  df <- list()
  df$edges <- data.frame(from = data_plot[,which(colnames(data_plot)=="var1")], to = data_plot[,which(colnames(data_plot)=="var2")], weight = cor.scaled/10)
  df$vertices <- nodes$label
  
  # THIS FUNCTION MUST BE REPLACED - dils is archived
  adj <- as.matrix(dils::AdjacencyFromEdgelist(df$edges, check.full = TRUE))
  colnames(adj[[1]])<- as.character(nodes$label[match(adj[[2]],nodes$ids)])
  rownames(adj[[1]])<- as.character(nodes$label[match(adj[[2]],nodes$ids)])
  mm.net <- network::network(adj[[1]], layout = "kamadakawai", directed = FALSE)
  
  if(!missing(node_coloring)){
    test <- mtm_get_stat_by_name(D, node_coloring)
    test <- test[match(adj[[2]], test$var),]
    map <- sign(test$statistic)*log10(test$p.value)
    mm.net %v% "strength" <- map
  } else {
    mm.net %v% "strength" <- 1
  }
  
  mm.col <- c("positive" = "#000000", "negative" = "#0000FF")
  x <- data_plot$statistic
  x[x>=0] <- "positive"
  x[x<0] <- "negative"
  mm.net %e% "pcor" <- abs(data_plot[,which(colnames(data_plot)=="statistic")])
  mm.net %e% "pos" <- x
  
  # add edge color (positive/negative)
  # add black circle around nodes
  p <- ggplot(mm.net, aes(x, y, xend = xend, yend = yend)) +
    geom_edges(color="black",aes(size=pcor)) +
    geom_nodes(aes(color = strength), size = 7) +
    scale_color_gradient2(low = "#0000FF", mid="white",high = "#FF0000") +
    geom_nodetext(color="grey50",aes(label = vertex.names),
                  size = 3, vjust = -0.6) +
    theme_blank() +
    theme(legend.position = "bottom")
  
  # if html_outfile given, save visnetwork to html
  if (!missing(html_outfile)) {
    # due to odd visSave path handling behavior, we need to export to a tmp file first and the move to final location
    tmpfile = sprintf("tmp_%s.html", uuid::UUIDgenerate())
    visNetwork::visSave(graph = p_vis, file = tmpfile)
    file.rename(tmpfile, html_outfile)
  }
  
  ## add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Correlation Network, aes: %s", stat_name),
      output = list(p),
      output2 = p_vis
    )
  ## return
  D
}





mti_funargs <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))
  
  for(i in dplyr::setdiff(names(formals), names(call)))
    call[i] <- list( formals[[i]] )
  
  # assemble results
  raw <-  as.list(match.call(sys.function(sys.parent(1)), call))
  if(typeof(raw[[1]])=="language") raw[[1]] <- as.character(raw[[1]])[[3]]
  if(FALSE){
    #if(any(sapply(raw[-1],typeof)=="symbol")){
    res <- list(
      fun = strsplit(gsub(".*::", "", as.character(raw[[1]])), '_')[[1]],
      args = sapply(raw[-1], function(x){if(is.symbol(x)&& as.character(x)!="."){eval.parent(x)}else{x}})
    )
  }else{
    res <- list(
      fun = strsplit(gsub(".*::", "", as.character(raw[[1]])), '_')[[1]],
      args = raw[-1]
    )
  }
  # remove "mt"
  res$fun = res$fun[res$fun!="mt"]
  # make sure D does not exist in args
  res$args$D = NULL
  # return
  res
  
}

mt_plots_box_scatter_new <- function(D,
                                     x,
                                     stat_name,
                                     plot_type,
                                     correct_confounder,
                                     feat_filter = p.value < 0.05,
                                     feat_sort = p.value,
                                     annotation = "{sprintf('P-value: %.1e', p.value)}",
                                     full_info = F,
                                     text_size = 3.88,
                                     jitter = "beeswarm",
                                     restrict_to_used_samples = T,
                                     ylabel=NULL,
                                     fit_line = T,
                                     fit_line_se = T,
                                     ggadd = NULL,
                                     pages = F,
                                     ...){
  
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(x)) stop("x must be provided")
  x <- dplyr::enquo(x)
  
  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("metab_filter" %in% dot_args) stop("You used the old MT naming convention metab_filter. Should be: feat_filter")
  if ("metab_sort" %in% dot_args) stop("You used the old MT naming convention metab_sort. Should be: feat_sort")
  if ("rows" %in% dot_args) stop("\"rows\" is no longer an acceptecd argument.")
  if ("cols" %in% dot_args) stop("\"cols\" is no longer an accepted argument.")
  if ("manual_ylab" %in% dot_args) stop("You used the old MT naming convention manual_ylab. Should be: ylabel")
  if ("manual_ylabel" %in% dot_args) stop("You used the old MT naming convention manual_ylabel. Should be: ylabel")
  if ("fitline" %in% dot_args) stop("You used the old MT naming convention fitline. Should be: fit_line")
  if ("fitline_se" %in% dot_args) stop("You used the old MT naming convention fitline_se. Should be: fit_line_se")
  
  # create dummy SE so original not changed
  Ds <- D
  
  ## CONFOUNDER
  if(!missing(correct_confounder)){
    #mti_logstatus(glue::glue("correcting for {correct_confounder}"))
    Ds <- mti_correctConfounder(Ds, correct_confounder)
  }
  
  ## rowData
  rd <- rowData(Ds) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(Ds))
  
  ## stat
  if(!missing(stat_name)){
    stat <- maplet::mtm_get_stat_by_name(Ds, stat_name) %>%
      dplyr::inner_join(rd, by = "var")
  }else{
    stat <- rd
    ### KC: ONLY IN BOXPLOT (Why?)
    restrict_to_used_samples <- F # not dependend on a stat
  }
  
  ## FILTER FEATURES
  ### KC: feat_filter is never not missing, should there be no default?
  if(!missing(feat_filter)){
    feat_filter_q <- dplyr::enquo(feat_filter)
    stat <- stat %>%
      dplyr::filter(!!feat_filter_q)
    #mti_logstatus(glue::glue("filter features: {feat_filter_q} [{nrow(stat)} remaining]"))
  }
  
  ## SORT FEATURES
  ### KC: feat_sort is never not missing, should there be no default?
  if(!missing(feat_sort)){
    feat_sort_q <- dplyr::enquo(feat_sort)
    stat <- stat %>%
      dplyr::arrange(!!feat_sort_q) %>%
      ## sort according to stat
      dplyr::mutate(name = factor(name, levels = unique(name)))
    #mti_logstatus(glue::glue("sorted features: {feat_sort_q}"))
  }
  
  ## CREATE PLOT
  dummy <- Ds %>%
    mti_format_se_samplewise() %>% # NOTE: No explosion of dataset size due to active restriction - 6/2/20, JK
    tidyr::gather(var, value, dplyr::one_of(rownames(Ds)))
  ## filter to groups?
  if(plot_type=="box"){
    if (restrict_to_used_samples) {
      filterto <- maplet::mtm_get_stat_by_name(Ds, stat_name, fullstruct=T)$samples.used
      dummy <- dummy[filterto,]
    }
  }
  
  # check x is a column in dataset
  mainvar <- x %>% dplyr::quo_name()
  if(mainvar %in% colnames(dummy) == F) stop(glue::glue("No column in plot data frame with name \"{mainvar}\"."))
  
  if(!full_info){
    # filter down only to the variables needed for plotting
    # need to parse x and ... list
    vars <- x %>% dplyr::quo_name()
    q <- dplyr::quos(...)
    if (length(q) > 0) {
      vars <- c(vars, q %>% lapply(function(x){x %>% as.character() %>% gsub("~","",.)}) %>% unlist() %>% as.vector())
    }
    vars <- unique(vars)
    
    plottitle <- ifelse(missing(stat_name),"",stat_name)
    if(plot_type=="box"){
      # make sure the main outcome variable x is a factor
      mainvar <-x %>% dplyr::quo_name()
      dummy[[mainvar]] <- as.factor(dummy[[mainvar]])
      
      p <- dummy %>%
        dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
        ## add feature names, but only restricted subset from statistics table
        dplyr::inner_join(stat[,dplyr::intersect(colnames(stat),c('var','statistic','p.value','p.adj','name'))], by = "var") %>%
        dplyr::select(-var) %>%
        ## do plot
        ggplot() +
        geom_boxplot(aes(x = as.factor(!!x), y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
        labs(x = NULL, y = NULL) +
        ggtitle(plottitle)
    }else{
      
      p <- dummy %>%
        dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
        ## add feature names, but only restricted subset from statistics table
        dplyr::inner_join(stat[,c('var','statistic','p.value','p.adj','name')], by = "var") %>%
        dplyr::select(-var) %>%
        ## do plot
        ggplot() +
        ## add fit line?
        {if (fit_line) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fit_line_se, color = "black") else NULL} +
        geom_point(aes(x = !!x, y = value, ...)) +
        labs(x = dplyr::quo_name(x), y="feature") +
        ggtitle(plottitle)
      
    }
    
  }else{
    
    if(plot_type=="box"){
      # leave full info in
      # can create huge data.frames
      
      plottitle <- ifelse(missing(stat_name),"",stat_name)
      p <- dummy %>%
        ## add feature names
        dplyr::inner_join(stat, by = "var") %>%
        ## do plot
        ggplot() +
        geom_boxplot(aes(x = !!x, y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
        labs(x = NULL, y = NULL) +
        ggtitle(plottitle)
      
    }else{
      # leave full info in
      # can create huge data.frames
      
      #
      plottitle <- ifelse(missing(stat_name),"",stat_name)
      p <- dummy %>%
        ## add feature names
        dplyr::inner_join(stat, by = "var") %>%
        ## do plot
        ggplot() +
        ## add fit line?
        {if (fit_line) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fit_line_se, color = "black") else NULL} +
        geom_point(aes(x = !!x, y = value, ...)) +
        labs(x = dplyr::quo_name(x), y="feature") +
        ggtitle(plottitle)
      
    }
    
  }
  
  ### BOX PLOT SPECIFIC
  if(plot_type=="box"){
    ## add ylabel
    if (!is.null(ylabel)) {
      p <- p + ylab(ylabel)
    } else {
      # add label if this is logged data
      r <- Ds %>% maplet::mtm_res_get_entries(c("pre","trans","log"))
      if (length(r)>0) {
        p <- p + ylab(r[[1]]$logtxt) # log text contains e.g. "log2"
      }
    }
    
    ## ADD JITTER
    if(jitter=="beeswarm"){
      p <- p +
        ggbeeswarm::geom_beeswarm(aes(x = !!x, y = value, ...))
    }else if(jitter=="jitter"){
      p <- p +
        geom_jitter(aes(x = !!x, y = value, ...))
    }
  }
  
  
  ### COMMON TO BOTH PLOTS
  ## ADD ANNOTATION
  if(!missing(annotation)){
    data_annotate <- stat %>%
      dplyr::mutate(annotate = glue::glue(annotation)) %>%
      dplyr::distinct(name, annotate)
    p <- p + geom_text(data = data_annotate,
                       aes(label = annotate),
                       x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size=text_size )
  }
  
  if (!is.null(ggadd)) p <- p+ggadd
  
  funargs <- mti_funargs()
  
  
  
  ## SPLIT TO MULTIPLE PAGES
  # if there is no plot, create a single empty page
  if (length(unique(stat$name))==0) {
    print("this 1")
    
    p <- ggplot() + geom_text(aes(x = 0, y = 0, label = 'no plots'), size = 10)
    output2 <- NULL
    
    D %<>%
      mti_generate_result(
        funargs = funargs,
        logtxt = sprintf("Feature ",ifelse(plot_type=="box", "boxplots", "scatter plots"),", aes: %s", mti_dots_to_str(...)),
        output = p,
        output2 = output2
      )
    
    interactive_plot <- plotly::ggplotly(p)
    
    # Return the interactive plot
    interactive_plot
    
  } else if(pages) {
    
    print("this 2")
    npages <- ceiling(length(unique(stat$name))/8)
    tmp <- lapply(1:npages, function(x){
      p + ggforce::facet_wrap_paginate(.~name, scales = "free", ncol = 2, nrow = 4, page=x)
    })
    p <- ggpubr::ggarrange(plotlist = tmp, ncol=1)
    p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(plottitle))
    output2 <- ceiling(length(unique(stat$name))/2)
    
    D %<>%
      mti_generate_result(
        funargs = funargs,
        logtxt = sprintf("Feature ",ifelse(plot_type=="box", "boxplots", "scatter plots"),", aes: %s", mti_dots_to_str(...)),
        output = p,
        output2 = output2
      )
    
    # Convert ggplot object to plotly object
    interactive_plot <- plotly::ggplotly(p)
    
    # Return the interactive plot
    interactive_plot
    
  }else {
    p <- p + facet_wrap(.~name, scales = "free", ncol=2)
    #p <- list(p)
    output2 <- ceiling(length(unique(stat$name))/2)
    
    D %<>%
      mti_generate_result(
        funargs = funargs,
        logtxt = sprintf("Feature ",ifelse(plot_type=="box", "boxplots", "scatter plots"),", aes: %s", mti_dots_to_str(...)),
        output = list(p),
        output2 = output2
      )
    
    # Convert ggplot object to plotly object
    interactive_plot <- plotly::ggplotly(p)
    
    # Return the interactive plot
    interactive_plot
  
  }
  
  
}

mti_correctConfounder <- function(D, formula){
  d <- D %>% mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
  d_cor <- rownames(D) %>%
    purrr::map_dfc(function(m){
      f   <- stats::update.formula(formula, stringr::str_c(m, "~."))
      mod <- stats::lm(f, data = d, na.action = na.exclude)
      res <- stats::resid(mod)
      res
    }) %>%
    stats::setNames(rownames(D)) %>%
    as.matrix() %>% t()
  colnames(d_cor) <- colnames(D)
  assay(D)        <- d_cor
  D
}

mti_format_se_samplewise <- function(D){
  # coldata and assay cannot have overlapping names
  # (this should be caugh earlier in the pipeline, but here is where it causes trouble)
  inters <- intersect(colnames(colData(D)), rownames(D))
  if (length(inters)>0) {
    stop(sprintf("There are metabolites and colData variables with the same name: %s", paste0(inters, collapse = ", ")))
  }
  # cbind
  cbind(colData(D),
        t(assay(D))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("merge.primary")
}


mt_plots_stats_pathway_bar_new <- function(D,
                                       stat_list,
                                       feat_filter = p.value < 1,
                                       group_col = "SUB_PATHWAY",
                                       color_col = NULL,
                                       y_scale = "fraction",
                                       sort_by_y = FALSE,
                                       assoc_sign_col,
                                       add_empty = FALSE,
                                       keep_unmapped = FALSE,
                                       outfile = NULL,
                                       ggadd = NULL,
                                       ...){
  
  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("stat_name" %in% dot_args) stop("You used the old MT naming convention stat_name. Should be: stat_list.")
  if ("metab_filter" %in% dot_args) stop("You used the old MT naming convention metab_filter. Should be: feat_filter.")
  if ("aggregate" %in% dot_args) stop("You used the old MT naming convention aggregate Should be: group_col.")
  if ("colorby" %in% dot_args) stop("You used the old MT naming convention colorby Should be: color_col.")
  if ("sort" %in% dot_args) stop("You used the old MT naming convention sort. Should be: sort_by_y.")
  if ("yscale" %in% dot_args) stop("You used the old MT naming convention yscale. Should be: y_scale.")
  if ("assoc_sign" %in% dot_args) stop("You used the old MT naming convention assoc_sign. Should be: assoc_sign_col.")
  if ("keep.unmapped" %in% dot_args) stop("You used the old MT naming convention keep.unmapped. Should be: keep_unmapped.")
  if ("output.file" %in% dot_args) stop("You used the old MT naming convention output.file. Should be: outfile.")
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_list) & !missing(feat_filter))
    stop("stat_list must be given for feat_filter to work.")
  if(missing(stat_list) & !missing(assoc_sign_col))
    stop("stat_list must be given for assoc_sign_col to work.")
  if(!(group_col %in% colnames(rowData(D))))
    stop(sprintf("group_col column '%s' not found in rowData", group_col))
  if(!is.null(color_col))
    if(!(color_col %in% colnames(rowData(D))))
      stop(sprintf("color_col column '%s' not found in rowData", color_col))
  
  
  print("WOFFFFFFFFFFFFFF HRERERER")
  
  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))
  # set the nulls to unknown
  if(keep_unmapped){
    rd[[group_col]][which(rd[[group_col]]=="NULL")] <- "Unmapped"
  } else{
    rd <- rd[which(rd[[group_col]]!="NULL"), ]
  }
  perc <- rd[[group_col]] %>%
    unlist %>% table(exclude = NULL) %>% as.data.frame()
  colnames(perc) <- c("name","count")
  
  flag_filter <- ifelse((!missing(feat_filter)), T,F)
  flag_sign <- ifelse((!missing(assoc_sign_col)), T,F)
  
  data <- lapply(stat_list %>% {names(.)=.;.}, function(ss){
    ## subselect variables
    # if(flag_filter) {
    #   feat_filter_q <- dplyr::enquo(feat_filter)
    #   sel <- maplet::mtm_get_stat_by_name(D=D,name=ss) %>%
    #     dplyr::filter(!!feat_filter_q) %>%
    #     dplyr::filter(var %in% rd$var)
    #   rd <- rd %>%
    #     dplyr::filter(var %in% sel$var)
    #   
    # }
    
    #my changes
    if (flag_filter) {
      feat_filter_q <- enquo(feat_filter)
      sel <- mtm_get_stat_by_name(D = D, name = ss) %>%
        filter(!!feat_filter_q) %>%
        filter(var %in% rd$var)
      
      print(sel)  # Print the sel object to check its content
      
      if (nrow(sel) > 0) {
        rd_filtered <- rd %>%
          filter(var %in% sel$var)
        
        if (nrow(rd_filtered) > 0) {
          rd <- rd_filtered
        }
      }
    }
    
    
    
    
    print("GGGYGYGYGYGYGYG HRERERER")
    
    # if filtering gives an empty matrix, produce an empty df
    if(nrow(rd)==0) {
      data_plot <- data.frame(name=as.character(),
                              count = as.numeric())
      anno <- data.frame()
      
    } else {
      # if assoc_sign_col given, include in data
      if(flag_sign){
        if(!(assoc_sign_col %in% colnames(sel))) {
          stop(sprintf("Could not find column called %s in the statistical results called %s", assoc_sign_col, stat_list))
        } else {
          # reorder sel according to rd
          sel <- sel[match(sel$var,rd$var),] %>%
            dplyr::mutate(association=ifelse(sign(!!sym(assoc_sign_col))>0, "positive", "negative")) %>%
            dplyr::select(var,association)
          
          data_plot <- data.frame(name=rd[[group_col]] %>% unlist %>% as.vector,
                                  association=rep(sel$association, times= (rd[[group_col]] %>% sapply(length)))) %>%
            table(exclude = NULL) %>% as.data.frame()
          colnames(data_plot) <- c("name","association","count")
          
        }
      } else {
        # reorder sel according to rd
        sel <- sel[match(sel$var,rd$var),] %>%
          dplyr::select(var)
        
        data_plot <- rd[[group_col]] %>%
          unlist %>% table(exclude = NULL) %>% as.data.frame()
        colnames(data_plot) <- c("name","count")
        
      }
      
      if(add_empty){
        # check which group_col entries are not included
        agg <- rowData(D) %>% as.data.frame() %>% .[[group_col]] %>% unlist %>% unique
        agg_empty <- agg[which(!(agg %in% unique(as.character(data_plot$name))))]
        
        # create data frame
        empty <- data.frame(name = agg_empty, count = rep(0, times=length(agg_empty)))
        
        if("association" %in% colnames(data_plot)){
          empty$association <- "positive"
        }
        
        # add to data
        data_plot %<>% dplyr::full_join(empty, by=colnames(data_plot))
      }
      
      # add number of features in each pathway
      perc <- data_plot %>% dplyr::select(name) %>%
        dplyr::left_join(perc, by="name")
      data_plot <- data_plot %>%
        # add fraction variable
        dplyr::mutate(fraction= count/perc$count)
      
      # add color column if not given
      if(is.null(color_col)) {
        color_col <- paste(group_col,"color", collapse = "_")
        rd[[color_col]] <- "pathway"
      }
      
      # create dictionary between group_col and color_col variables
      dict <- rd %>% dplyr::select(!!sym(group_col),!!sym(color_col)) %>% tidyr::unnest_longer(col=group_col) %>% as.data.frame()
      dict <- dict[!duplicated(dict[[group_col]]),]
      
      # add color to data_plot
      data_plot <- data_plot %>%
        dplyr::left_join(dict, by=c("name"=group_col)) %>%
        dplyr::rename(color=sym(color_col))
      
      # create annotation data
      anno <- data.frame(name = rep(rd$name, times=sapply(rd[[group_col]], length) %>% as.vector()),
                         var = rep(rd$var, times=sapply(rd[[group_col]], length) %>% as.vector()),
                         pathway = unlist(rd[[group_col]]),
                         color = if(!is.null(color_col)){rep(rd[[color_col]], times=sapply(rd[[group_col]], length) %>% as.vector())}else{"pathway"}) %>%
        dplyr::left_join(maplet::mtm_get_stat_by_name(D=D,name=ss) , by="var") %>%
        dplyr::select(-var)
      
      print("WVUVUVUVUVOOOO HRERERER")
      
      # if pathway mapping exists in the metadata, use the names provided there
      x <- D %>% metadata
      if ("pathways" %in% names(x)){
        if (group_col %in% names(x$pathways)) {
          # add pathway names to dataframe
          data_plot %<>%
            dplyr::left_join(x$pathways[[group_col]][,c("ID","pathway_name")], by=c("name"="ID"))
          anno %<>%
            dplyr::left_join(x$pathways[[group_col]][,c("ID","pathway_name")], by=c("pathway"="ID")) %>%
            dplyr::select(name,pathway_name,pathway,color,everything()) %>%
            dplyr::rename(pathway_id=pathway, pathway=pathway_name)
          
          # set Unknown pathway names to Unknown
          if(length(which(is.na(data_plot$pathway_name)))>0){
            data_plot$pathway_name[which(is.na(data_plot$pathway_name))] <- "Unknown"
          }
          # substitute codes for names and remove extra column
          data_plot %<>%
            dplyr::mutate(name=pathway_name) %>%
            dplyr::select(-pathway_name)
        } else{
          mti_logwarning(sprintf("%s field not found in the metadata",group_col))
        }
      }
      # create labels for plotting
      data_plot %<>% dplyr::mutate(label=sprintf("%s [%d]", name, perc$count))
      
      # convert labels to factor to sort alphabetically
      data_plot$label <- as.factor(data_plot$label)
      
      # add comparison name to df
      data_plot$comp <- ss
      anno$comp <- ss
      
    }
    list(dt = data_plot, anno = anno)
  })
  
  # function to revert string structure
  revert_list_str_4 <- function(ls) {
    # get sub-elements in same order
    x <- lapply(ls, `[`, names(ls[[1]]))
    # stack and reslice
    apply(do.call(rbind, x), 2, as.list)
  }
  
  data <- revert_list_str_4(data)
  
  print("CYCYCYCY HRERERER")
  
  # if there is at least one result, produce plot, otherwise output empty plot
  if((sapply(data$dt, function(ss){dim(ss)[1]}) %>% sum()) >0) {
    # merge list into a single dataframe
    data_plot <- do.call(rbind, data$dt) %>% as.data.frame()
    # get common colnames in stat table
    tt <- sapply(data$anno, colnames) %>% unlist %>% table %>% as.data.frame
    colnames(tt)[1] <- "var"
    tt %<>%
      dplyr::filter(Freq == max(Freq)) %>%
      dplyr::pull(var) %>%
      as.character
    anno <- lapply(data$anno, function(x){
      x %>%
        dplyr::select(tidyselect::any_of(tt))
    }) %>% {do.call(rbind, .)} %>% as.data.frame()
    
    # optional sorting (only for single statistical results)
    if (sort_by_y){
      data_plot$label <- reorder(data_plot$label, -data_plot[[y_scale]])
    }
    # sort comp so that facets appear in the same order given by the user
    data_plot$comp <- factor(data_plot$comp,levels=stat_list)
    
    # convert count to numeric
    data_plot$count %<>% as.numeric
    
    ## CREATE PLOT
    p <- ggplot(data_plot, aes(label)) +
      (if("association" %in% colnames(data_plot)) {geom_bar(data = subset(data_plot, association == "positive"), aes(y = !!sym(y_scale), fill = color), stat = "identity", position = "dodge", color="black", size=0.4)}) +
      (if("association" %in% colnames(data_plot)) {geom_bar(data = subset(data_plot, association == "negative"), aes(y = -!!sym(y_scale), fill = color), stat = "identity", position = "dodge", color="black", size=0.4)} else{geom_bar(aes(x=label, y=!!sym(y_scale), fill=color), stat = "identity", color="black", size=0.4)}) +
      (if(y_scale=="fraction") {ggtitle(sprintf("Fraction of pathway affected, %s", gsub("~", "", rlang::expr_text(dplyr::enquo(feat_filter)))))}else{ggtitle(sprintf("Number of hits per pathway, %s", gsub("~", "", rlang::expr_text(dplyr::enquo(feat_filter)))))}) +
      (if(y_scale=="count" & "association" %in% colnames(data_plot)) {expand_limits(y=c(-max(data_plot$count, na.rm = T)*1.7, max(data_plot$count, na.rm = T)*1.7))}) +
      (if(y_scale=="count" & !("association" %in% colnames(data_plot))) {expand_limits(y=c(0, max(data_plot$count, na.rm = T)*1.7))}) +
      (if(y_scale=="fraction" & "association" %in% colnames(data_plot)) {expand_limits(y=c(-1, 1))}) +
      geom_hline(yintercept = 0,colour = "black", size=0.4) +
      labs(x="",fill = color_col) +
      theme(plot.title = element_text(hjust = 0.4)) +
      scale_x_discrete(limits = rev(levels(data_plot$label)))
    
    # add phenotype labels to x axis
    if("association" %in% colnames(data_plot) & length(stat_list)==1){
      d <- maplet::mtm_get_stat_by_name(D, stat_list, fullstruct=T)
      if ("groups" %in% names(d) && length(d$groups)==2) {
        # get breaks
        ggbld <- ggplot2::ggplot_build(p)
        yticks = ggbld$layout$panel_params[[1]]$y$minor_breaks # using minor_breaks because sometimes breaks would not work
        # edit labels to include groups
        ytlabs = yticks
        ytlabs[1] <- sprintf("%s\n%s", yticks[1], sprintf("high in %s", d$groups[1]))
        ytlabs[length(ytlabs)] <- sprintf("%s\n%s", yticks[length(yticks)], sprintf("high in %s", d$groups[2]))
        # apply new labels
        p <- p +
          scale_y_continuous(breaks = yticks, labels = ytlabs)
      }
    }
    
    print("WOOOO HRERERER")
    
    # flip axes and add annotations on bars
    p <- p +
      coord_flip() +
      (if(y_scale=="count" & !("association" %in% colnames(data_plot))) {geom_text(data=data_plot, aes(label, !!sym(y_scale), label= sprintf("%.2f%%", fraction*100)),
                                                                                   position = position_dodge(width=0.9), hjust = -0.1, size=2.5)}) +
      (if(y_scale=="count" & "association" %in% colnames(data_plot)) {geom_text(data=dplyr::filter(data_plot, association=="positive"), aes(label, !!sym(y_scale), group= association,label= sprintf("%.2f%%", fraction*100)),
                                                                                position = position_dodge(width=0.9), hjust = -0.1, size=2.5)}) +
      (if(y_scale=="count" & "association" %in% colnames(data_plot)) {geom_text(data=data_plot %>% dplyr::filter(association=="negative"), aes(label, -!!sym(y_scale), group= association,label= sprintf("%.2f%%", fraction*100)),
                                                                                position = position_dodge(width=0.9), hjust = 1.1, size=2.5)}) +
      facet_wrap(~comp)
    
    # add custom elements?
    if (!is.null(ggadd)) p <- p + ggadd
    
    # save plot parameters to be passed to the html generator for dynamical plot height
    re <- p %>%
      ggplot2::ggplot_build() %>%
      magrittr::extract2('layout') %>%
      magrittr::extract2('layout')
    
    nr <- data_plot$name %>% unique %>% length # number of pathways
    ncol <- re$COL %>% max() # number of panel columns
    nrow <- re$ROW %>% max() # number of panel rows
    
  } else {
    p <- ggplot() +
      geom_text(aes(x=0,y=0, label="No significant results"), size=10)
    
    # save plot parameters to be passed to the html generator for dynamical plot height
    nr <- 0 # number of pathways
    ncol <- NULL
    nrow <- NULL
    
  }
  
  print("REACHED HRERERER")
  
  if(!is.null(outfile)){
    if(exists("data_plot")){
      wb = openxlsx::createWorkbook()
      sheet = openxlsx::addWorksheet(wb, "Parameters")
      if(is.null(color_col)){color_col <- 'none'}
      openxlsx::writeData(wb, sheet=sheet, list(comparisons = stat_list, feat_filter = gsub("~", "", rlang::expr_text(dplyr::enquo(feat_filter))), group_col = group_col, coloredby = color_col))
      sheet = openxlsx::addWorksheet(wb, "AggregatedPathways")
      openxlsx::writeData(wb, sheet=sheet, data_plot, rowNames = F, colNames = T)
      sheet = openxlsx::addWorksheet(wb, "IndividualResults")
      openxlsx::writeData(wb, sheet=sheet, anno, rowNames = F, colNames = T)
      openxlsx::saveWorkbook(wb, outfile, overwrite = T)
    } else {
      warning("mt_plots_statsbarplot: No significant results. outfile ignored.")
    }
  }
  
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = ifelse(exists("stat_list"), sprintf("bar plot for comparison %s, by %s, filtered for %s, using %s", paste(stat_list,collapse = ", "), group_col, gsub("~", "", rlang::expr_text(dplyr::enquo(feat_filter))), y_scale),
                      sprintf("bar plot by %s using %s", group_col, y_scale)),
      output = list(p),
      output2 = list(nr = nr, npancol = ncol, npanrow = nrow)
    )
  
  # Convert ggplot object to plotly object
  interactive_plot <- plotly::ggplotly(p)
  
  # Return the interactive plot
  interactive_plot
  
}

mti_settings_list <- function() {
  list(
    dummy = list(class="numeric", default=5),
    save_all_assays = list(class="logical", default=F)
  )
}

mti_add_to_list <- function(lst, element, oname) {
  # init if needed
  if (is.null(lst) || length(lst)==0) lst=c()
  # add and return
  lst[[length(lst)+1]] = element
  names(lst)[length(lst)] = oname
  
  lst
}

mti_logmsg <- function(msg) { logging::loginfo(mti_escape_percent(msg), logger="mt"); msg }
mti_logstatus <- function(msg) { logging::loginfo(mti_escape_percent(msg), logger="mts"); msg }
mti_logwarning <- function(msg) { logging::loginfo(mti_escape_percent(sprintf("WARNING: %s",msg)), logger="mtw"); msg }

# produce sprintf-safe copy of string (create for new loginfo() behaviour, that interprets strings as formats)
mti_escape_percent <- function(txt) gsub('%','%%',txt)
mti_dots_to_str <- function(...) {
  l = eval(substitute(alist(...)))
  paste(sapply(names(l), function(k){sprintf('%s=%s',k,as.character(l[[k]]))}), collapse = ', ')
}


mti_generate_result <- function(
    D,
    funargs,
    logtxt="",
    output=NULL,
    output2=NULL
) {
  
  # ensure structure of funargs
  stopifnot("fun" %in% names(funargs))
  stopifnot("args" %in% names(funargs))
  
  this.uuid = uuid::UUIDgenerate()
  # check if assay updated
  save_assays = mti_get_setting(D, "save_all_assays")
  if(save_assays){
    if(assays(D) %>% length() != 0){
      metadata(D)$assays %<>% mti_assay_ptr(res_assay=assay(D))
      assay_ptr <- metadata(D)$assays$head_ptr
    }else{
      assay_ptr <- NULL
    }
  }else{
    assay_ptr <- NULL
  }
  
  # assemble list
  metadata(D)$results %<>% mti_add_to_list(
    list(
      fun=funargs$fun,
      args=funargs$args,
      logtxt=mti_logmsg(logtxt),
      uuid=this.uuid,
      output=output,
      output2=output2,
      assay_ptr=assay_ptr
    ),
    oname = paste(paste(funargs$fun,collapse = "_"), this.uuid, sep = ".")
  )
  
  D
}





mt_stats_cormat_genenet_new = function(D, stat_name, samp_filter) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))
  
  ## stat_name
  if(missing(stat_name))
    stop("stat_name must be given")
  ## check for NA and throw an error if yes
  if(any(is.na(X)))
    stop("the data matrix contains NAs")
  
  ## FILTER SAMPLES
  if(!missing(samp_filter)) {
    
    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,]
    
  } else {
    samples.used = rep(T, ncol(D))
  }
  
  # filter
  X <- X[samples.used,]
  
  # compute partial correlation using GeneNet
  pcor_GeneNet <- GeneNet::ggm.estimate.pcor(as.matrix(X), method = "dynamic", verbose=FALSE)
  pval_GeneNet <- GeneNet::network.test.edges(pcor_GeneNet, plot=FALSE)
  
  # create result variables
  node1 <- colnames(pcor_GeneNet)[pval_GeneNet$node1]
  node2 <- colnames(pcor_GeneNet)[pval_GeneNet$node2]
  var <- paste0(node1,"_",node2, sep="")
  
  # create result table
  tab <- data.frame("var"=var, "statistic"=pval_GeneNet$pcor, "p.value"=pval_GeneNet$pval, "var1"=node1, "var2"=node2)
  
  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = 'GeneNet partial correlation',
      output = list(
        table = tab,
        name = stat_name,
        lstobj = NULL,
        samples.used = samples.used,
        outcome = NA
      )
    )
  
  # return
  tab
  
}



#' Generates markdown-based HTML output from \code{SummarizedExperiment}
#'
#' @description
#' Generates a fully automated report version of an entire (linear) pipeline.
#'
#' @description
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output HTML filename.
#' @param title Title of RMD document. Default: 'RMD output'.
#' @param output_calls Output detailed info on function calls? Default: F.
#' @param number_sections Number sections and sub-sections? Default: F.
#' @param start_after Tag name or UUID of pipeline step AFTER which to start. Default: NA (i.e. output entire pipeline).
#' @param use_plotly EXPERIMENTAL. Output interactive plotly plots? WARNING: Setting this argument to TRUE can significantly
#'    increase the runtime of the function. If plotly plots are not needed for a report, we recommend keeping this argument set to
#'    FALSE. Default: F.
#' @param keep_tmp Keep temporary files? Can be used to manually edit RMD afterwards. Default: F.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @author JK
#'
#' @export
mt_reporting_html_new <- function(D,
                              file,
                              title = 'RMD output',
                              output_calls=F,
                              number_sections=F,
                              start_after=NA,
                              use_plotly=F,
                              keep_tmp=F,
                              enrichment_results,
                              plotly_file_path,
                              mod3_file_path
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # unique string
  ustr <- uuid::UUIDgenerate()
  
  # define file names
  rmdfile <- sprintf("tmp_%s.RMD", ustr)
  rdsfile <-  sprintf("tmp_%s.rds", ustr) # only used of keep_tmp==T
  
  # generate RMD
  D %>% reporting_generateMD(
    file = rmdfile,
    readfrom = rdsfile,
    title = title,
    output.calls = output_calls,
    number.sections = number_sections,
    start.after=start_after,
    use.plotly = use_plotly,
    keep_tmp = keep_tmp,
    outfile = file,
    enrichment_results = enrichment_results,
    plotly_file_path = plotly_file_path,
    mod3_file_path = mod3_file_path)
  
  # save temp file that will be input for the RMD?
  if (keep_tmp){
    save(D, file=rdsfile)
    # keep temp files, move them to their final destination name
    file.rename(rmdfile, paste0(tools::file_path_sans_ext(file),'.rmd'))
    file.rename(rdsfile, paste0(tools::file_path_sans_ext(file),'.rds'))
    rmdfile <- paste0(tools::file_path_sans_ext(file),'.rmd')
  }
  
  # knit
  rmarkdown::render(rmdfile, params=list(D=D))
  
  # clean up
  if (!keep_tmp) {
    # rename to correct name
    file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), file)
    # no temp files left behind
    file.remove(rmdfile)
    # .rds does not need to be deleted because it was never generated
  }
  
  # return document, in case pipeline is supposed to keep running
  D
}


#' Markdown-based report generator
#'
#' Generates a fully automated report version of an entire (linear) pipeline.
#'
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' @param D \code{SummarizedExperiment} input
#' @param file File to be generated
#' @param readfrom Name of R data file data will be loaded and is supposed to contain SummarizedExperiment "D". Will not actually be loaded in this function, but while knitting the RMD later.
#' @param title Title of RMD document
#' @param firstheading Name of first heading
#' @param use.plotly Output interactive plotly plots? (experimental)
#' @param output.calls Output full information about all parameters of each function call into RMD?
#' @param number.sections Number sections and sub-sections? (default: F)
#' @param start.after UUID of pipeline step AFTER which to start (default: none, i.e. output entire pipeline)
#' @param keep_tmp Keep the intermediate rmarkdown file? Default: F.
#' @param outfile Name of html file to be created. Required if keep_tmp==T.
#'
#' @returns nothing
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @author JK, KC
#'
#' @noRd
reporting_generateMD <- function(
    D,
    file = 'MT.RMD',
    readfrom = 'mt.rds',
    title = 'RMD output',
    firstheading='Output',
    use.plotly=F,
    output.calls=F,
    number.sections=F,
    start.after=NA,
    keep_tmp=F,
    outfile,
    enrichment_results,
    plotly_file_path,
    mod3_file_path
) {
  
  
  #### helper functions
  out <- function(str){writeLines(str, h)}
  
  writechunk <- function(code, params='') {
    if (nchar(params)>0) params=paste0(' ', params)
    out(sprintf("```{r%s}", params))
    out(code)
    out("```\n")
  }
  
  
  #### initialize output file
  h <- file(file, open='wt')
  
  #### markdown header and first heading
  out(glue::glue('
---
title: {title}
output:
  html_document:
    toc: true
    toc_float: TRUE
    {if(number.sections){"number_sections: true"}else{""}}
params:
  D: NA
  enrichment_results: NA
---

    '))
  
  
  
  #### determine where to start from (first pipeline step, or after one)
  start.from = 1 # by default
  if (!is.na(start.after)) {
    # extract all UUIDs
    allids <- D %>% metadata() %>% .$results %>% purrr::map("uuid") %>% unlist()
    
    # if tag name provided, get uuid
    tag_name_list <- maplet::mtm_res_get_entries(D, c("reporting", "tag")) %>% purrr::map("output") %>% unlist()
    tag_name_idx <- match(start.after, tag_name_list)
    if(!is.na(tag_name_idx)){
      start.after <- tag_name_list[tag_name_idx] %>% names() %>% gsub("reporting_tag.", "", .)
    }
    
    # find the one to start after
    start.from <- which(allids==start.after)
    # error-check
    if (length(start.from)==0)
      stop(sprintf("Could not find pipeline step to start after: '%s'", start.after))
    if (D %>% metadata() %>% .$results %>% length() == start.from)
      stop(sprintf("Cannot start after pipeline step '%s', because it's the last entry of the pipeline", start.after))
    start.from <- start.from + 1
  }
  
  
  #### global chunk options
  writechunk("# default chunk options\nknitr::opts_chunk$set(warning=F,echo=F,results='hide',message=F)", params = "echo=F")
  
  #### chunk that loads libraries
  if (!use.plotly) {
    writechunk('# load libraries\nlibrary(maplet)\n')
  } else  {
    writechunk('# load libraries\nlibrary(maplet)\nlibrary("plotly")')
  }
  
  #### chunk that assigns list r from data
  if(keep_tmp){
    rds_file <- paste0(tools::file_path_sans_ext(outfile),".rds")
    writechunk(glue::glue('load("{rds_file}")\nr <- metadata(D)$results'))
  }else{
    writechunk(glue::glue('r <- metadata(D)$results'))
  }
  
  #### start of output
  out(glue::glue('# {firstheading}'))
  
  #### loop over results
  lvl=2 # current level of results = 2
  
  # loop over results
  r <- metadata(D)$results
  
  for (i in start.from:length(r)) {
    # to ignore?
    if (length(r[[i]]$fun)<2 || r[[i]]$fun[2]!="void") { # ignore void
      
      # reporting step?
      if (r[[i]]$fun[1]!="reporting") {
        # not reporting, actual pipeline step
        
        ## header
        out(glue::glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))
        
        ## detailed arguments?
        if (output.calls) {
          L <- r[[i]]$args
          out("*Function arguments:*<br/>")
          out(names(L) %>% lapply(function(x){sprintf("%s=%s",x, toString(L[[x]]))}) %>% paste0(collapse = "<br/>"))
          out("")
        }
        
        ## log text
        out(glue::glue('*Log text:*<br/>{r[[i]]$logtxt}\n\n'))
        
        ## plot?
        if (r[[i]]$fun[1]=="plots" || (r[[i]]$fun[1]=="wrapper"&!is.null(r[[i]]$output2))) {
          
          # special parameters?
          extraparams <- ""
          if(r[[i]]$fun[1]=="wrapper"){
            # no dynamic height for now
            height <- 3
            width <- 7
            extraparams <- sprintf(",fig.width=%f,fig.height=%f", width, height)
          } else if (r[[i]]$fun[2]=="stats"&r[[i]]$fun[3]=="pathway"&r[[i]]$fun[4]=="bar") {
            # dynamic height
            if(r[[i]]$output2$nr!=0){
              # set plot height
              height <- (62+(r[[i]]$output2$nr*r[[i]]$output2$npanrow*23.7143)+84)/2304*12 # manually curated using pixel measurements on example
              width <- 5+3*r[[i]]$output2$npancol
            } else{
              # if empty plot, set height to 3
              height <- 3
              width <- 8
            }
            extraparams <- sprintf(",fig.width=%f,fig.height=%f", width, height)
          } else if(length(r[[i]]$fun)>=3){
            if(r[[i]]$fun[2]=="box"&(r[[i]]$fun[3]=="scatter"|r[[i]]$fun[3]=="special")){
              # dynamic height
              if(!is.null(r[[i]]$output2)){
                # set plot height
                height <- (23+(r[[i]]$output2*160.65)+33)/2304*32 # manually curated using pixel measurements on example
                width <- 7
              } else{
                # if empty plot, set height to 3
                height <- 3
                width <- 7
              }
              extraparams <- sprintf(",fig.width=%f,fig.height=%f", width, height)
            }
          }
          
          
          # plot
          if (!use.plotly) {
            # use output2 for mt_plots_net
            if(r[[i]]$fun[1]=="wrapper" || r[[i]]$fun[2]=="net"){
              writechunk(glue::glue("r[[{i}]]$output2"), params = extraparams)
            }else{
              writechunk( glue::glue("r[[{i}]]$output"), params = extraparams)
            }
          } else {
            writechunk( glue::glue("
plotlist = r[[{i}]]$output %>% lapply(ggplotly)
htmltools::tagList(setNames(plotlist, NULL))
                             "), params=paste0('results="show"', extraparams))
          }
          # }
          
        }
        
        ## statistical result table?
        if (r[[i]]$fun[1]=="stats" || (r[[i]]$fun[1]=="wrapper"&!is.null(r[[i]]$output))) {
          # write warning if df too large
          if(nrow(r[[i]]$output$table) > 1000){
            out(glue::glue('WARNING: Large data frame ({nrow(r[[i]]$output$table)} rows). Displaying first
                           1000 rows.'))
          }
          # write out datatable
          writechunk(glue::glue('
# extract result table
df<-r[[{i}]]$output$table
# add feature names
rd <- rowData(D)
df <- cbind(name=as.data.frame(rd)$name[match(df$var, rownames(rd))], df) %>%
  dplyr::arrange(p.value)
# subset large data frames
if(nrow(df) > 1000) df <- df[1:1000, ]
# output
DT::datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 20, lengthMenu = c(10*(2^(0:3)), nrow(df)), autoWidth = TRUE, width = 1200, dom = "Bitlrp", buttons = c("copy", "csv", "excel", "pdf", "print")), class = "cell-border stripe", extensions = "Buttons")  %>% DT::formatStyle(columns = c(1:ncol(df)), fontSize = "80%", target= "row", lineHeight="80%")'),
                     params = "results='asis'")
          
        }
        
        # empty line as spacer
        out("")
        
        
      } else {
        
        # data or stats reporting step
        if(r[[i]]$fun[2]=="data" || r[[i]]$fun[2]=="stats"){
          ## header
          out(glue::glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))
          ## detailed arguments?
          if (output.calls) {
            L <- r[[i]]$args
            out("*Function arguments:*<br/>")
            out(names(L) %>% lapply(function(x){sprintf("%s=%s",x, toString(L[[x]]))}) %>% paste0(collapse = "<br/>"))
            out("")
          }
          ## log text
          out(glue::glue('*Log text:*<br/>{r[[i]]$logtxt}\n\n'))
        }
        
        
        # special reporting step
        if (r[[i]]$fun[2]=="heading") {
          # add extra heading
          out(glue::glue('{strrep("#",r[[i]]$output$lvl)} {r[[i]]$output$title}'))
          # result level is this heading +1
          lvl = r[[i]]$output$lvl + 1
        }
        
        # special reporting step
        if (r[[i]]$fun[2]=="text") {
          # enter text
          out(r[[i]]$output$text)
          out("")
        }
      }
    }
  }
  
  # Render the enrichment_results table (added here outside the loop)
  if (!is.null(enrichment_results)) {
    df_enrichment_results <- as.data.frame(enrichment_results)
    out(glue::glue("```{{r output_enrichment_results, results='asis'}}
# Render the enrichment_results table
DT::datatable(enrichment_results, rownames = FALSE)
```"))
  }
  
  # Add the Plotly plot using an iframe with the file path of the downloaded plot
  out(glue::glue("```{{=html}}
<iframe src='{plotly_file_path}' height='600px' width='1000px'></iframe>
```"))


  out(glue::glue("```{{=html}}
<iframe src='{mod3_file_path}' height='600px' width='1000px'></iframe>
```"))

#   out(glue::glue("```{{=html}}
# {plotly_file_path()}
# {mod3_file_path()}
# ```"))
  
  # clean up
  close(h)
  
}
