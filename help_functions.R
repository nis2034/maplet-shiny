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
mod3_plots_pls <- function (D, title = "PLS", 
                            subgroupvar,hover, ...) {
                            
 
 

 
  p <- NULL
  
  Dplot <- D
  # subgroupvar <- "GROUP_ID"

  
  pls_D <- autonomics:::pls(Dplot, subgroupvar,color =!!ensym(subgroupvar))
  
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
    
      if((replace(dt[[i]]$fun[2], is.na(dt[[i]]$fun[2]),"") =="volcano") || (replace(dt[[i]]$fun[3],is.na(dt[[i]]$fun[3]),"") =="scatter")||(replace(dt[[i]]$fun[3], is.na(dt[[i]]$fun[3]),"") =="pathway"))
      {
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
  
  
  #funargs <- mti_funargs()
  #D %<>%
  #  mti_generate_result(
  #    funargs = funargs,
  #    logtxt = sprintf('added pathway annotations using the %s pathway database', pwdb_name)
  #  )
  
  D
}


mt_stats_pathway_enrichment_new <- function(D, stat_name, pw_col, cutoff = 0.05) {
  
  stopifnot("SummarizedExperiment" %in% class(D))
  
  meta_D <- metadata(D)
  
  if(!"pathways" %in% names(meta_D)) stop("'pathways' does not exist in current SummarizedExperiment input")
  
  # Check if given pathway column actually exists
  if (!pw_col %in% names(meta_D$pathways)) stop(sprintf("'%s' not found in metabolite annotations.", pw))
  
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
                     n_pw_sig = sum(significant)) %>%
                     #mean_fc = mean(fc)) %>%
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
                     p_value_adjusted = stats::p.adjust(p_value, method = "fdr")) %>%
                     #mean_foldchange = mean_fc) %>%
    dplyr::arrange(p_value)
  
  metadata(D)$pathways$enrichment_results <-
    dplyr::as_tibble(enrichment_results)
  
  
  # funargs <- mti_funargs()
  # D %<>%
  #   mti_generate_result(
  #     funargs = funargs,
  #     logtxt = sprintf("performed pathway enrichment on %s pathways using Fihser's exact test",
  #                      nrow(enrichment_results))
  #   )
  
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