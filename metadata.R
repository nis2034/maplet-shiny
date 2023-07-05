library(maplet)

mt_anno_pathways_hmdb_ch <- function(D,
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
  
  print(colnames(pwdb))
  
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




file_data <- system.file("extdata", "example_data/sampledata.xlsx", package = "maplet")

# Preprocess dataset ------------------------------------------------------
D_pre <-
  mt_load_metabolon_v1(file=file_data, sheet="OrigScale") %>%
  mt_anno_pathways_hmdb_ch(in_col = "HMDb_ID",
                        out_col = "kegg_db",
                        pwdb_name = "KEGG",
                        db_dir = system.file("extdata", "precalc/hmdb", package = "maplet")) %>%
  mt_anno_pathways_remove_redundant(feat_col = "HMDb_ID", pw_col = "kegg_db") %>%
  mt_pre_filter_missingness(feat_max=0.2) %>%
  mt_pre_filter_missingness(samp_max=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batch_col = "BATCH_MOCK") %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # logging
  mt_pre_trans_log() %>%
  # KNN imputation
  mt_pre_impute_knn() %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group,
    samp_filter = (Group %in% c("treatment1","treatment2")),
    stat_name         = "comp",
    n_cores     = 1
  ) %>%
  # add fold changes to result tables
  mt_post_fold_change(stat_name = "comp") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "comp", method = "BH")




# Apply mt_stats_pathway_enrichment ---------------------------------------

D_pw <-
  D_pre %>%
  mt_stats_pathway_enrichment(pw_col = "kegg_db",
                              stat_name = "comp",
                              cutoff = 0.4)


# Show results ---------------------------------------
metadata(D_pw)$pathways$enrichment_results