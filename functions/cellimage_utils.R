## Like filter_immunomodulator_expression_df but for multiple genes

multi_filter_immunomodulator_expression_df <- function(
  df, id_col, filter_col, expression_col, filter_values){
  
  df %>% 
    get_complete_df_by_columns(c(
      id_col, 
      filter_col, 
      expression_col)) %>% 
    dplyr::select(
      FILTER = filter_col, 
      COUNT = expression_col,
      ID = id_col) %>% 
    dplyr::filter(FILTER %in% filter_values) %>% 
    dplyr::mutate(LOG_COUNT = log10(COUNT + 1)) %>% 
    dplyr::select(ID, FILTER, LOG_COUNT)
}

## This function is based on build_immunomodulator_expression_df

build_multi_immunomodulator_expression_df <- function(
  group_df, ## fmx_df row filtered based on availability of sample groups
  genes_needed,  ## gene choices 
  group_col, ## the fmx_df column for the group
  expression_df = panimmune_data$im_expr_df,
  expression_filter_col = "Symbol",
  expression_col = "normalized_count",
  id_col = "ParticipantBarcode"){
  
  expression_df <- multi_filter_immunomodulator_expression_df(
    expression_df, 
    id_col, 
    expression_filter_col,
    expression_col,
    genes_needed)
  
  group_df <- group_df %>% 
    get_complete_df_by_columns(c(group_col, id_col)) %>% 
    select(GROUP = group_col, ID = id_col)
  
  result_df <- 
    dplyr::inner_join(group_df, expression_df, by = "ID") %>%
    dplyr::select(ID,GROUP, FILTER, LOG_COUNT)

}
