
get.needed.variables <- function(unique.image.variable.ids,variable.annotations,source){
  variable.annotations %>% 
    filter(ImageVariableID %in% unique.image.variable.ids & Source==source) %>%
    pluck("FeatureLabel") %>% unique()
}

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

# cellcontent functions -------------------------------------------------------

build_cellcontent_df <- function(df, group_column){
  assert_df_has_columns(df, c(group_column, "Stromal_Fraction", "leukocyte_fraction"))
  long_df <- df %>% 
    dplyr::select(
      GROUP = group_column,
      "Stromal_Fraction", 
      "leukocyte_fraction") %>% 
    tidyr::drop_na()
  
  if(nrow(long_df) == 0) return(long_df)
  
  result_df <- long_df %>% 
    dplyr::mutate(Tumor_Fraction = 1 - Stromal_Fraction) %>% 
    tidyr::gather(fraction_type, fraction, -GROUP)
  assert_df_has_columns(result_df, c("GROUP", "fraction_type", "fraction"))
  return(result_df)
}

build_cell_fraction_df <- function(df, group_column, value_columns){
  assert_df_has_columns(df, c(group_column, value_columns))
  result_df <- df %>% 
    dplyr::select(GROUP = group_column, value_columns) %>% 
    tidyr::gather(fraction_type, fraction, -GROUP) %>% 
    tidyr::drop_na()
  assert_df_has_columns(result_df, c("GROUP", "fraction_type", "fraction"))
  return(result_df)
}
