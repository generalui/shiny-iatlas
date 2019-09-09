cellimage_UI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    titleBox("iAtlas Explorer — Cellular Image"),
    textBox(
      width = 12,
      p(stringr::str_c(
        "Explore the level of checkpoint proteins and cells",
        sep = " "
      ))  
    ),
    
    sectionBox(
      title = "Here we go",
      
      messageBox(
        width = 12,
        p("Here is what you must do.")
      ),
      
#      fluidRow(
#        optionsBox(
#          width = 4,
#          uiOutput(ns("survplot_opts")),
          
#          selectInput(
#            ns("timevar"),
#            "Survival Endpoint",
#            c("Overall Survival" = "OS_time", "Progression Free Interval" = "PFI_time_1"),
#            selected = "OS_time"
#          ),
          
#     ),
        
        plotBox(
          width = 8,
          plotOutput(ns("cellPlot"), height = 600) %>%
            shinycssloaders::withSpinner()
        )
 #     )
      
    )
    
  )
  
}

cellimage <- function(
    input, 
    output, 
    session, 
    group_display_choice, 
    group_internal_choice, 
    sample_group_df,
    subset_df, 
    plot_colors
){
  
  data_df <- reactive({
    subset_df() %>% 
      dplyr::select(
        x = group_internal_choice(), 
        "ParticipantBarcode") %>% 
      dplyr::inner_join(panimmune_data$im_expr_df, by = "ParticipantBarcode") %>% 
      dplyr::rename(label = "ParticipantBarcode")
  })
  
  
#  output$survplot_opts <- renderUI({
#    group_choice <- magrittr::set_names(list(group_internal_choice()), ss_choice())
#    var_choices <- c(
#      list("Current Sample Groups" = group_choice),
#      get_feature_df_nested_list())
#    selectInput(
#      ns("var1_surv"),
#      "Variable",
#      var_choices,
#      selected = group_internal_choice()
#    )
#  })
  
  
  output$cellPlot <- renderPlot({
    plot(0,0)
  })
  
}