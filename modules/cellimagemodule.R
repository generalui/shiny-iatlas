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

      fluidRow(
        optionsBox(
          column(
            width = 6,
            uiOutput(ns("select_ui"))
          )
      ),
            
        plotBox(
          width = 8,
          plotOutput(ns("cellPlot"), height = 600) %>%
            shinycssloaders::withSpinner()
        )
      )
      
    )
    
  )
  
}

cellimage <- function(
    input, 
    output, 
    session, 
    group_display_choice, 
    group_internal_choice,
    study_subset_choice,
    sample_group_df,
    subset_df, 
    plot_colors
){
    
    ns <- session$ns

    output$select_ui <- renderUI({
        
        req(
            panimmune_data$sample_group_df,
            group_internal_choice()
        )
        
        if(group_internal_choice() == "Subtype_Curated_Malta_Noushmehr_et_al"){
            req(study_subset_choice())
        }
        
        sample_group_vector <-  panimmune_data$sample_group_df %>% 
            dplyr::filter(sample_group ==  group_internal_choice()) %>% 
            `if`(
                group_internal_choice() == "Subtype_Curated_Malta_Noushmehr_et_al",
                dplyr::filter(., `TCGA Studies`== study_subset_choice()),
                .
            ) %>%
            dplyr::pull(FeatureValue)
        
        selectInput(
            ns("tbd_method"),
            "Select Group",
            choices = sample_group_vector
        )

    })
    
    
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
    cat("my choice ",input$tbd_method,"\n")
    cat("group_display_choice",group_display_choice(),"\n")
    cat("group_internal_choice",group_internal_choice(),"\n")
    cat("group membs")
    plot(0,0)
  })
  
  
  
}