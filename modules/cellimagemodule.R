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
    
    
#  data_df <- reactive({
#    subset_df() %>% 
#      dplyr::select(
#        x = group_internal_choice(), 
#        "ParticipantBarcode") %>% 
#      dplyr::inner_join(panimmune_data$im_expr_df, by = "ParticipantBarcode") %>% 
#      dplyr::rename(label = "ParticipantBarcode")
#  })
  
  
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

    ## Annotations of image objects
    ## Variable annotations are ImageVariableID, FeatureLabel, Source, ColorScale
    variable.annotations <- readr::read_tsv('data/cell_image_id_annotations.tsv') 
    ## Image obects, in order, labeled in terms of ImageVariableID
    image.object.labels <- read.table('data/cell_image_object_ids.txt',as.is=T)$V1
    ## must be met 
    ## image.object.labels %in% variable.annotations$ImageVariableID
    unique.image.variable.ids <- unique(image.object.labels)
    
    ##
    ## Needed cellular data
    ##
    
    
    cois <- get.data.variables(unique.image.variable.ids,variable.annotations,'fmx_df')
    dfc <- build_cellcontent_df(subset_df(),cois,group_internal_choice())
#    dfc <- build_cellcontent_df(group_df,cois,group_col) 
    dfc <- dfc %>% dplyr::rename(Group=GROUP,Variable=fraction_type,Value=fraction)
    ## Note that ParticipantBarcode is gone.  Each Group,Variable combo simply has instances
    
    
    ##
    ## Needed gene expression data
    ##
    
    ## input unique image variable IDs, get genes with IDs as in expression matrix
    gois <- get.data.variables(unique.image.variable.ids,variable.annotations,'im_expr_df')
    dfg <- build_multi_imageprotein_expression_df(subset_df(),gois,group_internal_choice())  ## dfg$FILTER is the Gene column 
    dfg <- dfg %>% dplyr::select(Group=GROUP,Variable=FILTER,Value=LOG_COUNT)
    ## Note that "ID" aka ParticipantBarcode is gone.  Each Group,Variable combo simply has instances
    
    ### data frame of all values
    
    ## dfc$Variable is chr, dfg$Variable is factor
    dfv <- dplyr::bind_rows(dfc, dfg)
    
    #########################################################################
    ##
    ## Variables ranges and summary
    ##
    #########################################################################
    
    ## Mean Value per Group and Variable
    meanz <- dfv %>% dplyr::group_by(Group,Variable) %>% dplyr::summarize(Mean=mean(Value)) 
    ## Max Value for each Variable (includes avg over Group)
    maxz <- dfv %>% dplyr::group_by(Variable) %>% dplyr::summarize(Max=max(Value))
    ## Min Value for each Variable (includes avg over Group)
    minz <- dfv %>% dplyr::group_by(Variable) %>% dplyr::summarize(Min=min(Value)) 
    ## Vector versions
    minvec <- minz %>% purrr::pluck("Min") ; names(minvec) <- minz %>% purrr::pluck("Variable")
    maxvec <- maxz %>% purrr::pluck("Max") ; names(maxvec) <- maxz %>% purrr::pluck("Variable")
    
    #########################################################################
    ##
    ## Get image and convert to grid object
    ##
    #########################################################################
    
    pic <- grImport2::readPicture("data/tcell-svg-take3-cairo.svg")
    #pic <- grImport2::readPicture("data/tcell-start-cairo-edited.svg")
    w <- grImport2::pictureGrob(pic)
    gTree.name <- grid::childNames(w) ## label of overall gTree object
    pathlabels <- w$children[[gTree.name]]$childrenOrder ## labels and order of children 
    fill.color.start <- character(length(pathlabels)) ; names(fill.color.start) <- pathlabels
    for (s in pathlabels){
      fill.color.start[s] <- w$children[[gTree.name]]$children[[s]]$gp$fill 
    }
    wnew <- w
    fill.color.new <- character(length(pathlabels)) ; names(fill.color.new) <- pathlabels ## this is for editing
    
    
    #########################################################################
    ##
    ## Get New Colors
    ##
    #########################################################################
    
    #sois <- unique(data_df[[group_col]])
    #soi <- "BRCA.LumB"
    #soi <- sois[5]
    soi <- input$tbd_method
    
    for (ind in seq(1,length(image.object.labels))){
      ioa <- image.object.labels[ind]
      datavar <- variable.annotations %>% dplyr::filter(ImageVariableID==ioa) %>% purrr::pluck("FeatureLabel")
      colormap <-   variable.annotations %>% dplyr::filter(ImageVariableID==ioa) %>% purrr::pluck("ColorScale")
      fill.color.new[ind] <- getVarColor(datavar,soi,colormap,dfv,minvec,maxvec)
    }
    for (s in pathlabels ){
      wnew$children[[gTree.name]]$children[[s]]$gp$fill <- fill.color.new[s]
    }
    
    #########################################################################
    ##
    ## DRAW 
    ##
    #########################################################################
    
    grid::grid.draw(wnew)
    
  })
  
  
  
}