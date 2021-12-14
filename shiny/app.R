library(shiny)
library(readxl)
library(DT)
library(gprofiler2)
library(plotly)
library(heatmaply)
ui <- fluidPage(
    headerPanel('Split-ORF Results showcase'),
    tabsetPanel(
        tabPanel("Reports",
                 sidebarPanel(
                     selectInput('Mainselection', 'Choose the results you want to see:',
                                 list('Split-ORF main pipeline NMD', 'Split-ORF main pipeline RI',
                                      'Uniqueness determination NMD','Uniqueness determination RI',
                                      'Ribo-seq','Ribo-seq mouse')),
                     conditionalPanel(
                         "input.Mainselection == 'Ribo-seq'",
                         sliderInput('minreads', 'Choose the minimum number of overlapping reads:',
                                     min = 0,max = 5, value = 2)
                     ),
                     conditionalPanel(
                       "input.Mainselection == 'Ribo-seq mouse'",
                       sliderInput('minreadsMouse', 'Choose the minimum number of overlapping reads:',
                                   min = 0,max = 5, value = 2)
                     )
                 ),
                 mainPanel(
                     htmlOutput("report")
                 )
        ),
        tabPanel("Overview Tables",
                 sidebarPanel(width = 3,
                     selectInput('Mainselection2', 'Choose the results you want to see:',
                                 list('AML Maxquant', 'Overview NMD', 'Overview RI')),
                     conditionalPanel(
                         "input.Mainselection2 == 'Overview NMD'",
                         selectInput('NMDSheet', 'Choose the table you want to see:',
                                     list('Genes','Transcripts','ORFs'))
                     ),
                     conditionalPanel(
                         "input.Mainselection2 == 'Overview RI'",
                         selectInput('RISheet', 'Choose the table you want to see:',
                                     list('Genes','Transcripts','ORFs'))
                     ),
                     conditionalPanel(
                         "input.Mainselection2 == 'AML Maxquant'",
                         selectInput('amlSheet', 'Choose the sheet you want to see:',
                                     list('All Rows (exact)','Rows of interest (exact)','All Rows (flanking)','Rows of interest (flanking)'))
                     )
                 ),
                 mainPanel(
                     conditionalPanel(
                         "input.Mainselection2 == 'AML Maxquant'",
                         tagList(
                             DTOutput("amltable"),
                             checkboxInput("aml_sel", "select/deselect all"),
                             actionButton("amlbutton","show heatmap"),
                             plotlyOutput("amlmap"),
                            # downloadButton('downloadImage', 'Download Heatmap')
                         )
                     ),
                     conditionalPanel(
                         "input.Mainselection2 == 'Overview NMD'",
                         tagList(
                             DT::DTOutput("NMDtable"),
                             checkboxInput("dt_sel", "select/deselect all"),
                             actionButton("NMDbutton","show functional analysis"),
                             plotlyOutput("NMDfunction"),
                             plotOutput("NMDfunction2")
                         )
                     ),
                     conditionalPanel(
                         "input.Mainselection2 == 'Overview RI'",
                         tagList(
                             DT::DTOutput("RItable"),
                             checkboxInput("dti_sel", "select/deselect all"),
                             actionButton("RIbutton","show functional analysis"),
                             plotlyOutput("RIfunction"),
                             plotOutput("RIfunction2")
                         )
                     )
                 )
        ),
        tabPanel("UCSC Links",
                 shiny::fluidRow(
                   shinydashboard::box(title = "NMD human", "View the transcripts, ORFs, unique regions and ribo-seq data in the UCSC genome browser", 
                                       shiny::actionButton(inputId='ab1', label="Open UCSC", 
                                                           icon = icon("align-center"), 
                                                           onclick ="window.open('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.customText=https://raw.githubusercontent.com/JustinMurtagh29/Split-Orf_data/main/NMD.txt')")
                   )
                 ),
                 shiny::fluidRow(
                   shinydashboard::box(title = "RI human", "View the transcripts, ORFs, unique regions and ribo-seq data in the UCSC genome browser", 
                                       shiny::actionButton(inputId='ab1', label="Open UCSC", 
                                                           icon = icon("align-center"), 
                                                           onclick ="window.open('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.customText=https://raw.githubusercontent.com/JustinMurtagh29/Split-Orf_data/main/RI.txt')")
                   )
                 ),
                 shiny::fluidRow(
                   shinydashboard::box(title = "NMD mouse", "View the transcripts, ORFs, unique regions and ribo-seq data in the UCSC genome browser", 
                                       shiny::actionButton(inputId='ab1', label="Open UCSC", 
                                                           icon = icon("align-center"), 
                                                           onclick ="window.open('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.customText=https://raw.githubusercontent.com/JustinMurtagh29/Split-Orf_data/main/NMD_mouse.txt')")
                   )
                 ),
                 shiny::fluidRow(
                   shinydashboard::box(title = "RI mouse", "View the transcripts, ORFs, unique regions and ribo-seq data in the UCSC genome browser", 
                                       shiny::actionButton(inputId='ab1', label="Open UCSC", 
                                                           icon = icon("align-center"), 
                                                           onclick ="window.open('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.customText=https://raw.githubusercontent.com/JustinMurtagh29/Split-Orf_data/main/RI_mouse.txt')")
                   )
                 )
        )
    )
    
)

addResourcePath("tmpuser", getwd())
server <- function(input, output) {
    getPage<-function() {
        if(input$Mainselection == "Ribo-seq"){
            #rmarkdown::render(input = "./RiboSeqReportShiny.Rmd", 
            #                  output_file = "./RiboSeqReport.html", 
            #                  #params = list(args = c('G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-13.21.38_new_NMD/BOWTIE/','G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-15.10.13_new_RI/BOWTIE/',input$minreads)),
            #                  params = list(args = c('./data/NMD/BOWTIE','./data/RI/BOWTIE',input$minreads)),
            #                  envir = new.env(parent = globalenv())
            #)
            #return(includeHTML("./RiboSeqReport.html"))
            fname=paste0("./RiboSeqReport_min_", input$minreads, ".html")
            return(includeHTML(fname))
        }
        if(input$Mainselection == "Ribo-seq mouse"){
          #rmarkdown::render(input = "./MouseRiboSeqReport.Rmd", 
          #                  output_file = "./MouseRiboSeqReport.html", 
          #                  params = list(args = c('./data/mouse/NMD/BOWTIE','./data/mouse/RI/BOWTIE',input$minreadsMouse)),
          #                  envir = new.env(parent = globalenv())
          #)
          #return(includeHTML("./MouseRiboSeqReport.html"))
          fmname=paste0("./MouseRiboSeqReport_min_", input$minreadsMouse, ".html")
          return(includeHTML(fmname))
        }
        if(input$Mainselection == "Split-ORF main pipeline NMD"){
            return(includeHTML("./data/NMD/Split-ORF_Report.html"))
        }
        if(input$Mainselection == "Split-ORF main pipeline RI"){
            return(includeHTML("./data/RI/Split-ORF_Report.html"))
        }
        if(input$Mainselection == "Uniqueness determination NMD"){
            return(includeHTML("./data/NMD/Uniqueness_Report.html"))
        }
        if(input$Mainselection == "Uniqueness determination RI"){
            return(includeHTML("./data/RI/Uniqueness_Report.html"))
        }
    }
    getSheet <- function(){
        if(input$Mainselection2 == "AML Maxquant"){
          return(read_xlsx(path = "./data/AML/ORF_search_unique_Exact.xlsx", sheet=input$amlSheet)) 
          # return(read_xlsx(path = "G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/shiny/splitORFdatabase/data/AML/ORF_search_unique_Exact.xlsx", sheet=input$amlSheet))
        }
        if(input$Mainselection2 == "Overview NMD"){
          return(read_xlsx(path = "./data/NMD/NMD_Overview.xlsx", sheet=input$NMDSheet)) 
          # return(read_xlsx(path = "G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/shiny/splitORFdatabase/data/NMD/NMD_Overview.xlsx", sheet=input$NMDSheet))
        }
        if(input$Mainselection2 == "Overview RI"){
          return(read_xlsx(path = "./data/RI/RI_Overview.xlsx", sheet=input$RISheet)) 
          # return(read_xlsx(path = "G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/shiny/splitORFdatabase/data/RI/RI_Overview.xlsx", sheet=input$RISheet))
        }
    }
    plotmap <- function(){
        cols=c(2,3,5:40)
        amldata <- getSheet()[input$amltable_rows_selected,cols]
        amldata=as.data.frame(amldata)
        row.names(amldata)=paste0(amldata[,1], " : ",amldata[,2])
        amldata=amldata[,3:ncol(amldata)]
        amldata <- replace(amldata, is.na(amldata), 0)
        amldata=as.matrix(amldata)
        mapheight = length(row.names(amldata))*10+450
        heatmaply(amldata) %>% layout(width=1400, height= mapheight)
    }
    output$amltable <- renderDT(getSheet(), filter = "top",selection = list(target = 'row'))
    aml_proxy <- DT::dataTableProxy("amltable")
    observeEvent(input$aml_sel, {
        if (isTRUE(input$aml_sel)) {
            DT::selectRows(aml_proxy, input$amltable_rows_all)
        } else {
            DT::selectRows(aml_proxy, NULL)
        }
    })
    output$amlmap <- renderPlotly({
        req(input$amlbutton)
        isolate({
            plotmap()
        })
    })

    output$NMDtable <- renderDT(getSheet(), filter = "top",selection = list(target = 'row'))
    dt_proxy <- DT::dataTableProxy("NMDtable")
    observeEvent(input$dt_sel, {
        if (isTRUE(input$dt_sel)) {
            DT::selectRows(dt_proxy, input$NMDtable_rows_all)
        } else {
            DT::selectRows(dt_proxy, NULL)
        }
    })
    output$NMDfunction <- renderPlotly({
        req(input$NMDbutton)
        isolate({
            functionality <- gost(as.vector(t(getSheet()[input$NMDtable_rows_selected,2])))
            gostplot(functionality) %>% layout(
              margin = list(b = 100, l = 100) # to fully display the x and y axis labels
            )
        })
    })

    tableheight<- reactive(
          if (is.null(input$NMDtable_rows_selected)){
            200
          }
          else{
            length(input$NMDtable_rows_selected)*5.3+200
          }
    )
    output$NMDfunction2 <- renderPlot({
      req(input$NMDbutton)
      isolate({
        functionality <- gost(as.vector(t(getSheet()[input$NMDtable_rows_selected,2])))
        publish_gosttable(functionality,use_colors = TRUE,show_columns = c("source", "term_name", "term_size", "intersection_size"),ggplot = TRUE)
      })
    },height = tableheight)
    output$RItable <- renderDT(getSheet(), filter = "top",selection = list(target = 'row'))
    dti_proxy <- DT::dataTableProxy("RItable")
    observeEvent(input$dti_sel, {
        if (isTRUE(input$dti_sel)) {
            DT::selectRows(dti_proxy, input$RItable_rows_all)
        } else {
            DT::selectRows(dti_proxy, NULL)
        }
    })
    #output$selected_rows <- renderPrint(print(input$dt_rows_selected))
    output$RIfunction <- renderPlotly({
        req(input$RIbutton)
        isolate({
            rifunctionality <- gost(as.vector(t(getSheet()[input$RItable_rows_selected,2])))
            gostplot(rifunctionality)  %>% layout(
              margin = list(b = 100, l = 100) # to fully display the x and y axis labels
            )
        })
    })
    ritableheight<- reactive(
      if (is.null(input$RItable_rows_selected)){
        200
      }
      else{
        length(input$RItable_rows_selected)*9.5+200
      }
    )
    output$RIfunction2 <- renderPlot({
      req(input$RIbutton)
      isolate({
        functionality <- gost(as.vector(t(getSheet()[input$RItable_rows_selected,2])))
        publish_gosttable(functionality,use_colors = TRUE,show_columns = c("source", "term_name", "term_size", "intersection_size"),ggplot = TRUE)
      })
    },height = ritableheight, width = 1400)
    output$report <-renderUI(getPage())
    
}

shinyApp(ui = ui, server = server)