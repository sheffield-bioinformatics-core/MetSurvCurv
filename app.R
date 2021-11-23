# This is MetSurvCurv
# Emily V Chambers, Mark J Dunning, Sheffield Bioinformatics Core

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggfortify)
library(partykit)
library(survival)
library(ggpmisc)
library(coin)

#data read in
expression <- read_csv("data/vdx_expression.csv",col_types = cols())
meta <- read_csv("data/vdx_meta.csv",col_types = cols())

ui <- dashboardPage(
    skin = "purple",

    header = dashboardHeader(title="MetSurvCurv"
            ),
    sidebar = dashboardSidebar(
            collapsed=F,
            selectizeInput(inputId = "genename", 
                    label ="Select Gene:", 
                    choices =NULL
                    )
            ),
    body = dashboardBody(
      fluidRow(
        box(
          width=8,
          plotOutput("survPlot")
        )
      ),
      fluidRow(
        box(
          width=8,
          textOutput("dataSelection")
        )
      ),
      tags$head(tags$style('.box-header{ display: none}')) 
    ),
    footer = dashboardFooter(left = "Sheffield Bioinformatics Core",right = "2021"
    )
)


# Data pre-processing ----

# Define server logic required to draw a survival plot
server <- function(input, output, session) {
    updateSelectizeInput(session, "genename", choices = colnames(expression)[-1], selected="GPN3",server = TRUE)
  
  output$dataSelection <- renderText({
    paste("Gene expression datasets published by Wang et al. [2005] and Minn et al. [2007] (VDX)")
  })
  
    output$survPlot <- renderPlot({
      req(input$genename)
      GName <- input$genename
      gene <- expression %>% select(ID,all_of(GName))
      surv_xfs <- Surv(meta$time/12, meta$event)
      curve_data <- meta %>% left_join(gene, by = "ID") %>% cbind(surv_xfs)
      fmula = as.formula(paste("surv_xfs ~ ", GName))
      ctree_xfs <- ctree(fmula, data = curve_data)
      splitrules <- partykit:::.list.rules.party(ctree_xfs, i=nodeids(ctree_xfs)) 
      if(length(splitrules)>1){
        split <-str_split(splitrules[2],pattern = " ")[[1]][3] %>% as.numeric() 
      }else{
        split <- median(curve_data[,GName])
      }
        curve_data$split <- split
        curve_data <- curve_data %>% mutate(threshold = ifelse(get(GName)<=split, paste(GName,intToUtf8(8804) , signif(split,3)), paste(GName,">", signif(split,3))))
        surv_xfs_dat <- survfit(surv_xfs~threshold,curve_data)
        surv_diff_dat <- survdiff(surv_xfs~threshold,curve_data)
        freqs <- curve_data %>% count(threshold)
        newPval <- signif(pchisq(surv_diff_dat$chisq, df = length(surv_diff_dat$n)-1, lower.tail=FALSE),3)
        p <- autoplot(surv_xfs_dat) +
              labs(x = "\n Time to Metastasis", y = "Probability of freedom from metastasis\n", title = paste(GName,"- p=", newPval)) +
              theme(plot.title = element_text(hjust = 0.5))+
              annotate(geom = "table", x = 0, y = 0.15, label = list(freqs), vjust = 1, hjust = 0) +
              scale_y_continuous(limits = c(0,1), labels = scales::percent)
        return(p)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
