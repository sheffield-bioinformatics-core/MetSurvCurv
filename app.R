#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ggfortify)
library(partykit)
library(survival)
#data read in
expression <- read_csv("data/vdx_expression.csv")
meta <- read_csv("data/vdx_meta.csv")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("MetSurvCurv"),

    sidebarLayout(position = "left",
        
        sidebarPanel(
            strong("Select Dataset"),
            actionButton(inputId = "vdx", 
                         label = "VDX", 
                    
                
            ),
            selectInput(inputId = "genename", 
                    label ="Select Gene", 
                    choices =colnames(expression),
                    selected="ABO"
                    ),
            width=3
            ),
        mainPanel(
            plotOutput("survPlot"),
        )
    )
)

# Data pre-processing ----




# Define server logic required to draw a survivaL plot
server <- function(input, output) {
    
    curve_data = reactive({
        GName <- input$genename
        gene <- expression %>% select(ID,all_of(GName))
        surv_xfs <- Surv(meta$time/12, meta$event)
        curve_data <- meta %>% left_join(gene, by = "ID") %>% cbind(surv_xfs)
        fmula = as.formula(paste("surv_xfs ~ ", GName))
        ctree_xfs <- ctree(fmula, data = curve_data)
        splitrules <- partykit:::.list.rules.party(ctree_xfs, i=nodeids(ctree_xfs))[2] %>% str_split(pattern = " ")
        split <- as.numeric(splitrules[[1]][3])
        #signif <- 1 - ctree_xfs@tree$criterion$maxcriterion
        
        curve_data <- curve_data %>% mutate(threshold = get(GName)<=split)

        curve_data
    })
    observe({print(curve_data())})
    output$survPlot <- renderPlot({
            surv_xfs_dat <- survfit(surv_xfs~threshold,curve_data())
            autoplot(surv_xfs_dat) +
            labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                 title = input$genename) + 
            theme(plot.title = element_text(hjust = 0.5), 
                  axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                  axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                  legend.title = element_text(face="bold", size = 10))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
