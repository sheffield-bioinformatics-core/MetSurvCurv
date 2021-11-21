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
expression <- read_csv("data/vdx_expression.csv")
meta <- read_csv("data/vdx_meta.csv")

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
            plotOutput("survPlot"),
            textOutput("dataSelection")
    ),
    footer = dashboardFooter(left = "Sheffield Bioinformatics Core",right = "2021"
    )
)


# Data pre-processing ----

# Define server logic required to draw a survival plot
server <- function(input, output, session) {
    
    updateSelectizeInput(session, "genename", choices = colnames(expression)[-1], selected="GPN3",server = TRUE)
    
    curve_data = reactive({
        GName <- input$genename
        gene <- expression %>% select(ID,all_of(GName))
        surv_xfs <- Surv(meta$time/12, meta$event)
        curve_data <- meta %>% left_join(gene, by = "ID") %>% cbind(surv_xfs)
        fmula = as.formula(paste("surv_xfs ~ ", GName))
        ctree_xfs <- ctree(fmula, data = curve_data)
        splitrules <- partykit:::.list.rules.party(ctree_xfs, i=nodeids(ctree_xfs))[2] %>% str_split(pattern = " ")
        split <- round(as.numeric(splitrules[[1]][3]),3)
        curve_data$split <- split
        #signif <- 1 - ctree_xfs@tree$criterion$maxcriterion
        curve_data <- curve_data %>% mutate(threshold = ifelse(get(GName)<=split, paste(input$genename,intToUtf8(8804) , split), paste(input$genename,">", split)))
        curve_data
    })
    
    output$dataSelection <- renderText({
       paste("Gene expression datasets published by Wang et al. [2005] and Minn et al. [2007] (VDX)")
    })
    
    #observe({print(curve_data())})
    output$survPlot <- renderPlot({
            plot_data <- curve_data()
            surv_xfs_dat <- survfit(surv_xfs~threshold,plot_data)
            surv_diff_dat <- survdiff(surv_xfs~threshold,plot_data)
            freqs <- plot_data %>% count(threshold)
            newPval <- signif(pchisq(surv_diff_dat$chisq, df = length(surv_diff_dat$n)-1, lower.tail=FALSE),3)
            autoplot(surv_xfs_dat) +
                labs(x = "\n Time to Metastasis", y = "Probability of freedom from metastasis\n", title = paste(input$genename,"- p=", newPval)) +
                theme(plot.title = element_text(hjust = 0.5))+
                annotate(geom = "table", x = 0, y = 0.15, label = list(freqs), vjust = 1, hjust = 0) +
                scale_y_continuous(limits = c(0,1), labels = scales::percent)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
