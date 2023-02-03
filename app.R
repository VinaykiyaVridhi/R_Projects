## Author: Vridhi Vinaykiya
## vridhiv@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
suppressPackageStartupMessages(library(tidyverse))
library(shiny)
library(ggplot2)
library(colourpicker)

# Define UI for application that draws a volcano plot
ui <- fluidPage(
    titlePanel("BF591 Assignment 7"),
    p("To use this application, download the CSV deseq_res.csv from the data directory of this app's repository."),
    sidebarLayout(
        sidebarPanel(
            fileInput("file", "Load differential expression results", accept=".csv"),
            radioButtons("xvar", "Choose the column for the x-axis",
                         choices = c("baseMean","log2FoldChange",
                                     "lfcSE", "stat", "pvalue", "padj"),
                         selected = "log2FoldChange"),
            radioButtons("yvar", "Choose the column for the y-axis",
                         choices = c("baseMean","log2FoldChange",
                                     "lfcSE", "stat", "pvalue", "padj"),
                         selected = "padj"),
            colourInput("col1", "Base point color", "#22577A"),
            colourInput("col2", "Highlight point color", "#FFCF56"),
            sliderInput(inputId = "slider", min = -300, max = 0,
                        label = "Select the magnitude of the p adjusted coloring:", 
                        value = -150, step = -30),
            actionButton("submit","Plot", icon = icon("chart-line"), width = "100%")
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Plot",
                         plotOutput("volcano")
                ),
                tabPanel("Table",
                         tableOutput("table")
                )
            )
        )
    )
)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is
    #' still a "function", but it will take no arguments. The `reactive({})` bit
    #' says "if any of my inputs (as in, input$...) are changed, run me again".
    #' This is useful when a user clicks a new button or loads a new file. In
    #' our case, look for the uploaded file's datapath argument and load it with
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
        if(is.null(input $ file $ datapath)) 
            {
            return(NULL)
        }
        dataf <- read.csv(input $ file $ datapath, header = TRUE)
        names(dataf)[1] <- "gene"
        return(dataf)
    })
    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles.
    #' Write a normal volcano plot using geom_point, and integrate all the above
    #' values into it as shown in the example app. The testing script will treat
    #' this as a normal function.
    #'
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
            if(is.null(dataf))
                return(NULL)
            plot <- ggplot(data=dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), 
                                           color = !!sym(y_name)<(1*10^(slider))))+
                geom_point() +
                theme_minimal() +
                scale_color_manual(values=c(color1, color2)) +
                ylab("-log10(padj)")+
                xlab(x_name)
            plot $ labels $ colour <- paste0(y_name," < 1 x 10^",toString(slider))
            return(plot)
        }
    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will
    #' evaluate it normally. Not only does this function filter the data frame to
    #' rows that are above the slider magnitude, it should also change the format
    #' of the p-value columns to display more digits. This is so that it looks
    #' better when displayed on the web page. I would suggest the function
    #' `formatC()`
    #'
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
        if(is.null(dataf))
            return(NULL)
        dataf <- dataf[dataf $ padj < (1*10^slider),]
        dataf $ pvalue <- formatC(dataf $ pvalue,
                                digits = 10)
        dataf $ padj <- formatC(dataf $ pad,
                              digits = 10)
        return(dataf)
    }
    #' These outputs aren't really functions, so they don't get a full skeleton,
    #' but use the renderPlot() and renderTabel() functions to return() a plot
    #' or table object, and those will be displayed in your application.
    observeEvent(input$submit,{
        volcano <- volcano_plot(load_data(), input $ xvar, input $ yvar,
                                input $ slider, input $ col1, input $ col2)
        dataf <- na.omit(load_data())
        table <- draw_table(dataf, input $ slider)
        output $ volcano <- renderPlot(volcano)
        output $ table <-  renderTable(table)
    })
}
# Run the application
shinyApp(ui = ui, server = server)