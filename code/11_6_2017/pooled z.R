library(shiny)
library(MASS)
library(mvtnorm)


ui <- fluidPage(
                titlePanel("Distributions of pooled z values under different correlation"),
                sidebarLayout(
                  sidebarPanel(
                    
                    # select number of variables m 
                    sliderInput("m", 
                                "number of variables",
                                min = 10, max = 1000,
                                value = 100, 
                                step = 10),
                    
                    # select number of obs n 
                    sliderInput("n", 
                                "number of obervations",
                                min = 10, max = 300,
                                value = 100, 
                                step = 5),
                    
                    # select corrletaion 
                    sliderInput("rho", 
                                "pair wise correlation",
                                min = 0, max = 1,
                                value = 0.5, 
                                step = 0.05)
                  ),
                  
                  # Output: histogram of pooled z
                  
                  mainPanel(
                    plotOutput(outputId = "dist.z", height = "300px")
                   # textOutput(outputId = "desc")
                    
                  )
                )
)




server <- function(input, output) {
  
  # Subset data
  df <- reactive({
    m <- input$m
    n <- input$n
    rho <- input$rho
    
    sigma <- diag(m)+matrix(rep(rho,m*m),nrow=m)-diag(rep(rho,m))
    
    df <- data.frame(mvrnorm(n = n, mu=rep(0,m), sigma))
    
    return(df)
    
  })
  
  
  # Create scatterplot object the plotOutput function is expecting
  output$dist.z <- renderPlot({
    df <- df()
    pooled.z <- unlist(df)
    hist(pooled.z,prob=TRUE,breaks=50)
    curve(dnorm(x, mean=0, sd=1), 
          col="darkblue", lwd=2, add=TRUE, yaxt="n")
    
    
  })
  
  # # Pull in description of trend
  # output$desc <- renderText({
  #   trend_text <- filter(trend_description, type == input$type) %>% pull(text)
  #   paste(trend_text, "The index is set to 1.0 on January 1, 2004 and is calculated only for US search traffic.")
  # })
}

# Create Shiny object
shinyApp(ui = ui, server = server)


##################################################################
