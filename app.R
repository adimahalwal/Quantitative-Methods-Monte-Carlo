# This is a Shiny web application. You can run the application by clicking 
# the 'Run App' button above.

library(shiny)
library(tidyverse)
library(MASS)

# Define UI
ui <- fluidPage(
  titlePanel("Monte Carlo Simulation of OLS Estimates with Multiple Parameter Sets (procedure 3)"),
  sidebarLayout(
    sidebarPanel(
      h4("Simulation Settings"),
      sliderInput("n", "Sample Size (n):", min = 10, max = 500, value = 50),
      sliderInput("num_sim", "Number of Simulations:", min = 100, max = 5000, value = 1000),
      
      h4("Parameter Set 1"),
      sliderInput("beta_0_1", "True Intercept (beta_0):", min = -10, max = 10, value = 2),
      sliderInput("beta_1_1", "True Slope (beta_1):", min = -10, max = 10, value = 1),
      sliderInput("beta_2_1", "True Slope (beta_2):", min = -10, max = 10, value = 1),
      sliderInput("rho_1", "Correlation (rho):", min = -1, max = 1, value = 0.1, step = 0.1),
      sliderInput("sigma_1", "Error Standard Deviation (sigma):", min = 0, max = 20, value = 10),
      
      h4("Parameter Set 2"),
      sliderInput("beta_0_2", "True Intercept (beta_0):", min = -10, max = 10, value = 3),
      sliderInput("beta_1_2", "True Slope (beta_1):", min = -10, max = 10, value = 1.5),
      sliderInput("beta_2_2", "True Slope (beta_2):", min = -10, max = 10, value = 0.5),
      sliderInput("rho_2", "Correlation (rho):", min = -1, max = 1, value = 0.2, step = 0.1),
      sliderInput("sigma_2", "Error Standard Deviation (sigma):", min = 0, max = 20, value = 10)
    ),
    mainPanel(
      verbatimTextOutput("results"),
      plotOutput("density_plots")
    )
  )
)

# Define server
server <- function(input, output) {
  
  # Function to run simulation based on parameter set inputs
  run_simulation <- function(beta_0, beta_1, beta_2, rho, sigma) {
    mean_vector <- c(0, 0)
    var_covar <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
    
    # Check if var_covar is positive definite
    if (det(var_covar) <= 0) {
      stop("The variance-covariance matrix is not positive definite.")
    }
    
    # Containers for results
    beta_0_estimates <- numeric(input$num_sim)
    beta_1_estimates <- numeric(input$num_sim)
    beta_2_estimates <- numeric(input$num_sim)
    r_squared_values <- numeric(input$num_sim)
    t_statistics_beta1 <- numeric(input$num_sim)
    p_values_beta1 <- numeric(input$num_sim)
    
    last_model <- NULL  # Placeholder for the last model
    
    # Monte Carlo Simulation
    for (i in 1:input$num_sim) {
      x <- mvrnorm(input$n, mu = mean_vector, Sigma = var_covar)
      epsilon <- rnorm(input$n, mean = 0, sd = sigma)
      y <- beta_0 + beta_1 * x[,1] + beta_2 * x[,2] + epsilon
      
      model <- lm(y ~ x[,1] + x[,2])
      model_summary <- summary(model)
      p_value <- model_summary$coefficients[3, "Pr(>|t|)"]
      
      # Check that p_value is a single value
      if (length(p_value) == 1 && p_value < 0.05) {
        beta_0_estimates[i] <- coef(model)[1]
        beta_1_estimates[i] <- coef(model)[2]
        beta_2_estimates[i] <- coef(model)[3]
        r_squared_values[i] <- model_summary$r.squared
        t_statistics_beta1[i] <- model_summary$coefficients[2, "t value"]
        p_values_beta1[i] <- model_summary$coefficients[2, "Pr(>|t|)"]
      } else {
        model <- lm(y ~ x[, 1])  # Adjust model to exclude x[,2]
        model_summary <- summary(model)
        beta_0_estimates[i] <- coef(model)[1]
        beta_1_estimates[i] <- coef(model)[2]
        r_squared_values[i] <- model_summary$r.squared
        t_statistics_beta1[i] <- model_summary$coefficients[2, "t value"]
        p_values_beta1[i] <- model_summary$coefficients[2, "Pr(>|t|)"]
      }
      
      last_model <- model  # Store the last model
    }
    
    list(
      beta_0_estimates = beta_0_estimates,
      beta_1_estimates = beta_1_estimates,
      beta_2_estimates = beta_2_estimates,
      r_squared_values = r_squared_values,
      t_statistics_beta1 = t_statistics_beta1,
      p_values_beta1 = p_values_beta1,
      rejection_rate_beta1 = mean(p_values_beta1 < 0.05),
      mean_beta_0 = mean(beta_0_estimates),
      mean_beta_1 = mean(beta_1_estimates),
      mean_beta_2 = mean(beta_2_estimates),
      mean_r_squared = mean(r_squared_values),
      sd_beta_0 = sd(beta_0_estimates),
      sd_beta_1 = sd(beta_1_estimates),
      sd_beta_2 = sd(beta_2_estimates),
      sd_r_squared = sd(r_squared_values),
      last_model = last_model  # Return the last model
    )
  }
  
  # Generate simulation results for each parameter set
  simulation_results <- reactive({
    list(
      set1 = run_simulation(input$beta_0_1, input$beta_1_1, input$beta_2_1, input$rho_1, input$sigma_1),
      set2 = run_simulation(input$beta_0_2, input$beta_1_2, input$beta_2_2, input$rho_2, input$sigma_2)
    )
  })
  
  # Display summary results
  output$results <- renderPrint({
    results <- simulation_results()
    
    cat("Parameter Set 1\n")
    cat("True Intercept (beta_0):", input$beta_0_1, "\n")
    cat("Mean of Estimated Intercepts:", results$set1$mean_beta_0, "\n")
    cat("Standard Deviation of Estimated Intercepts:", results$set1$sd_beta_0, "\n")
    cat("Mean of Estimated Slopes (beta_1, beta_2):", results$set1$mean_beta_1, ", ", results$set1$mean_beta_2, "\n")
    cat("Standard Deviation of Estimated Slopes (beta_1, beta_2):", results$set1$sd_beta_1, ", ", results$set1$sd_beta_2, "\n")
    cat("Mean R^2:", results$set1$mean_r_squared, "\n")
    cat("Standard Deviation of R^2:", results$set1$sd_r_squared, "\n")
    cat("Probability of Rejecting Null Hypothesis for beta_1:", results$set1$rejection_rate_beta1, "\n\n")
    
    cat("Parameter Set 2\n")
    cat("True Intercept (beta_0):", input$beta_0_2, "\n")
    cat("Mean of Estimated Intercepts:", results$set2$mean_beta_0, "\n")
    cat("Standard Deviation of Estimated Intercepts:", results$set2$sd_beta_0, "\n")
    cat("Mean of Estimated Slopes (beta_1, beta_2):", results$set2$mean_beta_1, ", ", results$set2$mean_beta_2, "\n")
    cat("Standard Deviation of Estimated Slopes (beta_1, beta_2):", results$set2$sd_beta_1, ", ", results$set2$sd_beta_2, "\n")
    cat("Mean R^2:", results$set2$mean_r_squared, "\n")
    cat("Standard Deviation of R^2:", results$set2$sd_r_squared, "\n")
    cat("Probability of Rejecting Null Hypothesis for beta_1:", results$set2$rejection_rate_beta1, "\n\n")
    
    cat("Summary of the last model from simulation for Parameter Set 1:\n")
    print(summary(results$set1$last_model))
    
    cat("\nSummary of the last model from simulation for Parameter Set 2:\n")
    print(summary(results$set2$last_model))
  })
  
  # Plot density plots with different curves for each parameter set
  output$density_plots <- renderPlot({
    results <- simulation_results()
    
    par(mfrow = c(1, 4), oma = c(0, 0, 2, 0), mar = c(5, 4, 4, 2))
    
    # Density plot for rejection rate of beta_1
    barplot(c(results$set1$rejection_rate_beta1, results$set2$rejection_rate_beta1),
            names.arg = c("Set 1", "Set 2"), col = c("blue", "red"),
            main = "Probability of Rejecting Null for beta_1", ylab = "Rejection Rate", ylim = c(0, 1))
    
    # Density plot for R^2 values
    plot(density(results$set1$r_squared_values), col = "blue", lwd = 2,
         main = "Density of R^2", xlab = "R^2")
    lines(density(results$set2$r_squared_values), col = "red", lwd = 2)
    legend("topright", legend = c("Set 1", "Set 2"), col = c("blue", "red"), lwd = 2)
    
    # Density plot for beta_1 estimates
    plot(density(results$set1$beta_1_estimates), col = "blue", lwd = 2,
         main = "Density of Slope Estimates (beta_1)", xlab = "Estimated beta_1")
    lines(density(results$set2$beta_1_estimates), col = "red", lwd = 2)
    legend("topright", legend = c("Set 1", "Set 2"), col = c("blue", "red"), lwd = 2)
    
    # Density plot for beta_2 estimates
    plot(density(results$set1$beta_2_estimates), col = "blue", lwd = 2,
         main = "Density of Slope Estimates (beta_2)", xlab = "Estimated beta_2")
    lines(density(results$set2$beta_2_estimates), col = "red", lwd = 2)
    legend("topright", legend = c("Set 1", "Set 2"), col = c("blue", "red"), lwd = 2)
  }, width = "auto", height = "auto")
}

# Run the application
shinyApp(ui = ui, server = server)

