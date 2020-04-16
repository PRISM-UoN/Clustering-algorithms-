#install.packages("survival")
#install.packages("survminer")
#install.packages("dplyr")

# Load packages
  library(survival)
  library(survminer)
  library(dplyr)
  library(foreign)
  library(haven)
  
# Load dataset
  # Set working directory
    setwd("R:/PRISM_DataSci/Stroke Medicine/ESUS Data")
    
  # Import data
    esus <- read_dta("esus_survival.dta")
    
    glimpse(esus)
    
    
  # Fit survival data using the Kaplan-Meier Method
    surv_object <- Surv(time = esus$time, event = esus$Recurrence)
    
    surv_object
    
    
  # Fit the Kaplan-Meier curve
    fit_pes_grp <- survfit(surv_object ~ pes_grp, data = esus)

    # Generate plot
      ggsurvplot(fit_pes_grp, data = esus, legend = "top",
                 legend.title = "Number of PES",
                 legend.labs = c("0 or 1 PES", "2 PES", "3 or more PES"),
                 risk.table = TRUE, risk.tab.y.text.col = TRUE,
                 break.time.by = 2,
                 xlim = c(0, 10), pval = TRUE)
      
      
      survminer::ggsurvplot(
        fit = survival::survfit(survival::Surv(time, Recurrence) ~ pes_grp, data = esus), 
        xlab = "Time (years)",
        ylab = "Overall survival probability",
        legend.title = "Number of PES",
        legend.labs = c("0 or 1 PES", "2 PES", "3 or more PES"),
        break.x.by = 2, xlim = c(0, 10),
        censor = FALSE,
        risk.table = TRUE,
        pval = TRUE,
        risk.table.y.text = FALSE)
      
      
  # Fit a Cox proportional hazards model
    fit.coxph <- coxph(surv_object ~ pes_grp,  data = esus)
    
    ggforest(fit.coxph, data = esus)
    
    
    
  # Kaplan-Meier plots stratified by clusters  
    survminer::ggsurvplot(
      fit = survival::survfit(survival::Surv(time, Recurrence) ~ cluster, data = esus), 
      xlab = "Time (years)",
      ylab = "Overall survival probability",
      legend.title = "Clusters",
      legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
      break.x.by = 2,
      break.y.by = 0.10,
      xlim = c(0, 10),
      ylim = c(0.2, 1.00),
      censor = FALSE,
      risk.table = TRUE,
      pval = TRUE,
      risk.table.y.text = FALSE)
      