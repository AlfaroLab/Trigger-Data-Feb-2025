 # =============================================================================
# Patches Variation Analysis 
# This code analyzes Trigger Fish patch spectral data to determine:
# 1. Whether we can tell different patch types apart (black, white, orange, blue)
# 2. How much variation exists within patches vs. across individuals
# =============================================================================

df <- read.csv("pablo_equation_1_transposed.csv")

# =============================================================================
# Defining the spectral angle calculation function
# =============================================================================
# This function calculates the angular difference between two spectral signatures
# Based on the formula from the piranha paper: α = cos⁻¹(ΣXY/√(Σ(X)²Σ(Y)²))
# Smaller angles indicate greater similarity between spectra
calc_spectral_angle <- function(x, y) 
{
  sum_xy <- sum(x * y)     
  sum_x2 <- sum(x^2)       
  sum_y2 <- sum(y^2)       
  alfa <- acos(sum_xy / sqrt(sum_x2 * sum_y2))
  
  return(alfa)  # Return angle in radians
}

# =============================================================================
# Identifying patch types and individuals for analysis
# =============================================================================
# Defining our patch types (colors) and individual identifiers (repeats)
patch_types <- c("black", "white", "orange", "blue")  # The four patch types
individuals <- c("01", "02", "03", "04")              # The four individuals (repeats)

# =============================================================================
# Organizing the data by patch type and individual
# =============================================================================
# Creating a nested list structure to easily access data:
# patch_columns$black$01 will contain the spectral data for black patch on individual 01
patch_columns <- list()
for(patch in patch_types) 
{
  # Creating a sub-list for each patch type
  patch_columns[[patch]] <- list()
  
  # For each individual, finding and storing the corresponding column data
  for(ind in individuals) 
  {
    col_name <- paste0(patch, "_", ind)  # Constructing column name (e.g., "black_01")
    
    # Only adding data if the column exists in our dataset
    if(col_name %in% names(df)) 
    {
      patch_columns[[patch]][[ind]] <- df[[col_name]]  # Storing the spectral data
    }
  }
}

# =============================================================================
# Calculating mean spectra for each patch type
# =============================================================================
patch_means <- data.frame(wavelength = df$wavelength)  # Starting with wavelength column
for(patch in patch_types) 
{
  # Combining all measurements for this patch type into a matrix
  patch_data <- do.call(cbind, patch_columns[[patch]])
  
  # Calculating row means (mean reflectance at each wavelength)
  patch_means[[patch]] <- rowMeans(patch_data, na.rm = TRUE)
}

# =============================================================================
# PLOT 1: Mean spectra for each patch type
# =============================================================================
# Defining colors for plotting that match the patch types
colors <- c("black", "gray", "orange", "blue")

# Starting with the first patch (black)
plot(df$wavelength, patch_means$black, 
     type="l",                                       # Line plot
     col=colors[1],                                  # black color
     lwd=2,                                          
     ylim=c(0, max(patch_means[,-1], na.rm=TRUE) * 1.1),  # Y-axis range with 10% margin
     xlab="Wavelength (nm)",                         
     ylab="Reflectance (%)",                         
     main="Mean Spectra by Patch Type")              

# Adding lines for the remaining patch types
for(i in 2:length(patch_types)) 
{
  lines(df$wavelength,                     # X values (wavelength)
        patch_means[[patch_types[i]]],     # Y values (reflectance)
        col=colors[i],                     # Color corresponding to patch type
        lwd=2)                             
}

# Legend to identify each line
legend("topright",                         
       patch_types,                        # Labels for legend
       col=colors,                         # Colors for legend
       lwd=2)                             

# =============================================================================
# Calculating spectral angles between different patch types
# =============================================================================
# This quantifies how different each patch type is from the others
cat("Spectral angles between patch types (in radians):\n")

# Creating a matrix to store the angles between each pair of patch types
patch_angles <- matrix(0,                          # Initializing with zeros
                       nrow=length(patch_types),   # Rows = number of patch types
                       ncol=length(patch_types))   # Columns = number of patch types

# Adding row and column names to the matrix for clarity
rownames(patch_angles) <- patch_types
colnames(patch_angles) <- patch_types

# Calculating the spectral angle between each pair of patch types
for(i in 1:length(patch_types)) 
{
  for(j in 1:length(patch_types)) 
  {
    if(i != j) # Skipping scenario comparing a patch with itself (angle would be 0)
    {  
      # Calculating spectral angle between mean spectra of two patch types
      patch_angles[i,j] <- calc_spectral_angle(
        patch_means[[patch_types[i]]],     # Mean spectrum of first patch type
        patch_means[[patch_types[j]]]      # Mean spectrum of second patch type
      )
    }
  }
}

# Printing the matrix of angles (rounded to 4 decimal places)
print(round(patch_angles, 4))

# =============================================================================
# Calculating within-patch variation
# =============================================================================
# This measures how much variation exists within the same patch type across different individuals
within_variation <- list()
for(patch in patch_types) 
{
  within_variation[[patch]] <- c()  # Initializing empty vector for this patch
  
  # Comparing all possible pairs of individuals for this patch type
  for(i in 1:(length(individuals)-1)) 
  {
    for(j in (i+1):length(individuals)) 
    {
      ind1 <- individuals[i]
      ind2 <- individuals[j]
      
      # Checking if we have data for both individuals for this patch
      if(!is.null(patch_columns[[patch]][[ind1]]) && !is.null(patch_columns[[patch]][[ind2]])) 
      {
        # Calculating spectral angle between the same patch type on different individuals
        angle <- calc_spectral_angle(
          patch_columns[[patch]][[ind1]],   # Patch on first individual
          patch_columns[[patch]][[ind2]]    # Same patch on second individual
        )
        
        # Storing the angle in our list
        within_variation[[patch]] <- c(within_variation[[patch]], angle)
      }
    }
  }
}

# =============================================================================
# Calculating between-patch variation
# =============================================================================
# This measures how different patch types are within the same individual
# (For example, how different is black from blue on individual 01)
between_variation <- list()
for(ind in individuals) 
{
  between_variation[[ind]] <- c()  # Initializing empty vector for this individual
  
  # Comparing all possible pairs of patch types for this individual
  for(i in 1:(length(patch_types)-1)) 
  {
    for(j in (i+1):length(patch_types)) 
    {
      patch1 <- patch_types[i]
      patch2 <- patch_types[j]
      
      # Checking if we have data for both patch types for this individual
      if(!is.null(patch_columns[[patch1]][[ind]]) && !is.null(patch_columns[[patch2]][[ind]])) 
      {
        # Calculating spectral angle between different patch types on the same individual
        angle <- calc_spectral_angle(
          patch_columns[[patch1]][[ind]],   # First patch type on this individual
          patch_columns[[patch2]][[ind]]    # Second patch type on this individual
        )
        
        # Storing the angle in our list
        between_variation[[ind]] <- c(between_variation[[ind]], angle)
      }
    }
  }
}

# =============================================================================
# PLOT 2: Boxplot of variation within vs. between patches
# =============================================================================
# Preparing data for boxplot by flattening the lists
within_all <- unlist(within_variation)    # Combining all within-patch angles
between_all <- unlist(between_variation)  # Combining all between-patch angles

# Creating a data frame for plotting
variation_data <- data.frame(
  Variation = c(rep("Within Patch", length(within_all)),     # Labels for within-patch data
                rep("Between Patches", length(between_all))),# Labels for between-patch data
  Angle = c(within_all, between_all)                         # All angle values
)

# Creating the boxplot
boxplot(Angle ~ Variation,               # Formula: Angle grouped by Variation type
        data=variation_data,             # Data to plot
        main="Variation Within vs. Between Patches",  
        ylab="Spectral Angle (radians)",  
        col=c("lightblue", "lightgreen"))  

# =============================================================================
# Calculating summary statistics for variation
# =============================================================================
# Printing summary stats for within-patch variation
cat("\nSummary of within-patch variation:\n")
within_summary <- summary(within_all)
print(within_summary)
cat("Standard deviation:", sd(within_all), "\n")

# Printing summary stats for between-patch variation
cat("\nSummary of between-patch variation:\n")
between_summary <- summary(between_all)
print(between_summary)
cat("Standard deviation:", sd(between_all), "\n")

# =============================================================================
# Statistical test: Are within-patch and between-patch variations different?
# =============================================================================
# Performing t-test to check if the difference is statistically significant
t_result <- t.test(within_all, between_all)
cat("\nT-test comparing within-patch vs. between-patch variation:\n")
print(t_result)

# =============================================================================
# PLOT 3: Spectral angles heatmap
# =============================================================================
# Creating a matrix of all pairwise spectral angle comparisons
# First, getting all column names that contain patch measurements
all_columns <- c()
for(patch in patch_types) 
{
  for(ind in individuals) 
  {
    col_name <- paste0(patch, "_", ind)  # Constructing column name (e.g., "black_01")
    if(col_name %in% names(df)) 
    {
      all_columns <- c(all_columns, col_name)  # Adding to our list if it exists
    }
  }
}

# Creating matrix to store all pairwise spectral angles
angle_matrix <- matrix(0,                         # Initializing with zeros
                       nrow=length(all_columns),  # Rows = number of columns
                       ncol=length(all_columns))  # Columns = number of columns

# Adding row and column names for clarity
rownames(angle_matrix) <- all_columns
colnames(angle_matrix) <- all_columns

# Calculating spectral angle between each pair of measurements
for(i in 1:length(all_columns)) 
{
  for(j in 1:length(all_columns)) 
  {
    if(i != j) # Skipping scenario comparing a column with itself (angle would be 0)
    {  
      # Calculating spectral angle between two measurements
      angle_matrix[i,j] <- calc_spectral_angle(
        df[[all_columns[i]]],   # First measurement
        df[[all_columns[j]]]    # Second measurement
      )
    }
  }
}

# Creating a heatmap visualization of all pairwise spectral angles
# This shows clusters of similar measurements
heatmap(angle_matrix, 
        main="Spectral Angles Between All Measurements",  
        xlab="Measurement",   
        ylab="Measurement",   
        col=heat.colors(100)) # Color palette (100 shades)

# =============================================================================
# PLOT 4: Individual variation by patch type
# =============================================================================
# Setting up a 2x2 grid for plotting the 4 patch types
par(mfrow=c(2,2))  # 2 rows, 2 columns of plots

# Creating one plot for each patch type
for(patch in patch_types) 
{
  # Creating an empty plot with appropriate limits
  plot(df$wavelength, df[[paste0(patch, "_01")]], 
       type="n",  # No points or lines initially (will add them in a sec)
       # Calculating appropriate y-axis limits
       ylim=c(0, max(sapply(individuals, function(ind) 
       {
         col <- paste0(patch, "_", ind)
         if(col %in% names(df)) max(df[[col]], na.rm=TRUE) else 0
       }) * 1.1)),  # Adding 10% margin
       xlab="Wavelength (nm)",   
       ylab="Reflectance (%)",  
       main=paste("Variation in", patch, "patch"))  
  
  # Adding lines for each individual's spectrum for this patch type
  for(ind in individuals) 
  {
    col_name <- paste0(patch, "_", ind)  # Constructing column name
    if(col_name %in% names(df)) 
    {
      # Plotting this individual's spectrum with a unique color
      lines(df$wavelength, df[[col_name]], 
            col=rainbow(length(individuals))[match(ind, individuals)])
    }
  }
  
  # Adding the mean spectrum as a thick black line
  lines(df$wavelength, patch_means[[patch]], col="black", lwd=2)
  
  # Legend to identify each line
  legend("topright",                                 
         c(individuals, "Mean"),                     
         col=c(rainbow(length(individuals)), "black"),  
         lwd=c(rep(1, length(individuals)), 2),     
         cex=0.7)                                   
}

# Resetting the plotting layout to default (1x1)
par(mfrow=c(1,1))

# =============================================================================
# Summary statistics by patch type and individual
# =============================================================================
# Calculating and displaying the average spectral angle within each patch type
cat("\nAverage spectral angle within each patch type:\n")
for(patch in patch_types) 
{
  if(length(within_variation[[patch]]) > 0) 
  {
    # If we have data, calculating and printing the mean
    cat(patch, ":", mean(within_variation[[patch]]), "\n")
  } 
  else 
  {
    # If no data is available, printing a message
    cat(patch, ": No data available\n")
  }
}

# Calculating and displaying the average spectral angle between patches for each individual
cat("\nAverage spectral angle between patches for each individual:\n")
for(ind in individuals) 
{
  if(length(between_variation[[ind]]) > 0) 
  {
    # If we have data, calculating and printing the mean
    cat("Individual", ind, ":", mean(between_variation[[ind]]), "\n")
  } 
  else 
  {
    # If no data is available, printing a message
    cat("Individual", ind, ": No data available\n")
  }
}
