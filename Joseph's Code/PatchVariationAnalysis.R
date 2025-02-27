# =============================================================================
# Patches Variation Analysis 
# This code analyzes Trigger Fish patch spectral data to determine:
# 1. Whether we can tell different patch types apart (black, white, orange, blue)
# 2. How much variation exists within patches vs. across individuals
# =============================================================================

# Load the dataset containing spectral measurements
# Each column represents a specific patch on a specific individual (e.g., black_01)
df <- read.csv("pablo_equation_1_transposed.csv")

# =============================================================================
# Define the spectral angle calculation function
# =============================================================================
# This function calculates the angular difference between two spectral signatures
# Based on the formula from the paper: α = cos⁻¹(ΣXY/√(Σ(X)²Σ(Y)²))
# Smaller angles indicate greater similarity between spectra
calc_spectral_angle <- function(x, y) {
  sum_xy <- sum(x * y)     # Calculate dot product of the two spectra
  sum_x2 <- sum(x^2)       # Sum of squares for first spectrum
  sum_y2 <- sum(y^2)       # Sum of squares for second spectrum
  
  # Calculate the spectral angle using the inverse cosine function
  # This is insensitive to brightness differences, focusing on spectral shape
  alfa <- acos(sum_xy / sqrt(sum_x2 * sum_y2))
  
  return(alfa)  # Return angle in radians
}

# =============================================================================
# Identify patch types and individuals for analysis
# =============================================================================
# Define our patch types (colors) and individual identifiers
patch_types <- c("black", "white", "orange", "blue")  # The four patch types
individuals <- c("01", "02", "03", "04")              # The four individuals

# =============================================================================
# Organize the data by patch type and individual
# =============================================================================
# Create a nested list structure to easily access data:
# patch_columns$black$01 will contain the spectral data for black patch on individual 01
patch_columns <- list()
for(patch in patch_types) {
  # Create a sub-list for each patch type
  patch_columns[[patch]] <- list()
  
  # For each individual, find and store the corresponding column data
  for(ind in individuals) {
    col_name <- paste0(patch, "_", ind)  # Construct column name (e.g., "black_01")
    
    # Only add data if the column exists in our dataset
    if(col_name %in% names(df)) {
      patch_columns[[patch]][[ind]] <- df[[col_name]]  # Store the spectral data
    }
  }
}

# =============================================================================
# Calculate mean spectra for each patch type
# =============================================================================
# For each patch type, calculate the average spectrum across all individuals
# This gives us a "typical" spectrum for each patch type
patch_means <- data.frame(wavelength = df$wavelength)  # Start with wavelength column
for(patch in patch_types) {
  # Combine all measurements for this patch type into a matrix
  patch_data <- do.call(cbind, patch_columns[[patch]])
  
  # Calculate row means (mean reflectance at each wavelength)
  patch_means[[patch]] <- rowMeans(patch_data, na.rm = TRUE)
}

# =============================================================================
# PLOT 1: Mean spectra for each patch type
# =============================================================================
# Define colors for plotting that match the patch types
colors <- c("black", "gray", "orange", "blue")

# Create a plot of the mean spectrum for each patch type
# Start with the first patch (black)
plot(df$wavelength, patch_means$black, 
     type="l",                                       # Line plot
     col=colors[1],                                  # Use black color
     lwd=2,                                          # Line width 2
     ylim=c(0, max(patch_means[,-1], na.rm=TRUE) * 1.1),  # Y-axis range with 10% margin
     xlab="Wavelength (nm)",                         # X-axis label
     ylab="Reflectance (%)",                         # Y-axis label
     main="Mean Spectra by Patch Type")              # Plot title

# Add lines for the remaining patch types
for(i in 2:length(patch_types)) {
  lines(df$wavelength,                     # X values (wavelength)
        patch_means[[patch_types[i]]],     # Y values (reflectance)
        col=colors[i],                     # Color corresponding to patch type
        lwd=2)                             # Line width 2
}

# Add a legend to identify each line
legend("topright",                         # Position in top-right corner
       patch_types,                        # Labels for legend
       col=colors,                         # Colors for legend
       lwd=2)                              # Line width for legend

# =============================================================================
# Calculate spectral angles between different patch types
# =============================================================================
# This quantifies how different each patch type is from the others
cat("Spectral angles between patch types (in radians):\n")

# Create a matrix to store the angles between each pair of patch types
patch_angles <- matrix(0,                         # Initialize with zeros
                       nrow=length(patch_types),   # Rows = number of patch types
                       ncol=length(patch_types))   # Columns = number of patch types

# Add row and column names to the matrix for clarity
rownames(patch_angles) <- patch_types
colnames(patch_angles) <- patch_types

# Calculate the spectral angle between each pair of patch types
for(i in 1:length(patch_types)) {
  for(j in 1:length(patch_types)) {
    if(i != j) {  # Skip comparing a patch with itself (angle would be 0)
      # Calculate spectral angle between mean spectra of two patch types
      patch_angles[i,j] <- calc_spectral_angle(
        patch_means[[patch_types[i]]],     # Mean spectrum of first patch type
        patch_means[[patch_types[j]]]      # Mean spectrum of second patch type
      )
    }
  }
}

# Print the matrix of angles, rounded to 4 decimal places
print(round(patch_angles, 4))

# =============================================================================
# Calculate within-patch variation
# =============================================================================
# This measures how much variation exists within the same patch type
# across different individuals - captures biological variation
within_variation <- list()
for(patch in patch_types) {
  within_variation[[patch]] <- c()  # Initialize empty vector for this patch
  
  # Compare all possible pairs of individuals for this patch type
  for(i in 1:(length(individuals)-1)) {
    for(j in (i+1):length(individuals)) {
      ind1 <- individuals[i]
      ind2 <- individuals[j]
      
      # Check if we have data for both individuals for this patch
      if(!is.null(patch_columns[[patch]][[ind1]]) && !is.null(patch_columns[[patch]][[ind2]])) {
        # Calculate spectral angle between the same patch type on different individuals
        angle <- calc_spectral_angle(
          patch_columns[[patch]][[ind1]],   # Patch on first individual
          patch_columns[[patch]][[ind2]]    # Same patch on second individual
        )
        
        # Store the angle in our list
        within_variation[[patch]] <- c(within_variation[[patch]], angle)
      }
    }
  }
}

# =============================================================================
# Calculate between-patch variation
# =============================================================================
# This measures how different patch types are within the same individual
# (For example, how different is black from blue on individual 01)
between_variation <- list()
for(ind in individuals) {
  between_variation[[ind]] <- c()  # Initialize empty vector for this individual
  
  # Compare all possible pairs of patch types for this individual
  for(i in 1:(length(patch_types)-1)) {
    for(j in (i+1):length(patch_types)) {
      patch1 <- patch_types[i]
      patch2 <- patch_types[j]
      
      # Check if we have data for both patch types for this individual
      if(!is.null(patch_columns[[patch1]][[ind]]) && !is.null(patch_columns[[patch2]][[ind]])) {
        # Calculate spectral angle between different patch types on the same individual
        angle <- calc_spectral_angle(
          patch_columns[[patch1]][[ind]],   # First patch type on this individual
          patch_columns[[patch2]][[ind]]    # Second patch type on this individual
        )
        
        # Store the angle in our list
        between_variation[[ind]] <- c(between_variation[[ind]], angle)
      }
    }
  }
}

# =============================================================================
# PLOT 2: Boxplot of variation within vs. between patches
# =============================================================================
# Prepare data for boxplot by flattening the lists
within_all <- unlist(within_variation)    # Combine all within-patch angles
between_all <- unlist(between_variation)  # Combine all between-patch angles

# Create a data frame for plotting
variation_data <- data.frame(
  Variation = c(rep("Within Patch", length(within_all)),     # Labels for within-patch data
                rep("Between Patches", length(between_all))), # Labels for between-patch data
  Angle = c(within_all, between_all)                         # All angle values
)

# Create the boxplot comparing within-patch and between-patch variation
boxplot(Angle ~ Variation,               # Formula: Angle grouped by Variation type
        data=variation_data,             # Data to plot
        main="Variation Within vs. Between Patches",  # Plot title
        ylab="Spectral Angle (radians)",  # Y-axis label
        col=c("lightblue", "lightgreen"))  # Colors for the boxes

# =============================================================================
# Calculate summary statistics for variation
# =============================================================================
# Print summary stats for within-patch variation
cat("\nSummary of within-patch variation:\n")
within_summary <- summary(within_all)
print(within_summary)
cat("Standard deviation:", sd(within_all), "\n")

# Print summary stats for between-patch variation
cat("\nSummary of between-patch variation:\n")
between_summary <- summary(between_all)
print(between_summary)
cat("Standard deviation:", sd(between_all), "\n")

# =============================================================================
# Statistical test: Are within-patch and between-patch variations different?
# =============================================================================
# Perform t-test to check if the difference is statistically significant
t_result <- t.test(within_all, between_all)
cat("\nT-test comparing within-patch vs. between-patch variation:\n")
print(t_result)

# =============================================================================
# PLOT 3: Spectral angles heatmap
# =============================================================================
# Create a comprehensive matrix of all pairwise spectral angle comparisons
# First, get all column names that contain patch measurements
all_columns <- c()
for(patch in patch_types) {
  for(ind in individuals) {
    col_name <- paste0(patch, "_", ind)  # Construct column name (e.g., "black_01")
    if(col_name %in% names(df)) {
      all_columns <- c(all_columns, col_name)  # Add to our list if it exists
    }
  }
}

# Create matrix to store all pairwise spectral angles
angle_matrix <- matrix(0,                     # Initialize with zeros
                       nrow=length(all_columns),  # Rows = number of columns
                       ncol=length(all_columns))  # Columns = number of columns

# Add row and column names for clarity
rownames(angle_matrix) <- all_columns
colnames(angle_matrix) <- all_columns

# Calculate spectral angle between each pair of measurements
for(i in 1:length(all_columns)) {
  for(j in 1:length(all_columns)) {
    if(i != j) {  # Skip comparing a column with itself (angle would be 0)
      # Calculate spectral angle between two measurements
      angle_matrix[i,j] <- calc_spectral_angle(
        df[[all_columns[i]]],   # First measurement
        df[[all_columns[j]]]    # Second measurement
      )
    }
  }
}

# Create a heatmap visualization of all pairwise spectral angles
# This shows clusters of similar measurements
heatmap(angle_matrix, 
        main="Spectral Angles Between All Measurements",  # Title
        xlab="Measurement",   # X-axis label
        ylab="Measurement",   # Y-axis label
        col=heat.colors(100)) # Color palette (100 shades)

# =============================================================================
# PLOT 4: Individual variation by patch type
# =============================================================================
# Set up a 2x2 grid for plotting the 4 patch types
par(mfrow=c(2,2))  # 2 rows, 2 columns of plots

# Create one plot for each patch type
for(patch in patch_types) {
  # Create an empty plot with appropriate limits
  plot(df$wavelength, df[[paste0(patch, "_01")]], 
       type="n",  # No points or lines initially (we'll add them next)
       # Calculate appropriate y-axis limits
       ylim=c(0, max(sapply(individuals, function(ind) {
         col <- paste0(patch, "_", ind)
         if(col %in% names(df)) max(df[[col]], na.rm=TRUE) else 0
       }) * 1.1)),  # Add 10% margin
       xlab="Wavelength (nm)",   # X-axis label
       ylab="Reflectance (%)",   # Y-axis label
       main=paste("Variation in", patch, "patch"))  # Title
  
  # Add lines for each individual's spectrum for this patch type
  for(ind in individuals) {
    col_name <- paste0(patch, "_", ind)  # Construct column name
    if(col_name %in% names(df)) {
      # Plot this individual's spectrum with a unique color
      lines(df$wavelength, df[[col_name]], 
            col=rainbow(length(individuals))[match(ind, individuals)])
    }
  }
  
  # Add the mean spectrum as a thick black line
  lines(df$wavelength, patch_means[[patch]], col="black", lwd=2)
  
  # Add a legend to identify each line
  legend("topright",                                 # Position
         c(individuals, "Mean"),                     # Labels
         col=c(rainbow(length(individuals)), "black"),  # Colors
         lwd=c(rep(1, length(individuals)), 2),     # Line widths
         cex=0.7)                                   # Text size
}

# Reset the plotting layout to default (1x1)
par(mfrow=c(1,1))

# =============================================================================
# Summary statistics by patch type and individual
# =============================================================================
# Calculate and display the average spectral angle within each patch type
cat("\nAverage spectral angle within each patch type:\n")
for(patch in patch_types) {
  if(length(within_variation[[patch]]) > 0) {
    # If we have data, calculate and print the mean
    cat(patch, ":", mean(within_variation[[patch]]), "\n")
  } else {
    # If no data is available, print a message
    cat(patch, ": No data available\n")
  }
}

# Calculate and display the average spectral angle between patches for each individual
cat("\nAverage spectral angle between patches for each individual:\n")
for(ind in individuals) {
  if(length(between_variation[[ind]]) > 0) {
    # If we have data, calculate and print the mean
    cat("Individual", ind, ":", mean(between_variation[[ind]]), "\n")
  } else {
    # If no data is available, print a message
    cat("Individual", ind, ": No data available\n")
  }
}