

# This script takes fsa files from sequencing
# It aligns the ladder to the sample trace
# converts from scan size to bp and converts to correct CAG size
# and exports the sizing data

# other sections in this script
# find bp with highest fluoresence
# scaling traces
# smoothing traces


# install packages if not already installed
# install.packages("BiocManager") 
# BiocManager::install("sangerseqR")
# install.packages("pracma") 
# install.packages("ggplot2")  
# install.packages("dplyr")    
# install.packages("readr")    
# install.packages("readxl")


# load libraries
library(sangerseqR)
library(pracma)
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)

# set working directory. use " " Use / Not \
setwd("")

fsa <- read.abif("file.fsa")

# Access a dye channel (blue trace is in data 1)
trace <- fsa@data$DATA.1

# view the raw fluorescence signal
plot(trace, type = "l", col = "blue", main = "Electropherogram Trace", xlab = "Scan", ylab = "Fluorescence")

# ladder trace
ladder_trace <- fsa@data$DATA.4

# view mapmarker 1000 ladder
plot(fsa@data$DATA.4, type = "l", col = "red", main = "Ladder Trace (DATA.4)")


###################################################################################################

#find peaks
peaks <- findpeaks(ladder_trace, threshold = 200, minpeakdistance = 20) # can change threshold/minpeakdistance

#mapmarker1000 - 23 peaks
known_sizes_bp <- c(50, 75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 )  # example ladder sizes

scan_positions <- peaks[,2]  # peak locations

#################################################################################
# if A MISMATCH WAS FOUND BETWEEN NUMBER OF
# PEAKS AND NUMBER OF KNOWN LADDER PEAKS
# CAN CHECK:
length(known_sizes_bp)
length(scan_positions)

# CHECK ITS CALLED RIGHT PEAKS
plot(ladder_trace, type = "l", main = "Ladder Trace")
points(scan_positions, ladder_trace[scan_positions], col = "red", pch = 19)


#if mismatch, need to correct it, otherwise go to next section

# IF THERE ARE EXTRA PEAKS - IDENTIFY THEIR SIZES:
sort(scan_positions)

# Manually choose the correct peak indices
# exclude the wrong ones
# correct_scan_positions <- scan_positions[2:17]
# or
# exclude the ones in bracket
correct_scan_positions <- scan_positions[!(scan_positions %in% c(1248, 1360, 2878))]

#Now both vectors have same length: check
length(correct_scan_positions)  # Should be 16(liz500) 
length(known_sizes_bp)          # Should be 16(liz500) 

#FIXED LADDER
plot(ladder_trace, type = "l", main = "Ladder Trace")
points(correct_scan_positions, ladder_trace[correct_scan_positions], col = "red", pch = 19)

#CORRECT THE SIZES if need to
#known_sizes_bp <- c(50, 75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 )  # example ladder sizes

#end of correcting ladder section
#############################################################################


# Plot ladder trace with peaks marked
plot(ladder_trace, type = "l", col = "orange", main = "FAM & Ladder Traces with Peaks",
     xlab = "Scan Position", ylab = "Fluorescence Intensity", ylim = range(c(trace, ladder_trace)))
points(correct_scan_positions, ladder_trace[correct_scan_positions], col = "red", pch = 19)
# Add fam trace on same plot but scaled and in orange
lines(trace, col = "blue")
legend("topright", legend = c("Ladder", "Ladder Peaks", "FAM Trace"),
       col = c("orange", "red", "blue"), lty = c(1, NA, 1), pch = c(NA, 19, NA))

###########################################################################

# order the peaks
correct_scan_positions_sorted <- sort(correct_scan_positions)

df_peaks <- data.frame(scan = correct_scan_positions_sorted, size = known_sizes_bp)


#######################################################################
#precise calibration for FRAGMENT ANALYSIS
scan_to_bp <- approxfun(correct_scan_positions_sorted, known_sizes_bp, method = "linear", rule = 2)

# Apply to entire traces
all_ladder_scans <- seq_along(ladder_trace)
all_ladder_bp <- scan_to_bp(all_ladder_scans)

all_fam_scans <- seq_along(trace)
all_fam_bp <- scan_to_bp(all_fam_scans)


#####################################
plot(all_ladder_bp, ladder_trace, type = "l", col = "orange",
     main = "Ladder and FAM Traces with Base Pair X-axis",
     xlab = "Base Pair (bp)", ylab = "Fluorescence",
     xlim = c(0, 500),
     ylim = range(c(trace, ladder_trace)))

points(scan_to_bp(df_peaks$scan), ladder_trace[df_peaks$scan], col = "red", pch = 19)

lines(all_fam_bp, trace, col = "blue")

legend("topright", legend = c("Ladder", "Ladder Peaks", "FAM Trace"),
       col = c("orange", "red", "blue"), lty = c(1, NA, 1), pch = c(NA, 19, NA))


######################################

#for cag - minus 110 and divide by 3 to get correct bp

# all scan positions in the trace
all_fam_scans <- seq_along(trace)

# Use interpolation to get raw bp estimates
all_fam_bp_raw <- scan_to_bp(all_fam_scans)

# Apply the transformation: (bp - 110) / 3
all_fam_bp_adjusted <- (all_fam_bp_raw - 110) / 3

############################
#PLOT
plot(all_fam_bp_adjusted, trace, type = "l", col = "blue",
     xlab = "Base Pair (adjusted)", ylab = "Fluorescence Intensity",
     main = "FAM Trace with Adjusted Base Pairs",
     xlim = c(100, 150),
     ylim = c(0, 9000))



# Create data frame for FAM
fam_data <- data.frame(
  Adjusted_BP = all_fam_bp_adjusted,
  Fluorescence = trace
)
# Create data frame for ladder
ladder_data <- data.frame(
  Scan_Position = all_ladder_scans,
  BP = all_ladder_bp,
  Fluorescence = ladder_trace
)

# Filter FAM trace for 50-100 bp - you can change range
fam_df_subset <- subset(fam_data, Adjusted_BP >= 50 & Adjusted_BP <= 100) 

# export csv
write.csv(fam_df_subset, "fam_data.csv", row.names = FALSE)
write.csv(ladder_data, "ladder_data.csv", row.names = FALSE)

######CAN OVERLAY TRACES IN GRAPHPAD#####


#find bp with highest fluoresence
# Define the bp range you want to analyze
# bp_min <- 400
#bp_max <- 600

# Find indices where bp is within the range 400-600
#indices_in_range <- which(all_fam_bp_adjusted >= bp_min & all_fam_bp_adjusted <= bp_max)

# Extract fluorescence values in that bp range
#fluorescence_in_range <- trace[indices_in_range]

# Find index of the maximum fluorescence in the range
#max_fluorescence_index <- which.max(fluorescence_in_range)

# Find corresponding bp for that max fluorescence
#bp_at_max <- all_fam_bp_adjusted[indices_in_range][max_fluorescence_index]

# Print result
#cat("Base pair with highest fluorescence between", bp_min, "and", bp_max, "bp is:", bp_at_max, "\n")

######################################################

#scaling traces
#data1 to be scaled
#data1 <- read_excel("filelocation/filename.xlsx")
data1 <- read.csv("filelocation/filename.csv")

#data2 to be used to scale the other one
# data2 <- read_excel("filelocation/filename.xlsx")
data2 <- read.csv("filelocation/filename.csv")


#get max fluorescence between 100-200 bp
range1 <- subset(data1, Adjusted_BP >= 100 & Adjusted_BP <= 200)
range2 <- subset(data2, Adjusted_BP >= 100 & Adjusted_BP <= 200)

max1 <- max(range1$Fluorescence, na.rm = TRUE)
max2 <- max(range2$Fluorescence, na.rm = TRUE)

scaling_factor <- max2 / max1
data1$RFU_scaled <- data1$Fluorescence * scaling_factor


plot(data2$Adjusted_BP, data2$Fluorescence, type = "l", col = "blue", lwd = 2,
     main = "Normalized Traces", xlab = "Base Pairs", ylab = "RFU", xlim = c(100, 200), ylim = c(0, 2500))
lines(data1$Adjusted_BP, data1$RFU_scaled, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Trace 2", "Trace 1 (scaled)"),
       col = c("blue", "red"), lty = c(1, 2))


# Export the normalized trace to a new CSV
write.csv(data1, "normalized_H01.csv", row.names = FALSE)


######################################################
# smoothing traces
# Read your csv file (change path as needed)
# need to be / not \

csv <- read.csv("normalized_H01.csv")

ggplot(csv, aes(x = Adjusted_BP)) +
  geom_line(aes(y = RFU_scaled), color = "blue", alpha = 0.5) +
  coord_cartesian(xlim = c(120, 160), ylim = c(0, 2500)) +
  labs(title = "Trace with Moving Average Smoothing",
       x = "Base Pair",
       y = "Fluorescence (RFU)") +
  theme_classic()

# Assuming you have a data frame called `csv` with columns: Adjusted_BP and RFU_scaled
# Filter data to range of interest
subset_data <- csv %>%
  filter(Adjusted_BP >= 119, Adjusted_BP <= 161,
         RFU_scaled >= 0, RFU_scaled <= 2500)


# Detect peaks with constraints
peaks <- findpeaks(
  subset_data$RFU_scaled,
  minpeakheight = 80,       # Ignore small peaks (tweak as needed)
  minpeakdistance = 25        # Skip peaks too close together (in index steps)
)

# Extract peak data
if (!is.null(peaks)) {
  peak_indices <- peaks[, 2]
  peak_data <- subset_data[peak_indices, ]
} else {
  peak_data <- data.frame()
}

# Plot with original trace + filtered peaks + line through peaks
ggplot() +
  geom_line(data = subset_data, aes(x = Adjusted_BP, y = RFU_scaled), color = "green") +
  geom_point(data = peak_data, aes(x = Adjusted_BP, y = RFU_scaled), color = "red", size = 2) +
  geom_line(data = peak_data, aes(x = Adjusted_BP, y = RFU_scaled), color = "blue", size = 1) +
  # geom_area(data = peak_data, aes(x = Adjusted_BP, y = RFU_scaled), fill = "blue", alpha = 0.4) +   # shaded area
  coord_cartesian(xlim = c(120, 160), ylim = c(0, 2500)) +
  labs(title = "Line Through Prominent Peak Tops", x = "Base Pair", y = "Fluorescence")

#have sort data in order
sorted_data <- peak_data %>%
  arrange(Adjusted_BP)


# change to sorted_data
ggplot() +
  geom_line(data = subset_data, aes(x = Adjusted_BP, y = RFU_scaled), color = "green") +
  geom_point(data = peak_data, aes(x = Adjusted_BP, y = RFU_scaled), color = "red", size = 2) +
  geom_line(data = sorted_data, aes(x = Adjusted_BP, y = RFU_scaled), color = "blue", size = 1) +
  # geom_area(data = peak_data, aes(x = Adjusted_BP, y = RFU_scaled), fill = "blue", alpha = 0.4) +   # shaded area
  coord_cartesian(xlim = c(120, 160), ylim = c(0, 2500)) +
  labs(title = "Line Through Prominent Peak Tops", x = "Base Pair", y = "Fluorescence")


# Then export
write.csv(sorted_data, "H01-smoothed.csv", row.names = FALSE)

