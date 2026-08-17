# ==============================================================================
# 15-min flow event selection
# ==============================================================================

# --------Main logic-------------------
# ---1. Data Processing---
# Read flow data, convert timestamps to POSIXct, and filter to the study period (1994–2021).
# Align with a complete time series to expose gaps and fill missing values via linear interpolation.
# Compute 10-point moving average and estimate baseflow.

# ---2. Threshold-Based Peak/Valley Selection---
# Detect local maxima and minima in raw and smoothed series.
# Apply Filter 1: mark maxima above the quantile-based threshold as peak candidates.

# ---3. Temporal Separation of Events---
# Apply Filter 2: separate peaks using a minimum temporal gap.
# In the main analysis, 48 x 15-min time steps = 12 h.
# Peaks separated by <= 12 h are assigned to the same event, whereas peaks
# separated by > 12 h are treated as separate events.

# ---4. Start/End Time Refinement---
# Compute delta_y = 4 × IQR for event significance.
# For each event:
# Select start point: pre-peak minimum with amplitude > delta_y and closest to baseflow.
# Select end point: post-peak minimum closest to baseflow before the next event peak.
# Update filter_3 and build events_data table (IDs, indices, times, months, seasons).

# ---5. Event Metrics---
# Calculate the Flashiness Index (FI) and Rate of Change (RoC) for each event.
# Assign meteorological season using the month of the event peak.

# ---6. Event Plotting and Saving---
# Save one event table per catchment.
# Optionally plot raw series, moving average, extrema, baseflow, and highlighted
# event intervals for diagnostic checking.

# --------Main datasets----------------
# input:
# catchment_selected_cv.csv:
#   List of selected catchments/stations used in the analysis.
#
# <station_id>.csv:
#   15-minute flow data for an individual catchment. Each file is expected to
#   contain:
#       datetime : timestamp in "%Y-%m-%d %H:%M:%S"
#       value    : discharge
#
# data_site:
#   Main dataset for a single catchment during the study period.
#
# full_time:
#   Complete continuous 15-minute time series, used to align data_site and
#   reveal missing timestamps.
#
# output:
# events_data:
#   Summary table of detected events, including event ID, event start/end
#   indices, peak index, FI, RoC, dates, month, and season.

# --------Main variables----------------
# data_site$value : Raw flow value.
# data_site$ma    : 10-point centred moving average.
# data_site$BF    : Estimated baseflow.
# data_site$R     : Difference between observed flow and estimated baseflow.
# threshold       : Quantile-based peak threshold.
# Pot_index       : Indices of potential peak points from Filter 1.
# min_diff        : Minimum event-separation threshold in 15-min time steps.
# delta_y         : Amplitude threshold for event significance (4 × IQR).
# filter_1        : Threshold exceedance indicator.
# filter_2        : Preliminary event grouping based on peak separation.
# filter_3        : Final event interval after start/end refinement.
# ==============================================================================

options(bitmapType = "cairo")

# Load packages
library(dplyr)
library(lubridate)
library(ggplot2)
library(zoo)
library(hydroEvents)
library(data.table)

# ==============================================================================
# User settings
# ==============================================================================

# ==============================================================================
# User-defined paths
# Replace the paths below with your own working paths
# ==============================================================================

CATCHMENT_FILE <- "your/path/to/catchment_selected.csv"
FLOW_DIR <- "your/path/to/flow_15min"
OUTPUT_ROOT <- "your/path/to/results/event_selection"
PLOT_ROOT <- "your/path/to/results/event_plots"

# Study period
STUDY_START <- as.POSIXct("1994-01-01 00:00:00", tz = "UTC")
STUDY_END   <- as.POSIXct("2021-12-31 23:59:59", tz = "UTC")

# Define the peak threshold quantile
q <- 0.99

# Set TRUE only if diagnostic event plots are required.
SAVE_EVENT_PLOTS <- FALSE


# ---------------- Data processing ----------------

# Read catchment information.
# The catchments used in the study were pre-screened for data completeness.
catchment_info <- read.csv(file = CATCHMENT_FILE, header = TRUE,stringsAsFactors = FALSE)

# Retrieve the station IDs
catchment_id <- catchment_info$id

# Create a complete 15-minute time series to expose missing timestamps in the
# raw data.
full_time <- data.frame(
  date = seq(
    from = STUDY_START,
    to   = STUDY_END,
    by   = "15 min"
  )
)


for (h in seq_along(q)) {
  
  # Create separate folders to store results for different threshold settings.
  folder_path <- file.path(OUTPUT_ROOT,paste0("Event_selection_", q[h]))
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  # Detect and extract events for each catchment.
  for (j in seq_along(catchment_id)) {
    
    output_file <- file.path(
      folder_path,
      paste0(catchment_id[j], ".csv")
    )
    
    print(
      paste(
        "The process: quantile ",
        round(h / length(q), 2) * 100,
        "%; catchment ",
        round(j / length(catchment_id) * 100, 2),
        "%"
      )
    )
    
    # Read raw 15-minute flow data.
    input_file <- file.path(
      FLOW_DIR,
      paste0(catchment_id[j], ".csv")
    )
    
    if (!file.exists(input_file)) {
      warning("Input file not found: ", input_file)
      next
    }
    
    data_site <- read.csv(
      file = input_file,
      header = TRUE,
      stringsAsFactors = FALSE
    )
    
    # Replace negative values with 0.
    data_site$value[data_site$value < 0] <- 0
    
    
    # --------Data Processing------------
    
    # Subset to the target time window and perform linear interpolation.
    row.names(data_site) <- NULL
    
    # Represent calendar timestamps (ISO 8601) as POSIXct in UTC.
    data_site$date <- as.POSIXct(
      data_site$datetime,
      format = "%Y-%m-%d %H:%M:%S",
      tz = "UTC"
    )
    
    # Define the study period.
    data_site <- data_site %>%
      filter(
        date >= STUDY_START,
        date <= STUDY_END
      )
    
    # Ensure a complete, regular 15-minute time grid by aligning to `full_time`.
    if (nrow(data_site) != nrow(full_time)) {
      data_site <- full_time %>%
        left_join(data_site, by = "date")
    }
    
    # Tag records with the current catchment ID.
    data_site$id <- catchment_id[j]
    
    # Interpolate missing values linearly.
    # rule = 2 retains the behaviour used in the original analysis at the
    # boundaries of the series.
    data_site$value <- na.approx(
      data_site$value,
      na.rm = FALSE,
      rule = 2
    )
    
    # Create a sequential index for downstream operations and plots.
    data_site$index <- seq_len(nrow(data_site))
    
    # Compute a centred 10-point moving average; edges are filled with NA.
    data_site$ma <- rollmean(
      data_site$value,
      k = 10,
      align = "center",
      fill = NA
    )
    
    
    # --------Threshold-based Peak and Valley Selection------------
    
    # Function to identify local maxima and minima in a numeric vector.
    find_local_extrema <- function(x) {
      
      n <- length(x)
      
      if (n < 3) {
        return(rep("0", n))
      }
      
      # Neighbor values (previous and next).
      prev <- shift(x, 1, type = "lag")
      x_next <- shift(x, 1, type = "lead")
      
      # Strict local minima/maxima: strictly less/greater than both neighbors.
      is_min <- (x < prev & x < x_next)
      is_max <- (x > prev & x > x_next)
      
      # Detect flat extrema.
      # Flat extrema are ties with exactly one neighbor, while remaining lower
      # or higher than the other neighbor.
      #
      # Example for a minimum:
      # x equals prev but is less than next, or equals next but is less than prev.
      is_min_flat <- (
        (x == prev & x < x_next) |
          (x == x_next & x < prev)
      )
      
      is_max_flat <- (
        (x == prev & x > x_next) |
          (x == x_next & x > prev)
      )
      
      # Initialize all points as "0", then assign extrema labels.
      result <- rep("0", n)
      result[which(is_min | is_min_flat)] <- "Min"
      result[which(is_max | is_max_flat)] <- "Max"
      
      return(result)
    }
    
    # Identify local extrema on the raw series.
    data_site$ext_value <- find_local_extrema(data_site$value)
    
    # Identify local extrema on the smoothed (moving-average) series.
    data_site$ext_ma <- find_local_extrema(data_site$ma)
    
    
    # Simulate / estimate baseflow.
    BF <- baseflowB(data_site$value,alpha = 0.995, passes = 3, r = 30)
    
    # Store the estimated baseflow in the dataset.
    data_site$BF <- BF$bf
    
    # Calculate the difference between observed flow and estimated baseflow.
    # Smaller R values indicate points closer to the estimated baseflow.
    data_site$R <- data_site$value - data_site$BF
    
    
    # ----- For raw (un-smoothed) data -----
    
    # Indices of local minima in the raw series.
    lmin_val <- which(data_site$ext_value == "Min")
    
    # Indices of local maxima in the raw series.
    lmax_val <- which(data_site$ext_value == "Max")
    
    # ----- For smoothed data -----
    
    # Indices of local minima in the smoothed series.
    lmin_ma <- which(data_site$ext_ma == "Min")
    
    # Indices of local maxima in the smoothed series.
    lmax_ma <- which(data_site$ext_ma == "Max")
    
    # Use the threshold to select candidate peaks in the raw series.
    threshold_quantile <- q[h]
    threshold <- quantile(
      data_site$value,
      threshold_quantile,
      na.rm = TRUE
    )
    
    # Identify potential peaks (candidate maxima) above the threshold.
    Pot_index <- which(
      data_site$value >= threshold &
        data_site$ext_value == "Max"
    )
    
    # If no eligible candidate peaks are present, save an empty output and move
    # to the next catchment.
    if (length(Pot_index) == 0 || length(lmax_ma) == 0) {
      warning("No eligible peak candidates for catchment ", catchment_id[j])
      next
    }
    
    # Exclude partial events outside the main analysis range.
    Pot_index <- Pot_index[Pot_index >= lmax_ma[1]]
    Pot_index <- Pot_index[Pot_index <= lmax_ma[length(lmax_ma)]]
    
    if (length(Pot_index) == 0) {
      warning("No peak candidates remained after edge filtering for catchment ",catchment_id[j])
      next
    }
    
    # Create a new column to store the result after the first filter.
    data_site$filter_1 <- 0
    
    # Mark points above the threshold as threshold exceedances.
    data_site$filter_1 <- ifelse(data_site$value > threshold,1,0)
    
    # --------Separation Based on Temporal Distance------------
    
    # Minimum separation between two events in time steps.
    # 48 × 15 min = 12 h.
    min_diff <- 48
    
    # Calculate the time lag between consecutive candidate peaks.
    lag <- diff(Pot_index)
    
    # Identify where the gap between peaks exceeds the minimum separation.
    # Therefore, peaks separated by <= 12 h remain within the same event.
    lag_index <- which(lag > min_diff)
    
    # Event end points: peaks followed by a gap > 12 h.
    end_index <- c(Pot_index[lag_index],Pot_index[length(Pot_index)])
    
    # Event start points: first peak, and peaks following a gap > 12 h.
    start_index <- c(Pot_index[1],Pot_index[lag_index + 1])
    
    # Filter 2: assign unique event IDs to each preliminary event segment.
    data_site$filter_2 <- 0
    event_rank <- seq_along(end_index)
    
    for (i in seq_along(end_index)) {
      data_site$filter_2[start_index[i]:end_index[i]] <- event_rank[i]
    }
    
    
    # --------Selection of Start/End Time------------
    
    # delta_y is set to 4 × IQR of the flow values.
    # IQR = Q3 - Q1 (75th percentile - 25th percentile).
    delta_y <- 4 * (quantile(data_site$value, 0.75, na.rm = TRUE) - quantile(data_site$value, 0.25, na.rm = TRUE))
    
    # Create a new column to store results after applying the start/end criteria.
    data_site$filter_3 <- 0
    
    # Create a data frame to store start, end, and peak indices for each event.
    events_data <- data.frame(
      event = numeric(0),
      str_index = numeric(0),
      end_index = numeric(0),
      max_index = numeric(0)
    )
    
    # Ensure date column is in POSIXct format (UTC timezone).
    data_site$date <- as.POSIXct(
      data_site$date,
      format = "%Y-%m-%d %H:%M:%S",
      tz = "UTC"
    )
    
    
    # Define start and end points for each event.
    for (i in seq_along(event_rank)) {
      
      # Current event subset as defined by filter_2.
      sub_data <- subset(data_site,data_site$filter_2 == event_rank[i])
      
      # Find the highest flow within the current preliminary event.
      max_val <- max(sub_data$value)
      
      # Index of the event peak; take the first if there are ties.
      sub_max_idx <- sub_data$index[sub_data$value == max_val]
      sub_max_idx <- sub_max_idx[1]
      
      
      # --------Start-point selection--------
      
      # Candidate start points before the peak based on minima in the smoothed
      # series within the preliminary event segment.
      sub_dat_min_back <- subset(sub_data,ext_ma == "Min" & index < sub_max_idx)
      
      # Index of the first candidate peak in the current preliminary event.
      p_str_idx <- min(sub_data$index)
      
      # Include smoothed minima occurring before that first candidate peak.
      p_str_idx <- lmin_ma[lmin_ma < p_str_idx]
      
      # Append smoothed minima found within the event segment before the main peak.
      p_str_idx <- c(p_str_idx,sub_dat_min_back$index)
      
      # Remove duplicates.
      p_str_idx <- unique(p_str_idx)
      
      # Ensure candidate starts are below the event peak value.
      p_str_idx <- p_str_idx[data_site$value[p_str_idx] < max_val]
      
      # If no candidate start points are found, record NA for this event and
      # continue to the next event.
      if (length(p_str_idx) == 0) {
        
        cat("→ [EMPTY] i =", i, "\n")
        
        events_data[i, ] <- data.frame(
          event = event_rank[i],
          str_index = NA,
          end_index = NA,
          max_index = sub_max_idx
        )
        
        next
      }
      
      # From the second event onward, the start point must occur after the
      # previous event's peak.
      if (i > 1) {
        
        pre_peak <- events_data$max_index[i - 1]
        
        p_str_idx <- p_str_idx[p_str_idx > pre_peak]
        
        # If no candidate start points remain, record NA and continue.
        if (length(p_str_idx) == 0) {
          
          events_data[i, ] <- data.frame(
            event = event_rank[i],
            str_index = NA,
            end_index = NA,
            max_index = sub_max_idx
          )
          
          next
        }
      }
      
      
      # Amplitude relative to the smoothed series:
      # event peak minus candidate local minimum.
      diff_back <- max_val - data_site$ma[p_str_idx]
      
      # Proximity to baseflow:
      # smaller R means closer to the estimated baseflow line.
      diff_bf <- data_site$R[p_str_idx]
      
      # Candidates whose amplitude exceeds delta_y.
      vally_index <- which(diff_back > delta_y)
      
      # If no candidate meets the delta_y criterion, take the candidate with the
      # largest amplitude drop.
      if (length(vally_index) == 0) {
        
        vally_index_min <- which.max(diff_back)
        str_index <- p_str_idx[vally_index_min]
        
      } else {
        
        # Among amplitude-qualified candidates, choose the one closest to
        # estimated baseflow.
        candidate_idx <- p_str_idx[vally_index]
        
        candidate_diff_bf <- diff_bf[vally_index]
        
        # Pick the last among ties for stability.
        best_bf_idx <- tail(which(candidate_diff_bf == min(candidate_diff_bf)),1)
        
        str_index_ma <- candidate_idx[best_bf_idx]
        
        
        # Refine start point on the raw series near the moving-average start
        # candidate (±5 samples).
        l_b <- str_index_ma - 5
        
        u_b <- ifelse(str_index_ma + 5 < sub_max_idx,str_index_ma + 5,sub_max_idx - 1)
        
        # Protect the lower bound at the start of the record.
        l_b <- max(1, l_b)
        
        raw_sub_data <- data_site[l_b:u_b,]
        
        # Prefer raw local minima if present; otherwise retain the whole window.
        if (nrow(subset(raw_sub_data,raw_sub_data$ext_value == "Min")) == 0) {
          raw_sub_data <- raw_sub_data
        } else {
          
          raw_sub_data <- subset(
            raw_sub_data,
            raw_sub_data$ext_value == "Min"
          )
        }
        
        # Choose the point closest to baseflow within the raw-data window.
        raw_sub_idx <- which.min(raw_sub_data$R)[1]
        str_index <- raw_sub_data$index[raw_sub_idx]
      }
      
      
      # --------End-point selection--------
      
      # Smoothed minima after the event peak.
      min_after_peak <- lmin_ma[lmin_ma > sub_max_idx]
      
      if (length(min_after_peak) != 0) {
        
        # Define the search window from the current peak to the next event's
        # peak, or to the end of the series for the final event.
        if (i <= length(event_rank) - 1) {
          next_data <- subset(data_site, data_site$filter_2 == event_rank[i + 1])
          max_val_next <- max(next_data$value)
          next_max_idx <- next_data$index[next_data$value == max_val_next][1]
          
        } else {
          
          next_max_idx <- nrow(data_site)
        }
        
        
        potential_end_data <- data_site[sub_max_idx:next_max_idx,]
        
        # Select the point closest to baseflow within this window.
        potential_end_idx <- which.min(potential_end_data$R)
        end_index_ma <- potential_end_data$index[potential_end_idx]
        
        
        # Refine end point on the raw series near the selected end point (±5
        # samples).
        l_b <- ifelse(end_index_ma - 5 > sub_max_idx, end_index_ma - 5, sub_max_idx + 1)
        
        if (
          end_index_ma <= nrow(data_site) - 5
        ) {
          
          u_b <- end_index_ma + 5
          
        } else {
          
          u_b <- nrow(data_site)
        }
        
        
        raw_sub_data <- data_site[l_b:u_b,]
        
        # Prefer raw local minima if present; otherwise retain the whole window.
        if (nrow(subset(raw_sub_data, raw_sub_data$ext_value == "Min")) == 0) {
          raw_sub_data <- raw_sub_data
        } else {
          
          raw_sub_data <- subset(raw_sub_data, raw_sub_data$ext_value == "Min")
        }
        
        # Choose the point closest to baseflow as the event end.
        raw_sub_idx <- which.min(raw_sub_data$R)[1]
        
        end_index <- raw_sub_data$index[raw_sub_idx]
        
        
        # Capture dates for diagnostics / QA.
        str_date <- data_site$date[str_index]
        
        max_date <- data_site$date[sub_max_idx]
        
        end_date <- data_site$date[end_index]
        
        
        # Persist event boundaries.
        events_data[i, ] <- data.frame(
          event = event_rank[i],
          str_index = str_index,
          end_index = end_index,
          max_index = sub_max_idx
        )
        
        # Stamp the final event segment into filter_3.
        data_site$filter_3[str_index:end_index] <- event_rank[i]
        
      } else {
        
        # If no post-peak minima exist, record NA end and continue.
        events_data[i, ] <- data.frame(
          event = event_rank[i],
          str_index = NA,
          end_index = NA,
          max_index = sub_max_idx
        )
        
        next
      }
    }
    
    
    # Remove rows with missing start/end values from the event summary table.
    events_data <- na.omit(events_data)
    
    # Reassign event IDs sequentially after cleaning.
    events_data$event <- seq_len(nrow(events_data))
    
    
    # --------Calculate FI------------
    
    events_data$FI <- 0
    
    for (i in seq_len(nrow(events_data))) {
      
      t_P <- events_data$max_index[i]
      t_i <- events_data$str_index[i]
      t_e <- events_data$end_index[i]
      
      delta_t <- t_e - t_i
      
      Q_P <- data_site$value[t_P]
      
      # A small constant is added to Q_i to avoid division by zero.
      Q_i <- data_site$value[t_i] + 0.00001
      
      events_data$FI[i] <- (((Q_P - Q_i) / Q_i) /((t_P - t_i) / delta_t))
    }
    
    
    # --------Calculate RoC------------
    
    events_data$RoC <- 0
    
    for (i in seq_len(nrow(events_data))) {
      
      t_P <- events_data$max_index[i]
      t_i <- events_data$str_index[i]
      
      Q_P <- data_site$value[t_P]
      Q_i <- data_site$value[t_i]
      
      # This preserves the calculation used in the original analysis.
      # The denominator is the number of 15-minute time steps between event
      # start and peak.
      events_data$RoC[i] <- ((Q_P - Q_i) /(t_P - t_i))
    }
    
    
    # --------Seasonal classification------------
    
    events_data$season <- 0
    
    # Extract the timestamp of the event peak from the main dataset.
    max_time <- events_data$max_index
    
    events_data$time <- data_site$date[max_time]
    
    # Extract peak month as a numeric value (1–12).
    events_data$month <- as.numeric(substr(events_data$time,6,7))
    
    # Extract start, end, and peak dates for each event.
    str_time <- events_data$str_index
    end_time <- events_data$end_index
    peak_time <- events_data$max_index
    
    events_data$str_date <- data_site$date[str_time]
    
    events_data$end_date <- data_site$date[end_time]
    
    events_data$peak_date <- data_site$date[peak_time]
    
    
    # Assign meteorological seasons based on the peak month.
    # DJF = Dec–Feb
    # MAM = Mar–May
    # JJA = Jun–Aug
    # SON = Sep–Nov
    for (i in seq_len(nrow(events_data))) {
      
      month <- events_data$month[i]
      
      if (month %in% c(12, 1, 2)) {
        
        events_data$season[i] <- "DJF"
        
      } else if (month %in% c(3, 4, 5)) {
        
        events_data$season[i] <- "MAM"
        
      } else if (month %in% c(6, 7, 8)) {
        
        events_data$season[i] <- "JJA"
        
      } else if (month %in% c(9, 10, 11)) {
        
        events_data$season[i] <- "SON"
      }
    }
    
    
    # Display a frequency table of events by season.
    print(table(events_data$season))
    
    events_data$time <- format(as.POSIXct(events_data$time),"%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    
    # --------Save event table------------
    
    write.csv(events_data,file = output_file,row.names = FALSE)
    
    
    # --------Optional event plots------------
    
    if (SAVE_EVENT_PLOTS) {
      
      # Create output directory for plots.
      save_dir <- file.path(PLOT_ROOT,as.character(catchment_id[j]))
      
      if (!dir.exists(save_dir)) {dir.create(save_dir,recursive = TRUE)}

      # Loop over all detected events.
      for (events in seq_len(nrow(events_data))) {
        
        # Event metadata.
        event_num <- events
        
        str_date <- as.POSIXct(events_data$str_date[event_num],tz = "UTC")
        
        end_date <- as.POSIXct(events_data$end_date[event_num],tz = "UTC")
        
        event_idx <- events_data$str_index[event_num]:events_data$end_index[event_num]
        
        
        # Ensure UTC timezone.
        str_date <- lubridate::with_tz(str_date,tzone = "UTC")
        
        end_date <- lubridate::with_tz(end_date,tzone = "UTC")
        
        
        # Subset data for ±1 day around the event period.
        data_p <- data_site[data_site$date >= str_date - days(1) & data_site$date <= end_date + days(1),]
        
        data_p$date <- as.POSIXct(data_p$date, format = "%Y-%m-%d %H:%M:%S",tz = "UTC")
        
        
        # Highlight event points in red; all others in black.
        data_p$color <- ifelse(data_p$index %in% event_idx,"red", "black")
        
        
        p <- ggplot(data_p, aes( x = date, y = value)) +
          
          # 10-point moving-average series.
          geom_line(aes(y = ma), color = "blue", linewidth = 1) +
          
          # Local minima in the smoothed series.
          geom_point(data = subset(data_p, ext_ma == "Min"),
            aes(y = ma),
            color = "purple",
            shape = 17,
            size = 2
          ) +
          
          # Raw flow series; event interval is highlighted.
          geom_line(aes(group = 1,color = color), linewidth = 1) +
          
          geom_point(aes(color = color),size = 1) +
          
          # Estimated baseflow.
          geom_line(aes(y = BF),color = "#00BFFF",linewidth = 1) +
          
          scale_color_identity() +
          
          scale_x_datetime(date_labels = "%b %d %H:%M") +
          
          labs(title = paste0("Event ", event_num),x = "Date", y = "Flow") +
          
          theme_minimal() +
          
          theme(
            plot.title = element_text(size = 20,hjust = 0.5),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12)
          )
        
        
        # Save plot to file.
        file_name <- file.path(save_dir,paste0("Event_",event_num,".png"))
        ggsave(filename = file_name, plot = p, width = 10,height = 6,dpi = 300)
      }
    }
  }
}

