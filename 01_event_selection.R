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


# --------Available functions----------------
# 1. find_local_extrema(x)
#    Identify local maxima and minima in a numeric time series.

# 2. empty_event_table()
#    Create an empty event table with the same structure as the normal
#    event-selection output.

# 3. save_event_plots(data_site, events_data, station_id, plot_dir, time_zone)
#    Generate and save diagnostic plots for selected events in one catchment.

# 4. select_events_one_catchment(...)
#    Run the complete event-selection workflow for a single catchment.
#    Returns the selected event table and can optionally return the processed
#    flow time series.

# 5. run_event_selection(...)
#    Run the event-selection workflow for all catchments listed in the
#    catchment information file and save one event table per catchment.


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

library(dplyr)
library(lubridate)
library(ggplot2)
library(zoo)
library(hydroEvents)
library(data.table)


# ==============================================================================
# 1. Helper functions
# ==============================================================================

# Identify local maxima and minima in a numeric vector.
find_local_extrema <- function(x) {
  
  n <- length(x)
  if (n < 3) return(rep("0", n))
  
  prev <- data.table::shift(x, 1, type = "lag")
  x_next <- data.table::shift(x, 1, type = "lead")
  
  # Strict local extrema.
  is_min <- (x < prev & x < x_next)
  is_max <- (x > prev & x > x_next)
  
  # Flat extrema: ties with exactly one neighbour while remaining lower/higher
  # than the other neighbour.
  is_min_flat <- ((x == prev & x < x_next) | (x == x_next & x < prev))
  is_max_flat <- ((x == prev & x > x_next) | (x == x_next & x > prev))
  
  result <- rep("0", n)
  result[which(is_min | is_min_flat)] <- "Min"
  result[which(is_max | is_max_flat)] <- "Max"
  
  return(result)
}


# Create an empty event table with the same fields as the normal output.
empty_event_table <- function() {
  
  data.frame(
    event = integer(0),
    str_index = integer(0),
    end_index = integer(0),
    max_index = integer(0),
    FI = numeric(0),
    RoC = numeric(0),
    season = character(0),
    time = character(0),
    month = integer(0),
    str_date = as.POSIXct(character(0), tz = "UTC"),
    end_date = as.POSIXct(character(0), tz = "UTC"),
    peak_date = as.POSIXct(character(0), tz = "UTC")
  )
}


# Save diagnostic event plots for one catchment.
save_event_plots <- function(data_site, events_data, station_id, plot_dir, time_zone = "UTC") {
  
  if (nrow(events_data) == 0) return(invisible(NULL))
  
  save_dir <- file.path(plot_dir, as.character(station_id))
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  for (event_num in seq_len(nrow(events_data))) {
    
    str_date <- as.POSIXct(events_data$str_date[event_num], tz = time_zone)
    end_date <- as.POSIXct(events_data$end_date[event_num], tz = time_zone)
    event_idx <- events_data$str_index[event_num]:events_data$end_index[event_num]
    
    str_date <- lubridate::with_tz(str_date, tzone = time_zone)
    end_date <- lubridate::with_tz(end_date, tzone = time_zone)
    
    # Plot +/-1 day around the selected event.
    data_p <- data_site[
      data_site$date >= str_date - days(1) &
        data_site$date <= end_date + days(1),
    ]
    
    data_p$color <- ifelse(data_p$index %in% event_idx, "red", "black")
    
    p <- ggplot(data_p, aes(x = date, y = value)) +
      geom_line(aes(y = ma), color = "blue", linewidth = 1, na.rm = TRUE) +
      geom_point(
        data = subset(data_p, ext_ma == "Min"),
        aes(y = ma),
        color = "purple",
        shape = 17,
        size = 2,
        na.rm = TRUE
      ) +
      geom_line(aes(group = 1, color = color), linewidth = 1) +
      geom_point(aes(color = color), size = 1) +
      geom_line(aes(y = BF), color = "#00BFFF", linewidth = 1) +
      scale_color_identity() +
      scale_x_datetime(date_labels = "%b %d %H:%M") +
      labs(title = paste0("Event ", event_num), x = "Date", y = "Flow") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
    
    file_name <- file.path(save_dir, paste0("Event_", event_num, ".png"))
    ggsave(filename = file_name, plot = p, width = 10, height = 6, dpi = 300)
  }
  
  invisible(NULL)
}


# ==============================================================================
# 2. Select events for ONE catchment
# ==============================================================================

select_events_one_catchment <- function(
    flow_file,
    station_id = NULL,
    datetime_col = "datetime",
    flow_col = "value",
    datetime_format = "%Y-%m-%d %H:%M:%S",
    time_zone = "UTC",
    study_start = "1994-01-01 00:00:00",
    study_end = "2021-12-31 23:59:59",
    time_step_min = 15,
    ma_window = 10,
    peak_quantile = 0.99,
    min_event_gap_hours = 12,
    refine_window_steps = 5,
    delta_iqr_multiplier = 4,
    baseflow_alpha = 0.995,
    baseflow_passes = 3,
    baseflow_r = 30,
    fi_epsilon = 0.00001,
    save_plots = FALSE,
    plot_dir = NULL,
    return_data_site = FALSE
) {
  
  # ---------------------------------------------------------------------------
  # Check inputs
  # ---------------------------------------------------------------------------
  
  if (!file.exists(flow_file)) stop(paste0("Flow file not found: ", flow_file))
  
  if (is.null(station_id)) {
    station_id <- tools::file_path_sans_ext(basename(flow_file))
  }
  
  if (time_step_min <= 0) stop("time_step_min must be > 0.")
  if (ma_window < 1) stop("ma_window must be >= 1.")
  if (peak_quantile <= 0 || peak_quantile >= 1) stop("peak_quantile must be between 0 and 1.")
  if (min_event_gap_hours <= 0) stop("min_event_gap_hours must be > 0.")
  if (refine_window_steps < 1) stop("refine_window_steps must be >= 1.")
  
  study_start <- as.POSIXct(study_start, format = datetime_format, tz = time_zone)
  study_end <- as.POSIXct(study_end, format = datetime_format, tz = time_zone)
  
  if (is.na(study_start) || is.na(study_end)) {
    stop("study_start or study_end could not be parsed. Check datetime_format.")
  }
  
  if (study_end <= study_start) stop("study_end must be later than study_start.")
  
  
  # ---------------------------------------------------------------------------
  # Data processing
  # ---------------------------------------------------------------------------
  
  data_site_raw <- read.csv(flow_file, header = TRUE, stringsAsFactors = FALSE)
  
  required_cols <- c(datetime_col, flow_col)
  missing_cols <- required_cols[!required_cols %in% names(data_site_raw)]
  
  if (length(missing_cols) > 0) {
    stop(paste0(
      "Missing required column(s) in ", flow_file, ": ",
      paste(missing_cols, collapse = ", "),
      ". Set datetime_col and flow_col to match your dataset."
    ))
  }
  
  # Standardise the two required input columns internally.
  data_site <- data.frame(
    datetime = data_site_raw[[datetime_col]],
    value = data_site_raw[[flow_col]]
  )
  
  if (!is.numeric(data_site$value)) {
    data_site$value <- suppressWarnings(as.numeric(data_site$value))
  }
  
  if (all(is.na(data_site$value))) {
    stop(paste0("The discharge column '", flow_col, "' contains no usable numeric values."))
  }
  
  # Replace negative flow values with 0, matching the original analysis.
  data_site$value[data_site$value < 0] <- 0
  
  # Convert timestamps.
  data_site$date <- as.POSIXct(data_site$datetime, format = datetime_format, tz = time_zone)
  
  if (all(is.na(data_site$date))) {
    stop(paste0(
      "All timestamps became NA. Check datetime_col, datetime_format, and time_zone."
    ))
  }
  
  # Remove rows whose timestamps could not be parsed.
  data_site <- data_site[!is.na(data_site$date), ]
  
  # Restrict to study period.
  data_site <- data_site %>%
    filter(date >= study_start, date <= study_end) %>%
    arrange(date)
  
  if (nrow(data_site) == 0) {
    stop("No observations remain within the requested study period.")
  }
  
  if (anyDuplicated(data_site$date)) {
    stop("Duplicate timestamps were detected. Resolve duplicate timestamps before running event selection.")
  }
  
  # Create the complete regular time series.
  full_time <- data.frame(
    date = seq(from = study_start, to = study_end, by = paste(time_step_min, "min"))
  )
  
  # Align to the complete regular grid to expose missing timestamps.
  data_site <- full_time %>%
    left_join(data_site %>% select(date, value), by = "date")
  
  data_site$id <- station_id
  
  # Linear interpolation of missing flow values.
  # rule = 2 preserves the boundary behaviour used in the original analysis.
  data_site$value <- zoo::na.approx(data_site$value, na.rm = FALSE, rule = 2)
  
  if (all(is.na(data_site$value))) {
    stop("No usable discharge data remain after interpolation.")
  }
  
  data_site$index <- seq_len(nrow(data_site))
  
  # Centred moving average.
  data_site$ma <- zoo::rollmean(data_site$value, k = ma_window, align = "center", fill = NA)
  
  
  # ---------------------------------------------------------------------------
  # Threshold-based peak and valley selection
  # ---------------------------------------------------------------------------
  
  data_site$ext_value <- find_local_extrema(data_site$value)
  data_site$ext_ma <- find_local_extrema(data_site$ma)
  
  # Estimate baseflow.
  BF <- hydroEvents::baseflowB(
    data_site$value,
    alpha = baseflow_alpha,
    passes = baseflow_passes,
    r = baseflow_r
  )
  
  data_site$BF <- BF$bf
  
  # Difference between observed flow and estimated baseflow.
  # Smaller R values indicate points closer to estimated baseflow.
  data_site$R <- data_site$value - data_site$BF
  
  lmin_ma <- which(data_site$ext_ma == "Min")
  lmax_ma <- which(data_site$ext_ma == "Max")
  
  threshold <- quantile(data_site$value, peak_quantile, na.rm = TRUE)
  
  # Candidate local maxima above the threshold.
  Pot_index <- which(
    data_site$value >= threshold &
      data_site$ext_value == "Max"
  )
  
  # Diagnostic threshold-exceedance field retained from the original workflow.
  data_site$filter_1 <- ifelse(data_site$value > threshold, 1, 0)
  
  if (length(Pot_index) == 0 || length(lmax_ma) == 0) {
    warning(paste0("No eligible peak candidates for catchment ", station_id, "."))
    events_data <- empty_event_table()
    
    if (return_data_site) {
      return(list(events_data = events_data, data_site = data_site))
    } else {
      return(events_data)
    }
  }
  
  # Exclude partial events outside the range bounded by smoothed local maxima.
  Pot_index <- Pot_index[Pot_index >= lmax_ma[1]]
  Pot_index <- Pot_index[Pot_index <= lmax_ma[length(lmax_ma)]]
  
  if (length(Pot_index) == 0) {
    warning(paste0("No peak candidates remained after edge filtering for catchment ", station_id, "."))
    events_data <- empty_event_table()
    
    if (return_data_site) {
      return(list(events_data = events_data, data_site = data_site))
    } else {
      return(events_data)
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Temporal separation of events
  # ---------------------------------------------------------------------------
  
  # Convert the requested event-separation time to number of input time steps.
  min_diff_raw <- min_event_gap_hours * 60 / time_step_min
  min_diff <- as.integer(round(min_diff_raw))
  
  if (abs(min_diff - min_diff_raw) > 1e-8) {
    warning("min_event_gap_hours is not an exact multiple of time_step_min; the number of steps was rounded.")
  }
  
  lag_peak <- diff(Pot_index)
  lag_index <- which(lag_peak > min_diff)
  
  # Event end points: peaks followed by a gap greater than the minimum separation.
  end_index <- c(Pot_index[lag_index], Pot_index[length(Pot_index)])
  
  # Event start points: first peak, and peaks following a large gap.
  start_index <- c(Pot_index[1], Pot_index[lag_index + 1])
  
  event_rank <- seq_along(end_index)
  
  data_site$filter_2 <- 0
  
  for (i in seq_along(event_rank)) {
    data_site$filter_2[start_index[i]:end_index[i]] <- event_rank[i]
  }
  
  
  # ---------------------------------------------------------------------------
  # Selection of event start/end time
  # ---------------------------------------------------------------------------
  
  # Event-significance amplitude threshold.
  delta_y <- delta_iqr_multiplier * (
    quantile(data_site$value, 0.75, na.rm = TRUE) -
      quantile(data_site$value, 0.25, na.rm = TRUE)
  )
  
  data_site$filter_3 <- 0
  
  events_data <- data.frame(
    event = event_rank,
    str_index = NA_integer_,
    end_index = NA_integer_,
    max_index = NA_integer_
  )
  
  
  for (i in seq_along(event_rank)) {
    
    # Current preliminary event as defined by filter_2.
    sub_data <- subset(data_site, filter_2 == event_rank[i])
    
    if (nrow(sub_data) == 0) next
    
    # Highest flow within the preliminary event.
    max_val <- max(sub_data$value, na.rm = TRUE)
    
    # Take the first if the maximum is tied.
    sub_max_idx <- sub_data$index[sub_data$value == max_val][1]
    events_data$max_index[i] <- sub_max_idx
    
    
    # -------------------------------------------------------------------------
    # Start-point selection
    # -------------------------------------------------------------------------
    
    # Smoothed minima inside the preliminary event before the main peak.
    sub_dat_min_back <- subset(sub_data, ext_ma == "Min" & index < sub_max_idx)
    
    # First candidate peak in this preliminary event.
    first_peak_index <- min(sub_data$index)
    
    # Include smoothed minima before the first candidate peak.
    p_str_idx <- lmin_ma[lmin_ma < first_peak_index]
    
    # Add smoothed minima inside the event before the main peak.
    p_str_idx <- c(p_str_idx, sub_dat_min_back$index)
    p_str_idx <- unique(p_str_idx)
    
    # Candidate start points must be below the event peak.
    p_str_idx <- p_str_idx[data_site$value[p_str_idx] < max_val]
    
    # From the second event onward, the start point must occur after the
    # preceding event peak.
    if (i > 1) {
      pre_peak <- events_data$max_index[i - 1]
      
      if (!is.na(pre_peak)) {
        p_str_idx <- p_str_idx[p_str_idx > pre_peak]
      }
    }
    
    # No valid candidate start point.
    if (length(p_str_idx) == 0) next
    
    # Amplitude from each candidate minimum to the main event peak.
    diff_back <- max_val - data_site$ma[p_str_idx]
    
    # Proximity to baseflow.
    diff_bf <- data_site$R[p_str_idx]
    
    # Candidates whose amplitude exceeds delta_y.
    valley_index <- which(diff_back > delta_y)
    
    if (length(valley_index) == 0) {
      
      # If no candidate passes the amplitude criterion, use the candidate with
      # the largest amplitude difference, matching the original analysis.
      valley_index_min <- which.max(diff_back)
      str_index <- p_str_idx[valley_index_min]
      
    } else {
      
      # Among amplitude-qualified candidates, choose the one closest to baseflow.
      candidate_idx <- p_str_idx[valley_index]
      candidate_diff_bf <- diff_bf[valley_index]
      
      # Pick the last among ties.
      best_bf_idx <- tail(
        which(candidate_diff_bf == min(candidate_diff_bf, na.rm = TRUE)),
        1
      )
      
      str_index_ma <- candidate_idx[best_bf_idx]
      
      # Refine the start point using the raw series within +/-N time steps.
      l_b <- max(1, str_index_ma - refine_window_steps)
      u_b <- min(str_index_ma + refine_window_steps, sub_max_idx - 1)
      
      if (u_b < l_b) next
      
      raw_sub_data <- data_site[l_b:u_b, ]
      
      # Prefer raw local minima when present.
      raw_minima <- subset(raw_sub_data, ext_value == "Min")
      if (nrow(raw_minima) > 0) raw_sub_data <- raw_minima
      
      # Choose the raw point closest to baseflow.
      raw_sub_idx <- which.min(raw_sub_data$R)[1]
      str_index <- raw_sub_data$index[raw_sub_idx]
    }
    
    
    # -------------------------------------------------------------------------
    # End-point selection
    # -------------------------------------------------------------------------
    
    # Check that at least one smoothed minimum exists after the event peak.
    min_after_peak <- lmin_ma[lmin_ma > sub_max_idx]
    
    if (length(min_after_peak) == 0) next
    
    # Search from the current peak to the next event peak, or to the end of the
    # record for the final event.
    if (i <= length(event_rank) - 1) {
      
      next_data <- subset(data_site, filter_2 == event_rank[i + 1])
      next_max_val <- max(next_data$value, na.rm = TRUE)
      next_max_idx <- next_data$index[next_data$value == next_max_val][1]
      
    } else {
      
      next_max_idx <- nrow(data_site)
    }
    
    potential_end_data <- data_site[sub_max_idx:next_max_idx, ]
    
    # Point closest to baseflow within the search window.
    potential_end_idx <- which.min(potential_end_data$R)
    end_index_ma <- potential_end_data$index[potential_end_idx]
    
    # Refine the end point using the raw series within +/-N time steps.
    l_b <- max(sub_max_idx + 1, end_index_ma - refine_window_steps)
    u_b <- min(nrow(data_site), end_index_ma + refine_window_steps)
    
    if (u_b < l_b) next
    
    raw_sub_data <- data_site[l_b:u_b, ]
    
    # Prefer raw local minima when present.
    raw_minima <- subset(raw_sub_data, ext_value == "Min")
    if (nrow(raw_minima) > 0) raw_sub_data <- raw_minima
    
    # Choose the point closest to baseflow as event end.
    raw_sub_idx <- which.min(raw_sub_data$R)[1]
    end_index <- raw_sub_data$index[raw_sub_idx]
    
    # Basic boundary check.
    if (is.na(str_index) || is.na(end_index) || str_index >= sub_max_idx || end_index <= sub_max_idx) {
      next
    }
    
    # Store final event boundaries.
    events_data$str_index[i] <- str_index
    events_data$end_index[i] <- end_index
    events_data$max_index[i] <- sub_max_idx
    
    # Final event mask.
    data_site$filter_3[str_index:end_index] <- event_rank[i]
  }
  
  
  # ---------------------------------------------------------------------------
  # Remove incomplete events
  # ---------------------------------------------------------------------------
  
  events_data <- na.omit(events_data)
  
  if (nrow(events_data) == 0) {
    warning(paste0("No complete events were retained for catchment ", station_id, "."))
    events_data <- empty_event_table()
    
    if (return_data_site) {
      return(list(events_data = events_data, data_site = data_site))
    } else {
      return(events_data)
    }
  }
  
  row.names(events_data) <- NULL
  events_data$event <- seq_len(nrow(events_data))
  
  
  # ---------------------------------------------------------------------------
  # Calculate FI
  # ---------------------------------------------------------------------------
  
  events_data$FI <- 0
  
  for (i in seq_len(nrow(events_data))) {
    
    t_P <- events_data$max_index[i]
    t_i <- events_data$str_index[i]
    t_e <- events_data$end_index[i]
    
    delta_t <- t_e - t_i
    
    Q_P <- data_site$value[t_P]
    
    # Small constant added to avoid division by zero.
    Q_i <- data_site$value[t_i] + fi_epsilon
    
    events_data$FI[i] <- (((Q_P - Q_i) / Q_i) / ((t_P - t_i) / delta_t))
  }
  
  
  # ---------------------------------------------------------------------------
  # Calculate RoC
  # ---------------------------------------------------------------------------
  
  events_data$RoC <- 0
  
  for (i in seq_len(nrow(events_data))) {
    
    t_P <- events_data$max_index[i]
    t_i <- events_data$str_index[i]
    
    Q_P <- data_site$value[t_P]
    Q_i <- data_site$value[t_i]
    
    # This preserves the calculation used in the original analysis.
    # The denominator is the number of input time steps between event start
    # and event peak.
    events_data$RoC[i] <- ((Q_P - Q_i) / (t_P - t_i))
  }
  
  
  # ---------------------------------------------------------------------------
  # Dates and seasonal classification
  # ---------------------------------------------------------------------------
  
  events_data$str_date <- data_site$date[events_data$str_index]
  events_data$end_date <- data_site$date[events_data$end_index]
  events_data$peak_date <- data_site$date[events_data$max_index]
  
  events_data$month <- lubridate::month(events_data$peak_date)
  
  events_data$season <- dplyr::case_when(
    events_data$month %in% c(12, 1, 2) ~ "DJF",
    events_data$month %in% c(3, 4, 5) ~ "MAM",
    events_data$month %in% c(6, 7, 8) ~ "JJA",
    events_data$month %in% c(9, 10, 11) ~ "SON",
    TRUE ~ NA_character_
  )
  
  # Retain the historical `time` field used in the original output.
  events_data$time <- format(events_data$peak_date, datetime_format, tz = time_zone)
  
  
  # ---------------------------------------------------------------------------
  # Optional diagnostic plots
  # ---------------------------------------------------------------------------
  
  if (save_plots) {
    
    if (is.null(plot_dir)) {
      stop("plot_dir must be supplied when save_plots = TRUE.")
    }
    
    save_event_plots(
      data_site = data_site,
      events_data = events_data,
      station_id = station_id,
      plot_dir = plot_dir,
      time_zone = time_zone
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  
  if (return_data_site) {
    return(list(events_data = events_data, data_site = data_site))
  } else {
    return(events_data)
  }
}


# ==============================================================================
# 3. Run event selection for ALL catchments
# ==============================================================================

run_event_selection <- function(
    catchment_file,
    flow_dir,
    output_root,
    plot_root = NULL,
    id_col = "id",
    datetime_col = "datetime",
    flow_col = "value",
    datetime_format = "%Y-%m-%d %H:%M:%S",
    time_zone = "UTC",
    study_start = "1994-01-01 00:00:00",
    study_end = "2021-12-31 23:59:59",
    time_step_min = 15,
    ma_window = 10,
    peak_quantile = 0.99,
    min_event_gap_hours = 12,
    refine_window_steps = 5,
    delta_iqr_multiplier = 4,
    baseflow_alpha = 0.995,
    baseflow_passes = 3,
    baseflow_r = 30,
    fi_epsilon = 0.00001,
    save_event_plots = FALSE,
    flow_file_extension = ".csv",
    return_results = FALSE
) {
  
  # ---------------------------------------------------------------------------
  # Check paths and catchment table
  # ---------------------------------------------------------------------------
  
  if (!file.exists(catchment_file)) stop(paste0("Catchment file not found: ", catchment_file))
  if (!dir.exists(flow_dir)) stop(paste0("Flow directory not found: ", flow_dir))
  
  catchment_info <- read.csv(catchment_file, header = TRUE, stringsAsFactors = FALSE)
  
  if (!id_col %in% names(catchment_info)) {
    stop(paste0(
      "Column '", id_col,
      "' was not found in catchment_file. Set id_col to match your dataset."
    ))
  }
  
  catchment_id <- as.character(catchment_info[[id_col]])
  
  # Create output directory.
  if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
  
  if (save_event_plots) {
    
    if (is.null(plot_root)) {
      plot_root <- file.path(output_root, "event_plots")
    }
    
    if (!dir.exists(plot_root)) dir.create(plot_root, recursive = TRUE)
  }
  
  
  # ---------------------------------------------------------------------------
  # Allow one or several peak thresholds
  # ---------------------------------------------------------------------------
  
  if (any(peak_quantile <= 0 | peak_quantile >= 1)) {
    stop("All peak_quantile values must be between 0 and 1.")
  }
  
  if (return_results) all_results <- list()
  
  
  for (h in seq_along(peak_quantile)) {
    
    q_current <- peak_quantile[h]
    q_label <- format(q_current, trim = TRUE, scientific = FALSE)
    
    # One output folder per threshold.
    folder_path <- file.path(output_root, paste0("Event_selection_", q_label))
    
    if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
    
    if (save_event_plots) {
      plot_folder <- file.path(plot_root, paste0("Event_selection_", q_label))
      if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)
    } else {
      plot_folder <- NULL
    }
    
    
    # -------------------------------------------------------------------------
    # Process each catchment
    # -------------------------------------------------------------------------
    
    for (j in seq_along(catchment_id)) {
      
      station_id <- catchment_id[j]
      
      message(
        paste0(
          "Threshold ", q_current,
          " | Catchment ", j, "/", length(catchment_id),
          " | ID: ", station_id
        )
      )
      
      input_file <- file.path(flow_dir, paste0(station_id, flow_file_extension))
      output_file <- file.path(folder_path, paste0(station_id, ".csv"))
      
      if (!file.exists(input_file)) {
        warning(paste0("Input file not found: ", input_file))
        next
      }
      
      events_data <- tryCatch(
        {
          select_events_one_catchment(
            flow_file = input_file,
            station_id = station_id,
            datetime_col = datetime_col,
            flow_col = flow_col,
            datetime_format = datetime_format,
            time_zone = time_zone,
            study_start = study_start,
            study_end = study_end,
            time_step_min = time_step_min,
            ma_window = ma_window,
            peak_quantile = q_current,
            min_event_gap_hours = min_event_gap_hours,
            refine_window_steps = refine_window_steps,
            delta_iqr_multiplier = delta_iqr_multiplier,
            baseflow_alpha = baseflow_alpha,
            baseflow_passes = baseflow_passes,
            baseflow_r = baseflow_r,
            fi_epsilon = fi_epsilon,
            save_plots = save_event_plots,
            plot_dir = plot_folder,
            return_data_site = FALSE
          )
        },
        error = function(e) {
          warning(paste0("Catchment ", station_id, " failed: ", conditionMessage(e)))
          return(NULL)
        }
      )
      
      if (is.null(events_data)) next
      
      # Save one event table per catchment.
      write.csv(events_data, file = output_file, row.names = FALSE)
      
      if (return_results) {
        all_results[[paste0("q", q_label, "_", station_id)]] <- events_data
      }
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Save the settings used in this run
  # ---------------------------------------------------------------------------
  
  settings <- data.frame(
    parameter = c(
      "id_col",
      "datetime_col",
      "flow_col",
      "datetime_format",
      "time_zone",
      "study_start",
      "study_end",
      "time_step_min",
      "ma_window",
      "peak_quantile",
      "min_event_gap_hours",
      "refine_window_steps",
      "delta_iqr_multiplier",
      "baseflow_alpha",
      "baseflow_passes",
      "baseflow_r",
      "fi_epsilon"
    ),
    value = c(
      id_col,
      datetime_col,
      flow_col,
      datetime_format,
      time_zone,
      study_start,
      study_end,
      time_step_min,
      ma_window,
      paste(peak_quantile, collapse = ","),
      min_event_gap_hours,
      refine_window_steps,
      delta_iqr_multiplier,
      baseflow_alpha,
      baseflow_passes,
      baseflow_r,
      fi_epsilon
    )
  )
  
  write.csv(
    settings,
    file.path(output_root, "event_selection_settings.csv"),
    row.names = FALSE
  )
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(output_root, "sessionInfo.txt")
  )
  
  
  message("Event selection completed.")
  
  if (return_results) {
    return(all_results)
  } else {
    invisible(NULL)
  }
}


# ==============================================================================
# 4. Examples
# ==============================================================================
#
# The functions above contain the full event-selection workflow.
# Users normally only need to change the arguments in the calls below.
#
#
# ------------------------------------------------------------------------------
# Example A: process ONE catchment without saving anything
# ------------------------------------------------------------------------------
#
# events <- select_events_one_catchment(
#   flow_file = "your/path/to/flow_15min/39001.csv",
#   station_id = "39001",
#   datetime_col = "datetime",
#   flow_col = "value"
# )
#
# head(events)
#
#
# If the user's columns are named differently:
#
# events <- select_events_one_catchment(
#   flow_file = "your/path/to/flow/39001.csv",
#   station_id = "39001",
#   datetime_col = "timestamp",
#   flow_col = "discharge"
# )
#
#
# ------------------------------------------------------------------------------
# Example B: run ALL catchments using the defaults from this study
# ------------------------------------------------------------------------------
#
# run_event_selection(
#   catchment_file = "your/path/to/catchment_selected.csv",
#   flow_dir = "your/path/to/flow_15min",
#   output_root = "your/path/to/results/event_selection",
#   id_col = "id",
#   datetime_col = "datetime",
#   flow_col = "value"
# )
#
#
# ------------------------------------------------------------------------------
# Example C: dataset with different column names
# ------------------------------------------------------------------------------
#
# run_event_selection(
#   catchment_file = "your/path/to/stations.csv",
#   flow_dir = "your/path/to/flow",
#   output_root = "your/path/to/results/event_selection",
#   id_col = "station_number",
#   datetime_col = "timestamp",
#   flow_col = "discharge"
# )
#
#
# ------------------------------------------------------------------------------
# Example D: save diagnostic event plots
# ------------------------------------------------------------------------------
#
# run_event_selection(
#   catchment_file = "your/path/to/catchment_selected.csv",
#   flow_dir = "your/path/to/flow_15min",
#   output_root = "your/path/to/results/event_selection",
#   plot_root = "your/path/to/results/event_plots",
#   id_col = "id",
#   datetime_col = "datetime",
#   flow_col = "value",
#   save_event_plots = TRUE
# )
#
#
# ------------------------------------------------------------------------------
# Example E: input data with a different temporal resolution
# ------------------------------------------------------------------------------
#
# IMPORTANT:
# Changing `time_step_min` automatically converts the 12-h event-separation
# criterion to the corresponding number of time steps. However, `ma_window`
# and `refine_window_steps` are also defined in numbers of observations and
# should be reconsidered for a different temporal resolution.
#
# Example for hourly data:
#
# run_event_selection(
#   catchment_file = "your/path/to/stations.csv",
#   flow_dir = "your/path/to/hourly_flow",
#   output_root = "your/path/to/results/event_selection",
#   id_col = "station_id",
#   datetime_col = "datetime",
#   flow_col = "flow",
#   time_step_min = 60,
#   ma_window = 10,
#   refine_window_steps = 5
# )
#
# ==============================================================================
