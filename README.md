# UK-flood-flashiness

This repository contains R code developed for a national-scale, event-based assessment of flood flashiness across UK catchments using 15-minute river discharge data from 1994 to 2021.

The analysis examines flood hydrograph characteristics, including the Flashiness Index (FI) and Rate of Change (RoC).

## R Scripts

### `Event_selection_process.R`

This script identifies flood events from continuous river discharge records and calculates event-level flashiness metrics.

The main steps include:

1. Preprocessing and interpolation of missing discharge observations.
2. Identifying flood peaks using a quantile-based threshold.
3. Separating flood events using a minimum temporal gap.
4. Determining event start and end times using hydrograph minima and estimated baseflow.
5. Calculating FI and RoC and assigning events to meteorological seasons.

The main analysis uses a 99th-percentile peak threshold and a 12-hour minimum separation between peak groups.

## Requirements

The script requires R and the following packages:

`dplyr`, `lubridate`, `ggplot2`, `zoo`, `hydroEvents`, `data.table`

## Input Data

The script requires:

- A catchment list in CSV format containing station IDs (`id`).
- Individual 15-minute discharge CSV files, named `<station_id>.csv`, containing `datetime` and `value` columns.

## Usage

After adjusting the input and output paths, the event-selection workflow can be run as follows:

```r
source("Event_selection_process.R")

run_event_selection(
  catchment_file = "path/to/catchment_selected.csv",
  flow_dir = "path/to/flow_15min",
  output_root = "path/to/results"
)
```

The script saves event-level CSV files for each catchment, together with the processing settings and R session information. Optional diagnostic hydrograph plots can also be generated.

## Data Availability

The original discharge data are not included in this repository. Users should obtain the relevant datasets from their respective data providers.
