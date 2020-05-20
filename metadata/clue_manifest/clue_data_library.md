Generate manifest file for uploading on the CLUE data library
================

The clue.io manifest file requires the following fields

    file_name   
    assay_protocol  
    data_level  
    level_code  
    level_desc  
    size
    md5

This notebook generates the manifest file.

All data files should be available locally in order to compute `size`
and `md5`.

``` r
library(tidyverse)
```

# Define batch and versioned github commit hash

``` r
rep_batch <- "2016_04_01_a549_48hr_batch1"
commit_hash <- "10d0fb44a8e87c076881a9c228fffb375c92d9d4"
```

# Define data levels

``` r
cell_painting_data_levels <- 
  tribble(~data_level, ~level_code, ~level_desc,
          "Level1", "CPLevel1", "Raw unprocessed images from microscopes",
          "Level2a", "CPLevel2a", "Per-cell level measurements stored across multiple CSVs", 
          "Level2b", "CPLevel2b", "Per-cell level measurements stored in a single backend file per assay plate", 
          "Level3", "CPLevel3", "Per-well level median-aggregated measurements",
          "Level4a", "CPLevel4a", "Morphological profiles computed using z-scores relative to the plate population",
          "Level4b", "CPLevel4b", "Feature selection applied to Level4b",
          "Level5a", "CPLevel5a", "Median-aggregated perturbation signatures",
          "Level5b", "CPLevel5b", "MODZ-aggregated perturbation signatures"
)
cell_painting_data_levels$assay_protocol <- "CellPaintingv2"
cell_painting_data_levels %>% show_table
```

| data\_level | level\_code | level\_desc                                                                     | assay\_protocol |
| :---------- | :---------- | :------------------------------------------------------------------------------ | :-------------- |
| Level1      | CPLevel1    | Raw unprocessed images from microscopes                                         | CellPaintingv2  |
| Level2a     | CPLevel2a   | Per-cell level measurements stored across multiple CSVs                         | CellPaintingv2  |
| Level2b     | CPLevel2b   | Per-cell level measurements stored in a single backend file per assay plate     | CellPaintingv2  |
| Level3      | CPLevel3    | Per-well level median-aggregated measurements                                   | CellPaintingv2  |
| Level4a     | CPLevel4a   | Morphological profiles computed using z-scores relative to the plate population | CellPaintingv2  |
| Level4b     | CPLevel4b   | Feature selection applied to Level4b                                            | CellPaintingv2  |
| Level5a     | CPLevel5a   | Median-aggregated perturbation signatures                                       | CellPaintingv2  |
| Level5b     | CPLevel5b   | MODZ-aggregated perturbation signatures                                         | CellPaintingv2  |

Specify suffixes for each data level

``` r
l2_suffix <- ".sqlite"
l3_suffix <- "_augmented.csv.gz"
l4a_suffix <- "_normalized.csv.gz"
l4b_suffix <- "_normalized_feature_select.csv.gz"
l5a_suffix <- "_consensus_median.csv.gz"
l5b_suffix <- "_consensus_modz.csv.gz"
```

# Create a list of plates in this experiment

``` r
all_rep_plates <- 
  read_csv(sprintf("../platemaps/%s/barcode_platemap.csv", rep_batch)) %>%
  distinct(Assay_Plate_Barcode)
```

    ## Parsed with column specification:
    ## cols(
    ##   Assay_Plate_Barcode = col_character(),
    ##   Plate_Map_Name = col_character(),
    ##   Batch_Number = col_double(),
    ##   Batch_Date = col_date(format = "")
    ## )

``` r
# remove some missing plates
missing_rep_plates <- 
  tribble(~Assay_Plate_Barcode,
          "SQ00015225",
          "SQ00015226",
          "SQ00015227",
          "SQ00015228")

rep_plates <-
  setdiff(all_rep_plates, missing_rep_plates)

rep_plates %>%
  head %>%
  show_table
```

| Assay\_Plate\_Barcode |
| :-------------------- |
| SQ00015201            |
| SQ00015202            |
| SQ00015200            |
| SQ00015204            |
| SQ00015205            |
| SQ00015203            |

``` r
rep_plates %>%
  count %>%
  show_table
```

|   n |
| --: |
| 136 |

# Get URLs to files

## Level 3 and Level 4

Generate local paths to Level 3 and Level 4 data

``` r
path_local <- paste0("../../profiles/", rep_batch)

level_3_4_files_local <- 
  rep_plates$Assay_Plate_Barcode %>%
  map_df(function(plate) {
    tibble(plate = plate,
           CPLevel3  = file.path(path_local, plate, paste0(plate, l3_suffix)),
           CPLevel4a = file.path(path_local, plate, paste0(plate, l4a_suffix)),
           CPLevel4b = file.path(path_local, plate, paste0(plate, l4b_suffix))
    )
  })
```

Generate remote paths to Level 3 and Level 4 data

``` r
path_url <- paste0(
  "https://media.githubusercontent.com/media/broadinstitute/lincs-cell-painting/",
  commit_hash,
  "/profiles/",
  rep_batch
)

level_3_4_files_url <- 
  rep_plates$Assay_Plate_Barcode %>%
  map_df(function(plate) {
    tibble(plate = plate,
           CPLevel3  = file.path(path_url, plate, paste0(plate, l3_suffix)),
           CPLevel4a = file.path(path_url, plate, paste0(plate, l4a_suffix)),
           CPLevel4b = file.path(path_url, plate, paste0(plate, l4b_suffix))
    )
  })
```

## Level 5

Generate local paths to Level 5 data

``` r
path_local <- "../../consensus"

level_5_files_local <-
  tibble(batch = rep_batch,
         CPLevel5a = file.path(path_local, rep_batch, paste0(rep_batch, l5a_suffix)),
         CPLevel5b = file.path(path_local, rep_batch, paste0(rep_batch, l5b_suffix)))
```

Generate remote paths to Level 5 data

``` r
path_url <- paste0(
  "https://media.githubusercontent.com/media/broadinstitute/lincs-cell-painting/",
  commit_hash,
  "/consensus"
)

level_5_files_url <-
  tibble(batch = rep_batch,
         CPLevel5a = file.path(path_url, rep_batch, paste0(rep_batch, l5a_suffix)),
         CPLevel5b = file.path(path_url, rep_batch, paste0(rep_batch, l5b_suffix)))
```

# Create manifest file

Note: this step requires the files to be locally available because it
checks the size and computes `md5`.

## Level 3 and 4

``` r
manifest_3_4_url <- 
  level_3_4_files_url %>% 
  pivot_longer(-plate, 
               names_to = "level_code",
               values_to = "file_name") %>% 
  inner_join(cell_painting_data_levels, by = "level_code") %>%
  select(plate, file_name, assay_protocol, data_level, level_code, level_desc) %>%
  arrange(plate, data_level) 

manifest_3_4_local <-
  level_3_4_files_local %>% 
  pivot_longer(-plate, 
               names_to = "level_code",
               values_to = "file_name") %>% 
  arrange(plate) %>%
  rowwise() %>%
  mutate(size = file.size(file_name)) %>%
  mutate(md5 = unname(tools::md5sum(file_name)))
  
manifest_3_4 <-
  inner_join(
    manifest_3_4_url,
    manifest_3_4_local %>% select(-file_name),
    by = c("plate", "level_code")) %>%
  select(-plate)
```

## Level 5

``` r
manifest_5_url <- 
  level_5_files_url %>% 
  pivot_longer(-batch, 
               names_to = "level_code",
               values_to = "file_name") %>% 
  inner_join(cell_painting_data_levels, by = "level_code") %>%
  select(batch, file_name, assay_protocol, data_level, level_code, level_desc) %>%
  arrange(batch, data_level) 

manifest_5_local <-
  level_5_files_local %>% 
  pivot_longer(-batch, 
               names_to = "level_code",
               values_to = "file_name") %>% 
  arrange(batch) %>%
  rowwise() %>%
  mutate(size = file.size(file_name)) %>%
  mutate(md5 = unname(tools::md5sum(file_name)))
  
manifest_5 <-
  inner_join(
    manifest_5_url,
    manifest_5_local %>% select(-file_name),
    by = c("batch", "level_code")) %>%
  select(-batch)
```

## Combine 3,4,5

``` r
manifest <- 
  bind_rows(
    manifest_3_4,
    manifest_5
  )
```

# Check if URLs exist

``` r
manifest_check <-
  manifest %>%
  rowwise() %>%
  mutate(
    url_exists = RCurl::url.exists(file_name),
    na_md5 = is.na(md5),
    na_size = is.na(size)
    ) %>%
  ungroup()
```

``` r
manifest_check %>%
  filter(!url_exists | na_md5 | na_size)
```

    ## # A tibble: 0 x 10
    ## # … with 10 variables: file_name <chr>, assay_protocol <chr>,
    ## #   data_level <chr>, level_code <chr>, level_desc <chr>, size <dbl>,
    ## #   md5 <chr>, url_exists <lgl>, na_md5 <lgl>, na_size <lgl>

# Write manifest

``` r
manifest %>%
  write_tsv("cell_painting_lincs_pilot_1_manifest.txt")
```
