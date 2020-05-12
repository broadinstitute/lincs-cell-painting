Generate manifest file for uploading on the CLUE data library
================

The manifest file requires the following fields

    file_name   
    assay_protocol  
    data_level  
    level_code  
    level_desc  
    size
    md5

``` r
library(tidyverse)
```

``` r
cell_painting_data_levels <- 
  tribble(~data_level, ~level_code, ~level_desc,
          "Level1", "CPLevel1", "Raw unprocessed images from microscopes",
          "Level2a", "CPLevel2a", "Per-cell level measurements stored across multiple CSVs", 
          "Level2b", "CPLevel2b", "Per-cell level measurements stored in a single backend file per assay plate", 
          "Level3", "CPLevel3", "Per-well level aggregated measurements",
          "Level4a", "CPLevel4a", "Morphological profiles computed using z-scores relative to the plate population",
          "Level4b", "CPLevel4b", "Feature selection applied to Level4b",
          "Level5", "CPLevel5", "Replicate-collapsed perturbation signatures"
)
cell_painting_data_levels$assay_protocol <- "CellPaintingv2"
cell_painting_data_levels %>% show_table
```

| data\_level | level\_code | level\_desc                                                                     | assay\_protocol |
| :---------- | :---------- | :------------------------------------------------------------------------------ | :-------------- |
| Level1      | CPLevel1    | Raw unprocessed images from microscopes                                         | CellPaintingv2  |
| Level2a     | CPLevel2a   | Per-cell level measurements stored across multiple CSVs                         | CellPaintingv2  |
| Level2b     | CPLevel2b   | Per-cell level measurements stored in a single backend file per assay plate     | CellPaintingv2  |
| Level3      | CPLevel3    | Per-well level aggregated measurements                                          | CellPaintingv2  |
| Level4a     | CPLevel4a   | Morphological profiles computed using z-scores relative to the plate population | CellPaintingv2  |
| Level4b     | CPLevel4b   | Feature selection applied to Level4b                                            | CellPaintingv2  |
| Level5      | CPLevel5    | Replicate-collapsed perturbation signatures                                     | CellPaintingv2  |

Create list of plates

``` r
all_rep_plates <- 
  read_csv("../platemaps/2016_04_01_a549_48hr_batch1/barcode_platemap.csv") %>%
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

``` r
# URL prefix for L3 and L4 data
path_url <- "https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/master/profiles/2016_04_01_a549_48hr_batch1"

path_local <- "../../profiles/2016_04_01_a549_48hr_batch1"

l3_suffix <- "_augmented.csv.gz"
l4a_suffix <- "_normalized.csv.gz"
l4b_suffix <- "_normalized_variable_selected.csv.gz"

level_files <- 
  rep_plates$Assay_Plate_Barcode %>%
  map_df(function(plate) {
    tibble(plate = plate,
           CPLevel3        = file.path(path_url,   plate, paste0(plate, l3_suffix)),
           CPLevel4a       = file.path(path_url,   plate, paste0(plate, l4a_suffix)),
           CPLevel4b       = file.path(path_url,   plate, paste0(plate, l4b_suffix))
    )
  })

level_files_local <- 
  rep_plates$Assay_Plate_Barcode %>%
  map_df(function(plate) {
    tibble(plate = plate,
           CPLevel3  = file.path(path_local, plate, paste0(plate, l3_suffix)),
           CPLevel4a = file.path(path_local, plate, paste0(plate, l4a_suffix)),
           CPLevel4b = file.path(path_local, plate, paste0(plate, l4b_suffix))
    )
  })
```

``` r
manifest <- 
  level_files %>% 
  pivot_longer(-plate, 
               names_to = "level_code",
               values_to = "file_name") %>% 
  inner_join(cell_painting_data_levels) %>%
  select(plate, file_name, assay_protocol, data_level, level_code, level_desc) %>%
  arrange(plate, data_level) 
```

    ## Joining, by = "level_code"

``` r
manifest_local <-
  level_files_local %>% 
  pivot_longer(-plate, 
               names_to = "level_code",
               values_to = "file_name") %>% 
  arrange(plate) %>%
  rowwise() %>%
  mutate(size = file.size(file_name)) %>%
  mutate(md5 = unname(tools::md5sum("environment.yml")))
  
manifest %<>%
  inner_join(
    manifest_local %>% 
      select(-file_name),
    by = c("plate", "level_code")) %>%
  select(-plate)

 manifest %>%
  write_tsv("cell_painting_lincs_pilot_1_manifest.txt")
```
