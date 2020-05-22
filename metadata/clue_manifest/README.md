# Manifest file for clue.io

Scripts to generate manifest file for clue.io.

`clue_data_library.Rmd` requires that the profiles are available in this `master` branch of this repository.

To reproduce the manifest file, execute the following command:

```bash
# Navigate to the clue_manifest directory
R -e "rmarkdown::render('clue_data_library.Rmd', output_file = 'clue_data_library.md')"
``` 
