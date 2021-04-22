# EpxGDM

This repository contains all scripts used to perform the analyses described in the manuscript entitled: "**_Epigenome Wide Association Study Reveals Methylation Loci Associated with Offspring Gestational Diabetes Mellitus Exposure and Maternal Methylome_**" from Canouil et al., 2021.

Scripts were executed using the following Docker images stored on [Docker Hub](https://hub.docker.com/repository/docker/umr1283/stat/):

- [umr1283/stat:R403](https://hub.docker.com/layers/umr1283/stat/R403/images/sha256-2b0fb490ba31e186f8b8f53c8e11e79177b6b7eeb8e8267b8c90f5f4c4d6c4f1?context=repo)
- [umr1283/stat:4.0.4](https://hub.docker.com/layers/umr1283/stat/4.0.4/images/sha256-774404184836d1dafeada1b4635adb7fa20d8520d27c38751867d52c0c2c09bd?context=repo)

Data sharing statement: see manuscript.

## Design

``` bash
nohup Rscript scripts/01-sample_sheet.R > logs/01.log &
```

## Quality-Control EPIC Array

``` bash
nohup Rscript -e 'rmarkdown::render(input = here::here("scripts", "02-qc_idats.Rmd"), output_file = "EpxGDM_EPIC_QC.html", output_dir = here::here("reports"), encoding = "UTF-8", intermediates_dir = here::here("outputs", "02-qc_idats"), params = list(output_directory = here::here("outputs", "02-qc_idats")))' > logs/02.log &

nohup Rscript -e 'rmarkdown::render(input = here::here("scripts", "02-qc_idats_snps.Rmd"), output_file = "EpxGDM_EPIC_QC_SNPs.html", output_dir = here::here("reports"), encoding = "UTF-8", intermediates_dir = here::here("outputs", "02-qc_idats_snps"), params = list(output_directory = here::here("outputs", "02-qc_idats_snps")))' > logs/02_snps.log &
```

## Statistical Analyses (FinnGeDi)

### Perform EWAS for GDM

``` bash
nohup Rscript scripts/08-ewas_finngedi.R > logs/08.log &

nohup Rscript scripts/08-ewas_finngedi_snps.R > logs/08_snps.log &
```

### Perform 'heritability' EWAS for GDM

``` bash
nohup Rscript scripts/09-ewas_finngedi_h.R > logs/09.log &

nohup Rscript scripts/09-ewas_finngedi_h_snps.R > logs/09_snps.log &
```

### Perform 'LMM' EWAS for GDM

``` bash
nohup Rscript scripts/10-ewas_finngedi_lmm.R > logs/10.log &

nohup Rscript scripts/10-ewas_finngedi_lmm_snps.R > logs/10_snps.log &
```

## Article

### Materials

``` bash
nohup Rscript scripts/12-article.R > logs/12.log &
```

### Rebuttal

#### Main Comments Checks

``` bash
nohup Rscript scripts/20210122-rebuttal.R > logs/20210122-rebuttal.log &
```

#### Perform 'LMM' EWAS for GDM with primipara and/or mumwtgain

``` bash
nohup Rscript scripts/15-ewas_finngedi_lmm_rebuttal.R > logs/15.log &
```

#### Perform 'heritability' EWAS for GDM

``` bash
nohup Rscript scripts/16-ewas_finngedi_h_rebuttal.R > logs/16.log &
```

#### Perform EWAS for GDM

``` bash
nohup Rscript scripts/17-ewas_finngedi_rebuttal.R > logs/17.log &
```

#### Article Materials

``` bash
nohup Rscript scripts/18-article_rebuttal.R > logs/18.log &
```

