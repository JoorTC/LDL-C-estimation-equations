# LDL-C Recalculation and Therapeutic Reclassification

This repository contains R code to calculate low-density lipoprotein cholesterol (LDL-C) using three commonly used estimation equations:

1. Friedewald equation
2. Martin-Hopkins equation
3. Sampson equation

The code was developed to support analyses evaluating how modern LDL-C estimation methods may affect therapeutic classification in patients monitored according to contemporary cardiovascular prevention guidelines.

## Required input

The input dataset must be an R data frame or tibble containing the following columns:

| Column | Description | Unit |
|---|---|---|
| `TC` | Total cholesterol | mg/dL |
| `HDL_C` | HDL cholesterol | mg/dL |
| `TG` | Triglycerides | mg/dL |

The dataset may optionally include:

| Column | Description | Unit |
|---|---|---|
| `non_HDL_C` | Non-HDL cholesterol | mg/dL |

If `non_HDL_C` is not provided, it is automatically calculated as:

`non_HDL_C = TC - HDL_C`

## Output

The function returns the original dataset with the following additional variables:

| Column | Description |
|---|---|
| `non_HDL_C` | Non-HDL cholesterol, calculated if missing |
| `LDL_C_Friedewald` | LDL-C estimated using the Friedewald equation |
| `MartinHopkins_factor` | Adjustable triglyceride divisor used in the Martin-Hopkins equation |
| `LDL_C_MartinHopkins` | LDL-C estimated using the Martin-Hopkins equation |
| `LDL_C_Sampson` | LDL-C estimated using the Sampson equation |

## Example

```r
library(dplyr)

source("R/ldl_equations.R")

example_data <- tibble::tibble(
  ID = 1:5,
  TG = c(80, 120, 180, 260, 410),
  TC = c(180, 160, 140, 210, 260),
  HDL_C = c(50, 45, 40, 35, 38)
)

example_results <- calculate_ldl_equations(example_data)

print(example_results)
