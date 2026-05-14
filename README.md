# LDL-C estimation equations
# Author: Jordi Tortosa-Carreres
#
# This script calculates LDL-C using:
# 1. Friedewald equation
# 2. Martin-Hopkins equation
# 3. Sampson equation
#
# Required input variables:
# - TC: total cholesterol, mg/dL
# - HDL_C: HDL cholesterol, mg/dL
# - TG: triglycerides, mg/dL
#
# Optional input variable:
# - non_HDL_C: non-HDL cholesterol, mg/dL
#   If not provided or missing, it is calculated as TC - HDL_C.
#
# Output:
# The original dataset plus:
# - non_HDL_C
# - LDL_C_Friedewald
# - MartinHopkins_factor
# - LDL_C_MartinHopkins
# - LDL_C_Sampson

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)

# Martin-Hopkins adjustable factor table
martin_hopkins_table <- read_csv("
TG_range,<100,100-129,130-159,160-189,190-219,>=220
7-44,3.3,3.2,3.1,3.0,2.9,2.7
45-49,3.8,3.5,3.5,3.4,3.3,3.2
50-53,4.0,3.8,3.6,3.6,3.5,3.3
54-56,4.2,3.9,3.9,3.7,3.6,3.4
57-59,4.2,4.1,3.9,3.8,3.7,3.6
60-61,4.4,4.1,4.0,4.0,3.8,3.8
62-64,4.5,4.2,4.1,4.0,3.9,3.7
65-66,4.6,4.3,4.1,4.1,4.1,3.9
67-68,4.5,4.5,4.2,4.2,3.9,3.8
69-70,4.7,4.4,4.3,4.1,4.1,3.8
71-72,4.7,4.5,4.2,4.2,4.2,4.0
73-75,4.9,4.6,4.4,4.3,4.2,4.1
76-77,4.8,4.5,4.5,4.3,4.2,4.1
78-79,4.9,4.6,4.4,4.3,4.3,4.2
80-81,5.0,4.7,4.5,4.4,4.3,4.2
82-83,5.1,4.8,4.6,4.4,4.4,4.2
84-85,5.0,4.7,4.7,4.5,4.4,4.3
86-87,5.1,4.8,4.6,4.5,4.5,4.3
88-89,5.2,4.9,4.7,4.4,4.4,4.2
90-91,5.3,5.0,4.7,4.6,4.5,4.3
92-93,5.2,4.9,4.8,4.6,4.4,4.4
94-95,5.3,5.0,4.8,4.5,4.5,4.3
96-97,5.3,5.1,4.8,4.6,4.6,4.4
98-100,5.4,5.2,4.9,4.7,4.5,4.3
101-102,5.6,5.1,4.9,4.6,4.6,4.4
103-104,5.5,5.2,5.0,4.7,4.7,4.5
105-106,5.5,5.3,5.0,4.8,4.6,4.6
107-109,5.6,5.3,5.0,4.7,4.7,4.5
110-111,5.8,5.3,5.0,4.8,4.6,4.6
112-114,5.7,5.4,5.1,4.9,4.7,4.5
115-116,5.8,5.5,5.2,4.8,4.8,4.6
117-119,5.9,5.4,5.1,5.0,4.8,4.6
120-121,6.0,5.5,5.2,5.0,4.8,4.6
122-124,5.9,5.6,5.3,5.0,4.8,4.6
125-127,6.0,5.7,5.3,5.0,4.8,4.6
128-130,6.1,5.7,5.4,5.1,4.9,4.7
131-133,6.0,5.7,5.3,5.1,4.9,4.7
134-136,6.2,5.8,5.4,5.2,5.0,4.7
137-140,6.3,5.8,5.5,5.2,5.0,5.0
141-143,6.2,5.9,5.5,5.3,5.1,4.9
144-147,6.4,5.9,5.6,5.3,5.1,4.8
148-151,6.5,6.0,5.6,5.4,5.1,4.8
152-155,6.6,6.1,5.7,5.4,5.2,4.9
156-159,6.6,6.1,5.8,5.4,5.2,4.9
160-164,6.7,6.2,5.8,5.5,5.2,4.9
165-169,6.8,6.3,5.9,5.5,5.3,5.0
170-174,6.9,6.4,6.0,5.6,5.3,5.0
175-180,7.0,6.4,6.0,5.6,5.4,5.0
181-187,7.1,6.5,6.1,5.8,5.4,5.1
188-194,7.3,6.6,6.2,5.8,5.5,5.2
195-202,7.4,6.7,6.3,5.9,5.6,5.2
203-211,7.5,6.8,6.4,6.0,5.6,5.3
212-221,7.6,7.0,6.5,6.1,5.7,5.4
222-233,7.9,7.2,6.6,6.2,5.8,5.4
234-248,8.1,7.3,6.7,6.3,5.9,5.5
249-267,8.4,7.5,6.9,6.4,6.0,5.6
268-292,8.6,7.7,7.1,6.6,6.2,5.6
293-330,9.2,8.1,7.3,6.8,6.3,5.8
331-399,9.9,8.6,7.8,7.2,6.7,6.0
400-13975,11.9,10.0,8.8,8.1,7.5,6.7
", show_col_types = FALSE)

martin_hopkins_table_long <- martin_hopkins_table %>%
  separate(TG_range, into = c("TG_min", "TG_max"), sep = "-", convert = TRUE) %>%
  pivot_longer(
    cols = -c(TG_min, TG_max),
    names_to = "non_HDL_C_range",
    values_to = "MartinHopkins_factor"
  ) %>%
  mutate(
    non_HDL_C_min = case_when(
      non_HDL_C_range == "<100" ~ 0,
      non_HDL_C_range == ">=220" ~ 220,
      TRUE ~ as.numeric(str_extract(non_HDL_C_range, "^[0-9]+"))
    ),
    non_HDL_C_max = case_when(
      non_HDL_C_range == "<100" ~ 99,
      non_HDL_C_range == ">=220" ~ Inf,
      TRUE ~ as.numeric(str_extract(non_HDL_C_range, "[0-9]+$"))
    )
  ) %>%
  select(TG_min, TG_max, non_HDL_C_min, non_HDL_C_max, MartinHopkins_factor)

get_martin_hopkins_factor <- function(TG, non_HDL_C) {
  pmap_dbl(
    list(TG, non_HDL_C),
    function(tg_i, non_hdl_i) {
      if (is.na(tg_i) || is.na(non_hdl_i)) {
        return(NA_real_)
      }

      factor_value <- martin_hopkins_table_long %>%
        filter(
          tg_i >= TG_min,
          tg_i <= TG_max,
          non_hdl_i >= non_HDL_C_min,
          non_hdl_i <= non_HDL_C_max
        ) %>%
        pull(MartinHopkins_factor)

      if (length(factor_value) == 0) {
        return(NA_real_)
      }

      factor_value[1]
    }
  )
}

calculate_ldl_equations <- function(data) {
  required_columns <- c("TC", "HDL_C", "TG")

  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "The dataset must contain the following columns: ",
        paste(required_columns, collapse = ", "),
        ". Missing: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  if (!"non_HDL_C" %in% names(data)) {
    data <- data %>%
      mutate(non_HDL_C = TC - HDL_C)
  }

  data %>%
    mutate(
      non_HDL_C = if_else(is.na(non_HDL_C), TC - HDL_C, non_HDL_C),

      LDL_C_Friedewald = TC - HDL_C - (TG / 5),

      MartinHopkins_factor = get_martin_hopkins_factor(TG, non_HDL_C),

      LDL_C_MartinHopkins = TC - HDL_C - (TG / MartinHopkins_factor),

      LDL_C_Sampson =
        (TC / 0.948) -
        (HDL_C / 0.971) -
        (
          (TG / 8.56) +
          ((TG * non_HDL_C) / 2140) -
          ((TG^2) / 16100)
        ) -
        9.44
    )
}

# Example:

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
