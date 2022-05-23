library(readr)
library(dplyr)
library(tidyr)
region_durations = read_csv("inst/extdata/Duration of brushing of each dental region.csv")


session_durations = read_csv("inst/extdata/Non effective brushing duration.csv")
colnames(session_durations) = c("Participant", "Session", "Session duration (samples)")

library(magrittr)
library(dplyr)
session_durations %<>% group_by(
  Participant) %>%
  mutate(
    Session = 1:n(),
    `Session duration (seconds)` = `Session duration (samples)` / 25,
    n_samples = `Session duration (samples)`
) %>%
  relocate(
    Participant, Session
  )

usethis::use_data(session_durations, overwrite = TRUE)

## rearrange data

to_pivot = grep(value = TRUE, "(Man)|(Max)", colnames(region_durations))
data2 = region_durations %>%
  pivot_longer(cols = all_of(to_pivot), names_to = "Region", values_to = "n_samples")

code1 = c("L" = "Lingual", "O" = "Occlusal", "B" = "Buccal")
strings1 = substr(data2$Region, 4,4)
code2 = c("R" = "Right", "L" = "Left", "A" = "Anterior")
data2 %<>%
  mutate(
    Surface = code1[substr(Region, 5,5)],
    Side = code2[substr(Region, 4,4)],
    Side = relevel(factor(Side), ref = "Right"),
    Jaw = relevel(factor(if_else(substr(Region, 1,3) == "Max", "Maxillar", "Mandibular")), ref = 'Mandibular'),
    `Duration (seconds)` = n_samples / 25
  )

session_durations2 = data2 %>%
  group_by(Participant, Session) %>%
  dplyr::summarize(
    .groups = "drop",
    `Session duration (samples)` = sum(n_samples)) %>%
  mutate(
    n_samples = `Session duration (samples)`,
    `Session duration (seconds)` = `Session duration (samples)` / 25
  )

data2 %<>%
  left_join(
    session_durations2 %>% select(-n_samples),
    by = c("Participant", "Session")
  )

region_durations = data2

usethis::use_data(region_durations, overwrite = TRUE)

usethis::use_data(session_durations2, overwrite = TRUE)



pressure_durations =
  read_csv("inst/extdata/Duration of brushing with excessive pressure on each dental region.csv") %>%
  pivot_longer(cols = all_of(to_pivot), names_to = "Region", values_to = "n_samples") %>%
  mutate(
    Surface = code1[substr(Region, 5,5)],
    Side = code2[substr(Region, 4,4)],
    Side = relevel(factor(Side), ref = "Right"),
    Jaw = relevel(factor(if_else(substr(Region, 1,3) == "Max", "Maxillar", "Mandibular")), ref = 'Mandibular'),
    `Duration (seconds)` = n_samples / 25
  ) %>%
  left_join(
    session_durations2 %>% select(-n_samples),
    by = c("Participant", "Session")
  )
  usethis::use_data(pressure_durations, overwrite = TRUE)
