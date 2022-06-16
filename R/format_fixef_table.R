#' Title
#'
#' @param data
#'
#' @return
#' @export
#' @importFrom magrittr %>% %<>%
#' @import dplyr


format_fixef_table = function(model, type = "cond")
{
  coef(summary(fit4))[[type]] %>%
    as_tibble() %>%
    select(-`z value`) %>%
    bind_cols(
      confint(model, method = 'Wald', oldNames=FALSE)[names(fixef(model))[[type]],] %>% as_tibble()
    ) %>%
    mutate(
      `Variable name` = names(fixef(model))[[type]],
      across(all_of(c("Estimate", "Std. Error", "Pr(>|z|)", "2.5 %", "97.5 %")), formatC, format = "f", digits = 3),
      `Pr(>|z|)` = if_else(`Pr(>|z|)` == "0.000", "<0.001", `Pr(>|z|)`)
    ) %>%
    relocate(`Variable name`, .before = everything())
}
