#' Title
#'
#' @param comp component to analyze
#' @param model input model
#' @param combine_stats whether to combine statistics
#'
#' @return a tibble
#' @export
#'
#' @importFrom broom.mixed tidy
format_ranef_table = function(model, comp = "cond", combine_stats = TRUE)
{
  to_return = model %>%
    broom.mixed::tidy(effects="ran_vals") %>%
    filter(if(!is.null(comp)){component  == comp} else {TRUE}, group == "Participant") %>%
    mutate(
      pval = 2 * pnorm(abs(estimate), sd = std.error, lower = FALSE),
      CIlow = estimate - 1.96 * std.error,
      CIhigh = estimate + 1.96 * std.error,
      across(all_of(c("estimate", "std.error", "pval", "CIlow", "CIhigh")), formatC, format = "f", digits = 3),
      pval = if_else(pval == "0.000", "<0.001", pval),
      all = paste(sep = "", estimate, " (P: ", pval, "; 95% CI: ", CIlow, ", ", CIhigh, ")")) %>%
    rename(`Participant #` = level)

  if(combine_stats)
  {
    to_return %<>% pivot_wider(
      id_cols = "Participant #",
      names_from = "term",
      values_from = "all")
  } else
  {

    to_return %<>% select(
      `Participant #`,
      Parameter  = term,
      Estimate = estimate,
      `Std. Error` = std.error,
      `Pr(>|z|)` = pval,
      `2.5 %` = CIlow,
      `97.5 %` = CIhigh)

  }

  return(to_return)
}
