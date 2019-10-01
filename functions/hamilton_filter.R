hamilton_filter <- function(x, h = 8, p = 4, ...) {
  neverHP <- neverhpfilter::yth_glm(x = x, h = h, p = p, ...)
  trend <- xts::as.xts(unname(neverHP$fitted.values), order.by = 
                         get(paste0("as.", class(index(x))))(names(neverHP$fitted.values)))
  names(trend) <- paste0(names(x), ".trend")
  cycle <- xts::as.xts(unname(neverHP$residuals), order.by = 
                         get(paste0("as.", class(index(x))))(names(neverHP$residuals)))
  names(cycle) <- paste0(names(x), ".cycle")
  random <- x - stats::lag(x, k = h, na.pad = TRUE)
  names(random) <- paste0(names(x), ".random")
  all <- merge(x, trend, cycle, random)
  names(all) <- c(names(x), paste0(names(x), ".", c("trend", "cycle", "random")))
  ret_list <- list(glm_result = neverHP, ts = all)
  return(ret_list)
}
