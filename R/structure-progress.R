# One owner per public call; callbacks let callers own their display instead.
.structure_progress <- function(progress, env = parent.frame()) {
  if (is.function(progress)) {
    return(progress)
  }
  checkmate::assert_flag(progress)
  if (!progress) {
    return(NULL)
  }
  id <- cli::cli_progress_bar(
    name = "Constructing structures",
    total = NA,
    auto_terminate = FALSE,
    format = "{cli::pb_name}: {cli::pb_status} {cli::pb_current}/{cli::pb_total}",
    .auto_close = TRUE,
    .envir = env
  )
  function(stage, current, total) {
    cli::cli_progress_update(
      id = id,
      set = current,
      total = total,
      status = stage
    )
  }
}

.structure_progress_stage <- function(progress, stage, total) {
  if (is.null(progress)) {
    return(NULL)
  }
  progress(stage, 0L, total)
  last <- proc.time()[[3L]]
  function(current) {
    now <- proc.time()[[3L]]
    if (current == total || now - last >= 0.1) {
      progress(stage, current, total)
      last <<- now
    }
    invisible(NULL)
  }
}

.structure_progress_map <- function(x, .f, progress, stage, map = lapply) {
  use_purrr <- identical(map, purrr::map)
  if (is.null(progress)) {
    if (use_purrr) {
      return(purrr::map(x, .f))
    }
    return(lapply(x, .f))
  }
  update <- .structure_progress_stage(progress, stage, length(x))
  i <- 0L
  wrapped <- function(element) {
    result <- .f(element)
    i <<- i + 1L
    update(i)
    result
  }
  if (use_purrr) {
    return(purrr::map(x, wrapped))
  }
  lapply(x, wrapped)
}
