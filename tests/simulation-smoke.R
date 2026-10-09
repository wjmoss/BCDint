# Run the real simulation driver over all 24 configurations in an isolated copy.
# Run from the repository root: Rscript --vanilla tests/simulation-smoke.R
run_smoke_test <- function() {
  repo <- normalizePath(".")
  scripts <- c("simulation.R", "generateModel.R", "ricf_int.R", "ricf_dg.R")
  stopifnot(all(file.exists(file.path(repo, scripts))))

  work <- tempfile("bcd-simulation-smoke-")
  dir.create(work)
  previous_repl <- Sys.getenv("BCD_REPL", unset=NA_character_)
  passed <- FALSE
  on.exit({
    setwd(repo)
    if (is.na(previous_repl)) Sys.unsetenv("BCD_REPL") else
      Sys.setenv(BCD_REPL=previous_repl)
    if (passed) unlink(work, recursive=TRUE)
  }, add=TRUE)

  stopifnot(all(file.copy(file.path(repo, scripts), work)))
  setwd(work)
  stopifnot(!dir.exists("data"))
  replicates <- 2L
  Sys.setenv(BCD_REPL=replicates)
  executable <- file.path(R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(executable, c("--vanilla", "simulation.R"),
    stdout="simulation.log", stderr="simulation-errors.log")
  check <- function(condition, message) {
    if (!isTRUE(condition))
      stop(paste0(message, "; outputs retained in ", work), call.=FALSE)
  }
  check(status == 0L, "Simulation driver failed")
  check(dir.exists("data"), "Simulation did not create data/")

  grid <- expand.grid(v=c(5,10), n=c(5,10), l=c(0,3,4), d=c(.2,.3))
  expected <- with(grid, paste0("v=", v, "n=", n, "l=", l,
                                "d=", d, "seed=1.rda"))
  actual <- list.files("data", pattern="[.]rda$")
  check(identical(sort(actual), sort(expected)),
        "Expected exactly 24 configuration archives")

  for (file in expected) {
    saved <- new.env(parent=emptyenv())
    load(file.path("data", file), envir=saved)
    check(all(c("res", "models", "res_int", "res_dg_agg") %in% ls(saved)),
          paste(file, "has missing saved objects"))
    check(identical(dim(saved$res), c(2L,4L)) && all(is.finite(saved$res)),
          paste(file, "has invalid summary metrics"))
    check(all(saved$res[,1] == replicates) && all(saved$res[,4] >= 0),
          paste(file, "did not converge or has invalid timing"))
    for (name in c("models", "res_int", "res_dg_agg"))
      check(length(saved[[name]]) == replicates,
            paste(file, "has the wrong replicate count in", name))
    for (fit in c(saved$res_int, saved$res_dg_agg))
      check(all(is.finite(fit$Lambdahat)) &&
            all(is.finite(fit$Omegahat)) && all(fit$Omegahat > 0) &&
            all(is.finite(fit$Sigmahat)), paste(file, "has invalid estimates"))
  }
  passed <- TRUE
  cat("Passed all 24 simulation configurations with 2 replicates each.\n")
  cat("Temporary outputs removed; repository data/ was not used.\n")
}

run_smoke_test()
