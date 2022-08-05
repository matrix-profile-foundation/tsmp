# nolint start
options(repos = c(CRAN = "https://cran.rstudio.org"))
if (.Platform$OS.type == "windows") {
  Sys.setenv(LC_CTYPE = "C")
}
source("renv/activate.R")

if (Sys.getenv("CI") == "") { # not CI

  if (Sys.getenv("RSTUDIO") == "") { # not in RStudio
    if (interactive()) {
      options(
        warnPartialMatchArgs = FALSE,
        warnPartialMatchDollar = FALSE,
        warnPartialMatchAttr = FALSE,
        usethis.protocol = "https",
        warn = 1 # warnings appear immediately, not in the end
        # error = recover
      )
      options(
        vsc.rstudioapi = TRUE,
        max.print = 1000,
        width = 200,
        vsc.show_object_size = TRUE,
        vsc.globalenv = TRUE,
        vsc.dev.args = list(width = 1000, height = 700)
      )

      options(languageserver.formatting_style = function(options) {
        style <- styler::tidyverse_style(scope = "tokens", indent_by = 2)
        style
      })

      suppressMessages(
        suppressWarnings({
          require("testthat", quietly = TRUE)
          require("devtools", quietly = TRUE)
          require("usethis", quietly = TRUE)
          require("conflicted", quietly = TRUE)
          require("here", quietly = TRUE)
          require("glue", quietly = TRUE)
        })
      )

      conflicted::conflict_prefer("filter", "dplyr")
      options(dplyr.summarise.inform = FALSE)

      if (.Platform$OS.type != "windows") {
        if (suppressMessages(requireNamespace("prettycode", quietly = TRUE))) {
          suppressMessages(prettycode::prettycode())
        }
      }
      if (Sys.getenv("RADIAN_VERSION") == "") {
        loadhistory() # if no file, no problem.

        # Cleaning up function
        .Last <- function() {
          savehistory() # comment this line if you don't want to save history
          cat("bye bye...\n") # print this so we see if any non-interactive session is lost here
        }
      }
    }
  } else { # in RStudio
        suppressMessages(
      suppressWarnings({
        require("testthat", quietly = TRUE)
        require("devtools", quietly = TRUE)
        require("usethis", quietly = TRUE)
        require("conflicted", quietly = TRUE)
        require("here", quietly = TRUE)
        require("glue", quietly = TRUE)
      })
    )

    conflicted::conflict_prefer("filter", "dplyr")
    options(dplyr.summarise.inform = FALSE)
  }
} else { # is CI
  message("Running .RProfile in CI")
}
# nolint end
