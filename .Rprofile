# nolint start
if (Sys.getenv("CI") == "") {
  options(repos = c(CRAN = "https://cran.rstudio.org"))

  if (.Platform$OS.type == "windows") {
    Sys.setenv(LC_CTYPE = "C")
  }

  # options(error = renv:::renv_error_handler_call())
  # options(warn = 2L)

  # Activate renv for all kind of environments
  source("renv/activate.R")
  cat("Renv activated from .Rprofile\n")

  # not CI
  if (interactive()) {
    if (Sys.getenv("RSTUDIO") == "") {
      # interactive and not RSTUDIO ENV
      # This is the default terminal environment
      # cat("INTERACTIVE; no RSTUDIO\n")
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

      vscode_init <- file.path(Sys.getenv("HOME"), ".vscode-R", "init.R")
      if (interactive() && file.exists(vscode_init)) {
        source(vscode_init)
        if (exists(".First.sys", envir = globalenv(), inherits = FALSE)) {
          .First.sys()
        }
      }


      # if httpgd is installed, let's use it
      # This breaks rendering video
      # if ("httpgd" %in% .packages(all.available = TRUE)) {
      #   options(vsc.plot = FALSE)
      #   options(device = function(...) {
      #     httpgd::hgd(silent = TRUE)
      #     .vsc.browser(httpgd::hgd_url(history = FALSE), viewer = "Beside")
      #   })
      # }

      suppressMessages(
        suppressWarnings({
          require("testthat", quietly = TRUE)
          require("devtools", quietly = TRUE)
          require("usethis", quietly = TRUE)
          if (require("conflicted", quietly = TRUE)) {
            suppressMessages({
              conflicted::conflict_prefer("filter", "dplyr")
              conflicted::conflict_prefer("box", "shinydashboard")
              conflicted::conflict_prefer("notificationItem", "shinydashboard")
            })
          }
          require("here", quietly = TRUE)
          require("glue", quietly = TRUE)
        })
      )

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
    } else {
      # interactive and RSTUDIO ENV
      # This is supposed to be the RSTUDIO terminal
      # cat("INTERACTIVE; RSTUDIO\n")
      # is RSTUDIO
      suppressMessages(
        suppressWarnings({
          require("testthat", quietly = TRUE)
          require("devtools", quietly = TRUE)
          require("usethis", quietly = TRUE)
          if (require("conflicted", quietly = TRUE)) {
            conflicted::conflict_prefer("filter", "dplyr")
          }
          require("here", quietly = TRUE)
          require("glue", quietly = TRUE)
        })
      )

      options(dplyr.summarise.inform = FALSE)
    }
  } else {
    if (Sys.getenv("RSTUDIO") == "") {
      # non-interactive and not RSTUDIO ENV
      # This is the default non-interactive environment
      cat("NON-INTERACTIVE; no RSTUDIO\n")
      invisible(NULL)
    } else {
      # non-interactive and RSTUDIO ENV
      # This is the supposed to be the background RSTUDIO environment
      cat("NON-INTERACTIVE; RSTUDIO\n")
      invisible(NULL)
    }
  }
} else {
  cat("CI ENV\n")
  invisible(NULL)
  # is CI
  # suppressMessages(
  #   suppressWarnings({
  #     require("here", quietly = TRUE)
  #     require("workflowr", quietly = TRUE)
  #     require("targets", quietly = TRUE)
  #     require("tarchetypes", quietly = TRUE)
  #     require("gittargets", quietly = TRUE)
  #   })
  # )
}

# nolint end
