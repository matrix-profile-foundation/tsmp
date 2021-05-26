# nolint start
source("renv/activate.R")

if (Sys.getenv("CI") == "") { # not CI

  if (Sys.getenv("RSTUDIO") == "") {
    suppressMessages(if (requireNamespace("languageserver")) {
      a <- try(suppressWarnings(source(file.path(
        Sys.getenv(if (.Platform$OS.type == "windows") {
          "USERPROFILE"
        } else {
          "HOME"
        }),
        ".vscode-R",
        "init.R"
      ))),
      silent = TRUE
      )

      rm(a)
    })
  }

  if (interactive()) {
    options(
      warnPartialMatchArgs = FALSE,
      warnPartialMatchDollar = FALSE,
      warnPartialMatchAttr = FALSE,
      usethis.protocol = "https",
      vsc.rstudioapi = TRUE
      # error = recover
    )

    suppressMessages(
      suppressWarnings({
        require("testthat", quietly = TRUE)
        require("devtools", quietly = TRUE)
        require("usethis", quietly = TRUE)
        require("conflicted", quietly = TRUE)
        require("here", quietly = TRUE)
      })
    )
    # suppressMessages(prettycode::prettycode())

    if (suppressMessages(requireNamespace("prompt", quietly = TRUE))) {
      prompt::set_prompt(function(...) {
        paste0(
          "[",
          prompt::git_branch(),
          prompt::git_dirty(),
          prompt::git_arrows(),
          "] ",
          prompt::prompt_runtime()
        )
      })
    }
  }
} else { # is CI
  message("Running .RProfile in CI")
}
# nolint end
