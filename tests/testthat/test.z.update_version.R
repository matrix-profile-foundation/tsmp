context("Increase Package version")

if (skip_on_cran() && skip_on_travis() && skip_on_appveyor()) {
  get.package.version <- function(package.location = ".") {
    ## Read DESCRIPTION file
    desc <- readLines(file.path(package.location, "DESCRIPTION"))

    ## Find the line where the version is defined
    version.line <- grep("^Version\\:", desc)

    ## Extract version number
    current.version <- gsub("^Version\\:\\s*", "", desc[version.line])

    ## Return the current number
    return(current.version)
  }

  update.package.version <- function(package.location = ".") {
    ## Read DESCRIPTION file
    desc <- readLines(file.path(package.location, "DESCRIPTION"))

    ## Find the line where the version is defined
    version.line <- grep("^Version\\:", desc)

    ## Extract version number
    current.version <- gsub("^Version\\:\\s*", "", desc[version.line])
    #|
    ## Split the version number into two; a piece to keep, a piece to increment
    version.number <- strsplit(current.version, "\\.")[[1]]
    version.parts <- length(version.number)
    version.number.keep <- paste(version.number[1:(version.parts - 1)], sep = "", collapse = ".")
    version.number.update <- version.number[version.parts]

    ## Replace old version number with new one (increment by 1)
    old.version <- as.numeric(version.number.update)
    new.version <- old.version + 1

    ## Build final version number
    version.final <- paste(version.number.keep, new.version, sep = ".")

    ## Update DESCRIPTION file (in R)
    desc[version.line] <- paste0("Version: ", version.final)

    ## Update the actual DESCRIPTION file
    writeLines(desc, file.path(package.location, "DESCRIPTION"))

    ## Return the updated version number to screen
    return(version.final)
  }

  curr <- get.package.version(Sys.getenv("R_PACKRAT_PROJECT_DIR"))
  new <- update.package.version(Sys.getenv("R_PACKRAT_PROJECT_DIR"))

  test_that("New version is not old version", {
    expect_false(isTRUE(all.equal(new, curr)))
  })
}
