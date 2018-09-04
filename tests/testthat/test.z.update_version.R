context("Increase Package version")

if (skip_on_cran() && skip_on_travis()) {
  get_package_version <- function(package_location = ".") {
    ## Read DESCRIPTION file
    desc <- readLines(file.path(package_location, "DESCRIPTION"))

    ## Find the line where the version is defined
    version_line <- grep("^Version\\:", desc)

    ## Extract version number
    current_version <- gsub("^Version\\:\\s*", "", desc[version_line])

    ## Return the current number
    return(current_version)
  }

  update_package_version <- function(package_location = ".") {
    ## Read DESCRIPTION file
    desc <- readLines(file.path(package_location, "DESCRIPTION"))

    ## Find the line where the version is defined
    version_line <- grep("^Version\\:", desc)

    ## Extract version number
    current_version <- gsub("^Version\\:\\s*", "", desc[version_line])
    #|
    ## Split the version number into two; a piece to keep, a piece to increment
    version_number <- strsplit(current_version, "\\.")[[1]]
    version_parts <- length(version_number)
    version_number_keep <- paste(version_number[1:(version_parts - 1)], sep = "", collapse = ".")
    version_number_update <- version_number[version_parts]

    ## Replace old version number with new one (increment by 1)
    old_version <- as.numeric(version_number_update)
    new_version <- old_version + 1

    ## Build final version number
    version_final <- paste(version_number_keep, new_version, sep = ".")

    ## Update DESCRIPTION file (in R)
    desc[version_line] <- paste0("Version: ", version_final)

    ## Update the actual DESCRIPTION file
    writeLines(desc, file.path(package_location, "DESCRIPTION"))

    ## Return the updated version number to screen
    return(version_final)
  }

  curr <- get_package_version(Sys.getenv("R_PACKRAT_PROJECT_DIR"))
  new <- update_package_version(Sys.getenv("R_PACKRAT_PROJECT_DIR"))

  test_that("New version is not old version", {
    expect_false(isTRUE(all.equal(new, curr)))
  })
}
