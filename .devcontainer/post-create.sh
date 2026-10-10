#!/usr/bin/env bash
set -euo pipefail

echo "[devcontainer] Ensuring renv is installed"
R -q -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")'

echo "[devcontainer] Checking whether renv is already active and ensuring renv infrastructure"
Rscript - <<'RSCRIPT'
activated <- FALSE
if (requireNamespace("renv", quietly = TRUE)) {
	if (exists("activated", envir = asNamespace("renv"), inherits = FALSE)) {
		activated <- renv::activated()
	} else {
		# fallback: check env var or library path for renv
		activated <- nzchar(Sys.getenv("RENV_PROJECT")) || any(grepl("renv/library", .libPaths()))
	}
}

if (activated) {
	message("[devcontainer] renv already active in this R session")
} else {
	need_init <- ( !file.exists("renv.lock") || !dir.exists("renv") || !file.exists(file.path("renv","activate.R")) )
	if (need_init) {
		message("[devcontainer] renv infrastructure missing — initializing (bare = TRUE)")
		if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")
		renv::init(bare = TRUE)
	} else {
		message("[devcontainer] renv infrastructure present but not active; activating by sourcing renv/activate.R")
		tryCatch(source(file.path("renv","activate.R")), error = function(e) message("[devcontainer] failed to source renv/activate.R: ", e$message))
	}
}
RSCRIPT

echo "[devcontainer] Installing prerequisite CRAN packages: languageserver, prompt, etc"
R -q -e 'renv::install("languageserver")'
R -q -e 'renv::install(c("prompt", "devtools", "conflicted", "here", "codemetar", "jsonld"))'
R -q -e 'renv::install(c("ggnetwork", "miniCRAN", "r-lib/pkgapi", "r-lib/revdepcheck", "rhub"))'

echo "[devcontainer] Installing nx10/httpgd via renv"
R -q -e 'renv::install("nx10/httpgd")'

echo "[devcontainer] Running renv::install() to install project dependencies"
R -q -e 'renv::install()'

echo "[devcontainer] Configuring Git for SSH commit signing"
git config --global user.name "Francisco Bischoff"
git config --global user.email "franzbischoff@gmail.com"
git config --global gpg.format ssh
git config --global commit.gpgsign true

# Get the SSH key from the forwarded agent and set it as signing key
SSH_KEY=$(ssh-add -L 2>/dev/null | head -1)
if [ -n "$SSH_KEY" ]; then
  git config --global user.signingkey "$SSH_KEY"
  echo "[devcontainer] SSH signing key configured: ${SSH_KEY:0:50}..."
else
  echo "[devcontainer] WARNING: No SSH key found in agent. Signing may not work."
fi

echo "[devcontainer] post-create steps completed"
