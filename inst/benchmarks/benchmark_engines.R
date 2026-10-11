# ==============================================================================
# Benchmark Suite: Motores Nativos C++ (matrixprofiler) vs Legado (tsmp)
# ==============================================================================

library(matrixprofiler)
library(tsmp)

cat("==================================================================\n")
cat(" BENCHMARK 1: mSTOMP (Multidimensional STOMP: C++ vs Legado)\n")
cat("==================================================================\n")

# Gerando dados sintéticos multivariados (3 dimensões, 1.500 pontos)
set.seed(42)
n_obs <- 1500
n_dims <- 3
w <- 50
mat_data <- matrix(rnorm(n_obs * n_dims), nrow = n_obs, ncol = n_dims)

# 1. mSTOMP nativo C++ (matrixprofiler::mstomp)
t_native_start <- Sys.time()
res_native <- matrixprofiler::mstomp(mat_data, window_size = w, n_workers = 1, progress = FALSE)
t_native <- as.numeric(difftime(Sys.time(), t_native_start, units = "secs"))

cat(sprintf("mSTOMP Nativo C++:  %.4f s\n", t_native))
cat(sprintf("Velocidade do motor: %.1f iterações/seg\n", (n_obs - w + 1) / t_native))

cat("\n==================================================================\n")
cat(" BENCHMARK 2: MPdist (C++ com Monotonic Deque vs R)\n")
cat("==================================================================\n")

ts_ref <- rnorm(2000)
ts_query <- rnorm(2000)
w_mpdist <- 50

t_mpdist_start <- Sys.time()
d_mpdist <- matrixprofiler::mpdist(ts_ref, ts_query, window_size = w_mpdist, n_workers = 1)
t_mpdist <- as.numeric(difftime(Sys.time(), t_mpdist_start, units = "secs"))

cat(sprintf("MPdist Nativo C++:   %.4f s (dist = %.4f)\n", t_mpdist, d_mpdist))

cat("\n==================================================================\n")
cat(" BENCHMARK 3: MASS Multithreading (1 vs 4 Workers)\n")
cat("==================================================================\n")

n_large <- 30000
ts_large <- rnorm(n_large)
w_mass <- 100

pre_large <- matrixprofiler::mass_pre(ts_large, window_size = w_mass, type = "normalized")

# v2 serial vs paralelo
t_v2_1_start <- Sys.time()
res_v2_1 <- matrixprofiler::mass(pre_large, ts_large, index = 1, version = "v2", n_workers = 1)
t_v2_1 <- as.numeric(difftime(Sys.time(), t_v2_1_start, units = "secs"))

t_v2_4_start <- Sys.time()
res_v2_4 <- matrixprofiler::mass(pre_large, ts_large, index = 1, version = "v2", n_workers = 4)
t_v2_4 <- as.numeric(difftime(Sys.time(), t_v2_4_start, units = "secs"))

cat(sprintf("MASS v2 (1 worker):  %.4f s\n", t_v2_1))
cat(sprintf("MASS v2 (4 workers): %.4f s (speedup: %.2fx)\n", t_v2_4, t_v2_1 / max(1e-6, t_v2_4)))

# v3 serial vs paralelo
t_v3_1_start <- Sys.time()
res_v3_1 <- matrixprofiler::mass(pre_large, ts_large, index = 1, version = "v3", n_workers = 1)
t_v3_1 <- as.numeric(difftime(Sys.time(), t_v3_1_start, units = "secs"))

t_v3_4_start <- Sys.time()
res_v3_4 <- matrixprofiler::mass(pre_large, ts_large, index = 1, version = "v3", n_workers = 4)
t_v3_4 <- as.numeric(difftime(Sys.time(), t_v3_4_start, units = "secs"))

cat(sprintf("MASS v3 (1 worker):  %.4f s\n", t_v3_1))
cat(sprintf("MASS v3 (4 workers): %.4f s (speedup: %.2fx)\n", t_v3_4, t_v3_1 / max(1e-6, t_v3_4)))

cat("\n==================================================================\n")
cat(" BENCHMARK SUITE FINALIZADA COM SUCESSO!\n")
cat("==================================================================\n")
