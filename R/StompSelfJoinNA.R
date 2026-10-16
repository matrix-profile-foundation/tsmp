StompSelfJoinNA <- function(time_series, window_size, exclusion_zone) {
#' This is the code for the paper "Admissible Motif Discovery with Missing
#' Data" by Yan Zhu, Abdullah Mueen and Eamonn Keogh.
#'
#' Input:
#'    time_series: Time Series with missing values
#'    window_size: Subsequence length, window size for motif discovery
#' Output:
#'   Matrix Profile: The lowerbound Matrix Profile
#'   MPindex: The lowerbound Matrix Profile Index
#'
#' Translated from Matlab to R by Christoper English
#' All mistakes are mine as their code worked and results seemed reliable
#' as executed in Octave
#' @param time_series 'vector' or 'matrix' time series with missing values
#' @param window_size 'int' sliding window size
#' @param exclusion_zone 'numeric' to prevent trivial matches tsmp default 1/2 Seismology example uses 1/4 (0.25)
#' @return Returns 'MatrixProfile', 'MPindex', 'Z' (time_series with NA replaced by 0)
# mode:1-FFT 2-no FFT only mode 1 implemented
mode=1;
time_series = time_series # for debugging
window_size = window_size # for debugging
# set trivial match exclusion zone
exclusion_zone = exclusion_zone 
exclusion_zone = round(window_size * exclusion_zone)

# check input
if (window_size > length(time_series)/2) {
error('Error: Time series is too short relative to desired window size')
}
if (window_size < 4) {
error('Error: window_size must be at least 4L')
}
num = length(time_series)
MatrixProfileLen = num - window_size +1
Vmax = rep.int(0, times = MatrixProfileLen)
Vmin = rep.int(0, times = MatrixProfileLen)

#' create Z, B, and X from time_series
#' switching out NA for 0 insulates for use with tsmp
#' might need to DO THE SAME for is.infinite for purposes of generality
Z = time_series
Z[is.na(Z)] = 0
for (s in 1:MatrixProfileLen) {   
Vmax[s] = max(Z[s]:Z[s+window_size -1], na.rm = FALSE)
Vmin[s] = min(Z[s]:Z[s+window_size -1], na.rm = FALSE)
}
B <- rep.int(1, times = length(time_series))
B[is.na(time_series)] = 0
X = Z^2
#
MatrixProfile = matrix(0, MatrixProfileLen,1) #from rep.int vec to matrix
MPindex = rep.int(0, times = MatrixProfileLen) # why 0 here and 1 at line 137
# Attempt to include progress bar
pb <- progress::progress_bar$new(
   format = 'StompNA [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta', clear = FALSE, total=MatrixProfileLen, width = 80)

Z_mass_pre = mass_pre_NA(Z, window_size=window_size)
B_mass_pre = mass_pre_NA(B, window_size=window_size)
X_mass_pre = mass_pre_NA(X, window_size=window_size)

tictac <- Sys.time()
# maximum variance value of all subsequences

max_sigma2 = max(Z_mass_pre$sigmax2[B_mass_pre$data_mean == 1])

distanceProfile = matrix(0, MatrixProfileLen, 1) 
updatePos = matrix(FALSE, MatrixProfileLen, 1) # from rep_len vec to matr

#preallocate
QZ = rep.int(0, times = MatrixProfileLen)
QB = rep.int(0, times = MatrixProfileLen)
BZ = rep.int(0, times = MatrixProfileLen)
ZB = rep.int(0, times = MatrixProfileLen)
BX = rep.int(0, times = MatrixProfileLen)
XB = rep.int(0, times = MatrixProfileLen)

# what follows here switched to tsmp::dist_profile
# and still need to calc all stats as is done in fastfindNNPre and fastfindNN
QZ_tsmp_dp <- tsmp::dist_profile(data=Z, query=Z, window_size=window_size, method='v2')
QB_tsmp_dp <- tsmp::dist_profile(data=B, query=B, window_size=window_size, method='v2')
BZ_tsmp_dp <- tsmp::dist_profile(data=Z, query=B, window_size=window_size, method='v2')
ZB_tsmp_dp <- tsmp::dist_profile(data=B, query=Z, window_size=window_size, method='v2')
BX_tsmp_dp <- tsmp::dist_profile(data=X, query=B, window_size=window_size, method='v2')
XB_tsmp_dp <- tsmp::dist_profile(data=B, query=X, window_size=window_size, method='v2') 

QZ = QZ_tsmp_dp$last_product
QB = QB_tsmp_dp$last_product
BZ = BZ_tsmp_dp$last_product
ZB = ZB_tsmp_dp$last_product
BX = BX_tsmp_dp$last_product
XB = XB_tsmp_dp$last_product

#first 
firstQZ = QZ_tsmp_dp$last_product
firstQB = QB_tsmp_dp$last_product
firstBZ = BZ_tsmp_dp$last_product
firstZB = ZB_tsmp_dp$last_product
firstBX = BX_tsmp_dp$last_product
firstXB = XB_tsmp_dp$last_product

# generate inputs in the (function) environment for use by CalcLBDist

meani_or = ZB/QB
meanj_or = BZ/QB
sigmai_or2 = XB/QB-meani_or^2
sigmaj_or2 = BX/QB-meanj_or^2
q= (QZ/QB - (meani_or * meanj_or)) /sqrt(abs(sigmai_or2 * sigmaj_or2))

# account for m in CalcLB_dist
m = window_size

distanceProfile = CalcLB_dist(1, num, m, Vmax, Vmin, QZ, QB, BZ, ZB, BX, XB, meanz = Z_mass_pre$data_mean, sigmaz = Z_mass_pre$data_sd, meanb = B_mass_pre$data_mean, sigmab = B_mass_pre$data_sd,meanx = X_mass_pre$data_mean,sigmax = X_mass_pre$data_sd, sigmai_or2, sigmaj_or2, q, max_sigma2, D = distanceProfile) 

#apply exlcusion zone 
if (exclusion_zone > 0) {
exc_st = max(1, 1-exclusion_zone)
exc_ed = min(MatrixProfileLen, 1+exclusion_zone)
distanceProfile[exc_st:exc_ed,1] <- Inf
}
# evaluate initial matrix profile
MatrixProfile[1:MatrixProfileLen,1] = distanceProfile[1:MatrixProfileLen,1] 
MPindex[1:length(MatrixProfileLen)] = 1
updatePos = (distanceProfile < MatrixProfile) # updatePos should all be FALSE now
from = 2

# and now, the long weird dance
for (i in from:MatrixProfileLen) {

# update dot products
QZ[2:(num-m +1)] = QZ[1:(num -m)] -Z[(i-1)] * Z[1:(num -m)] + Z[(i +m -1)] * Z[(m +1) : num]
QB[2:(num-m +1)] = QB[1:(num -m)] -B[(i-1)] * B[1:(num -m)] + B[(i +m -1)] * B[(m +1) : num]
BZ[2:(num-m +1)] = BZ[1:(num -m)] -B[(i-1)] * Z[1:(num -m)] + B[(i +m -1)] * Z[(m +1) : num]
ZB[2:(num-m +1)] = ZB[1:(num -m)] -Z[(i-1)] * B[1:(num -m)] + Z[(i +m -1)] * B[(m +1) : num]
BX[2:(num-m +1)] = BX[1:(num -m)] -B[(i-1)] * X[1:(num -m)] + B[(i +m -1)] * X[(m +1) : num]
XB[2:(num-m +1)] = XB[1:(num -m)] -X[(i-1)] * B[1:(num -m)] + X[(i +m -1)] * B[(m +1) : num]

QZ[1] = firstQZ[i]
QB[1] = firstQB[i]
BZ[1] = firstBZ[i]
ZB[1] = firstZB[i]
BX[1] = firstBX[i]
XB[1] = firstXB[i]
# update stats for CalcLB_dist
meani_or = ZB/QB
meanj_or = BZ/QB
sigmai_or2 = XB/QB-meani_or^2
sigmaj_or2 = BX/QB-meanj_or^2
q= (QZ/QB - (meani_or * meanj_or)) /sqrt(abs(sigmai_or2 * sigmaj_or2))


distanceProfile = try(CalcLB_dist(i, num, m, Vmax, Vmin, QZ, QB, BZ, ZB, BX, XB, meanz = Z_mass_pre$data_mean, sigmaz = Z_mass_pre$data_sd, meanb = B_mass_pre$data_mean, sigmab = B_mass_pre$data_sd, meanx = X_mass_pre$data_mean, sigmax = X_mass_pre$data_sd, sigmai_or2, sigmaj_or2, q, max_sigma2, D = distanceProfile))
# update exclusion_zone
exc_st[i] = max(i, i - exclusion_zone)
exc_ed[i] = min(MatrixProfileLen, i + exclusion_zone)
distanceProfile[exc_st[i]:exc_ed[i], 1] = Inf

# figure out and store nearest neighbor
updatePos = (distanceProfile < MatrixProfile)

MatrixProfile[updatePos] = distanceProfile[updatePos]
#if (any(is.na(MatrixProfile) == TRUE)) browser()
MPindex[which(updatePos)] = i # index, as it says, but not intuitive
# for some debugging, prior to progress::
#if (window_size == 160 & i %% 100 == 0) print(i) 
#if (window_size == 2000 & i %% 1000 == 0) print(i)
#if (i %% 300 == 0) browser()
pb$tick()
}
tictac <- Sys.time() - tictac
message(sprintf('Finished in %.2f %s', tictac, units(tictac)))
MatrixProfile = sqrt(abs(MatrixProfile))

return(list(MatrixProfile = MatrixProfile, MPindex = MPindex, data = Z))

} 

