#' This is the code for the paper "Admissible Motif Discovery with Missing
#' Data" by Yan Zhu, Abdullah Mueen and Eamonn Keogh.
#'
#' Translated from Matlab to R by Christoper English
#' All mistakes are mine as their code worked and results seem reliable
#' as executed in Octave
#' this function in intended to be called from within StompSelfJoinNA

#' @param i 'int' iterator of for loop
#' @param num 'int' length of time_series
#' @param m 'int' window_size
#' @param Vmax 'vector' at 'i' max of 'm' forward values - calculated once in StompSelfJoinNA
#' @param Vmin 'vector' at 'i' min of 'm' forward values - calculated once in StompSelfJoinNA
#' @param QZ 'vector' last_product as calc'd by tsmp::dist_profile
#' @param QB 'vector' last_product as calc'd by tsmp::dist_profile
#' @param BZ 'vector' last_product as calc'd by tsmp::dist_profile
#' @param ZB 'vector' last_product as calc'd by tsmp::dist_profile
#' @param BX 'vector' last_product as calc'd by tsmp::dist_profile
#' @param XB 'vector' last_product as calc'd by tsmp::dist_profile
#' @param meanz 'vector' as calc'd by mass_pre_NA for Z - See StompSelfJoinNA
#' @param sigmaz 'vector' as calc'd by mass_pre_NA for Z - See StompSelfJoinNA
#' @param meanb 'vector' as calc'd by mass_pre_NA for B - See StompSelfJoinNA
#' @param sigmab 'vector' as calc'd by mass_pre_NA for B - See StompSelfJoinNA
#' @param meanx 'vector' as calc'd by mass_pre_NA for X - See StompSelfJoinNA
#' @param sigmax 'vector' as calc'd by mass_pre_NA for X - See StompSelfJoinNA
#' @param sigmai_or2
#' @param sigmaj_or2
#' @param q
#' @param max_sigma2
#' @param D
#' @return Returns 'D' the distanceProfile
CalcLB_dist <- function(i,num,m,Vmax,Vmin,QZ,QB,BZ,ZB,BX,XB,meanz,sigmaz,meanb,sigmab,meanx,sigmax,sigmai_or2,sigmaj_or2,q,max_sigma2, D) {

for (j in 1:(num-m+1)){
     
if (identical(QB[j], m)) { # to identical() from == and back to == and back to identical

D[j] = 2*m*(1-q[j])       #QB[j] = m
#print('case_1')
} else {                    # QB[j] != m
if (meanb[i] > meanb[j]) {

sigma_or2 = sigmaj_or2[j]
sigma2 = sigmaz[i]^2
l_r = meanb[i]*m
mean_r = meanz[i]/meanb[i]
sigma_r2 = (sigmaz[i]^2 + meanz[i]^2)/meanb[i] - mean_r^2
} else {                   #meanb[i] !> meanb[j]
sigma_or2 = sigmaj_or2[j]
sigma2 = sigmaz[j]^2
l_r = meanb[j]*m
mean_r = meanz[j]/meanb[j]
sigma_r2 = (sigmaz[j]^2 + meanz[j]^2)/meanb[j] - mean_r^2
}
if (max(meanb[i], meanb[j]) == 1) { #to ==

if (meanb[i] > meanb[j]) { #Ti, m doesn't have missing value at i
sigma_or2 = sigmai_or2[j]
sigma2 = sigmaz[j]^2
}

if (q[j] <= 0) { # to identical and back to <=

D[j] = QB[j] * sigma_or2/sigma2
#print('case_2')
} else {
D[j] = QB[j] * sigma_or2/sigma2*(1-q[j]^2)
#print('case_3')
}
} else {

meanr_i = meanz[i]/meanb[i]
meanr_j = meanz[j]/meanb[j]
maxsigma2_i=meanb[i]*((sigmaz[i]^2+meanz[i]^2)/meanb[i]-meanr_i^2+Vmax[i]*Vmin[i]+meanr_i*(meanr_i-Vmax[i]-Vmin[i]))+(Vmax[i]-Vmin[i])^2/4
maxsigma2_j=meanb[j]*((sigmaz[j]^2+meanz[j]^2)/meanb[j]-meanr_j^2+Vmax[j]*Vmin[j]+meanr_j*(meanr_j-Vmax[j]-Vmin[j]))+(Vmax[j]-Vmin[j])^2/4
v=max(sigmai_or2[j]/maxsigma2_i,sigmaj_or2[j]/maxsigma2_j)
if (q[j] <= 0) {
D[j] = QB[j]*v 
#print('case_4')
} else {
D[j] = QB[j]*v*(1-q[j]^2)
#print('case_5')
}
#if (i %% 599 == 0 & j %% 599 == 0) browser()
}
}
}
return(D)
}

