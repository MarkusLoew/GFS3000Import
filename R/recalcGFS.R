#' Recalculate gas exchange parameters measured by a WALZ GFS3000 IRGA
#' 
#' Recalculate leaf gas exchange, assimilation, evaporation, gH2O, VPD and ci after a change of leaf area. Based on equations from the WALZ GFS3000 manual. 
#'
#' @param data Data frame with the gas exchange data
#' @param Area When Area = NA uses the vector Area in the data. If value is provided for Area it uses this value.
#' @param inplace Logical. When FALSE, only returns a data frame with the parameters that changed (default). When TRUE, returns the full data with the recalculated data in place.
#' @return data.frame with the recalculated values.
#' @export


WalzGFS3000recalc <- function(data, Area = NA, inplace = FALSE) {

# to re-calculate:
# Assimilation, Evaporation, gH2O, VPD and ci

# use area from data frame if no new area is provided as argument
if (is.na(Area)) {
  Area <- data$Area}
  
# Calc CO2 mole fraction in cuvette (ca)
# ca = CO2sam - dCO2ZP

# Calc differential CO2 mole fraction in measure mode (dCO2MP)
# dCO2MP = CO2sam - CO2abs

# Calc H2O mole fraction in cuvette (wa)
# wa = H2Osam - dH2OZP

# Calc of differential H2O mole fraction in MP mode (dH2OMp)
# dH2OMp = H2Osamp - H2Oabs

# Calc SVP(Tcuv) saturation vapour pressure at Tcuv according to Goff-Gratch
# Calc SVP(Tleaf) saturation vapour pressure at Tcuv according to Goff-Gratch

# SVPt: saturation vapour pressure over a plane surface of pure ordinary liquid water dependent on T [hPa]
# T: absolute thermodynamic temperature [K]
# Ts: steam point temperature (373.16 K)
# SPws: saturation pressure of pure ordinary liquid water at stem point temperate (1013.246 hPa, i.e. 1 standard atmosphere)

MySVPt <- function(t) {
# following Goff / Gratch equation from Walz GFS3000 manual,
# tested a few approaches
# t is temperature in degC

# SVPt in [hPa]
  Ts <- 373.16
  SPws <- 1013.246
  t <- t + 273.16
  T <- t
  logSVPt <- -7.90298 * ((Ts/t) - 1) + 5.02808 * log10(Ts/t) -
             (1.3816 * 10^-7) * (10^(11.344 *(1-t/Ts)) - 1) +
             (8.1328 * 10^-3) * (10^((-3.49149 * (Ts/t-1))-1)) +
             log10(SPws)
   # other approaches
	#  logSVPt <-  -7.90298 * (373.16/T-1) +
	#                     5.02808 * log(373.16/T) -
	#                     1.3816 * 10-7 * (10^(11.344 * (1-T/373.16))  -1) +
	#                     8.1328 * 10-3 * (10^(-3.49149 * (373.16/T-1))  -1) +
	#                    log(1013.246) 
	#  A <- -7.90298 * (373.16/T - 1)
	#  B <- 5.02808 * log10(373.16/T)
	#  C <- 1.3816*10^-7 * 10^(11.344 * (1- (T/373.16)) -1)
	#  D <- 8.1328 * 10^-3 * (10^-3.49149 * (373.16/T - 1) -1)
	#  E <- log10(1013.246)
	#  logSVPt <- A + B - C + D + E

  return(10^logSVPt)
}


# relative humidity
# rh = (wa * Pamb) / (SVP(Tcuv))

# calculation of transpiration Rate (E)
# E = (ue * (w0 - we)) / (LA * (1 - w0))
# in WALZ: ue = Flow
#         w0 - we = dH2OMP - dH2OZP
#         w0 = wa
#         LA = area
# therefore
#E <- (Flow * (dH2OMP - dH2OZP)) / Area * (1 - wa)
E <- (data$Flow * ((data$dH2OMP * 0.000001) - (data$dH2OZP * 0.000001))) / ((Area * 0.1) * (1 - (data$wa * 0.000001)))
E.div <- E / 1000
#E <- E * 10

# calculation of VPD
# where SVP(Tleaf) is saturation vapour pressure at Tleaf following Goff / Gratch
# and, following Walz GFS3000 manual: Pcuv = Pamb
# VPD = ( (SVP(Tleaf) / Pamb) - wa ) / (1 - ( (SVP(Tleaf)/Pamb) + wa) / 2 )
# VPD to be double-checked, regarding potential overestimation.
VPD <- ( (MySVPt(data$Tleaf) / data$Pamb) - data$wa ) / 
       (1 - (( (MySVPt(data$Tleaf) / 
       data$Pamb) + data$wa) / 2 )) *  10

# calculation of Water Vapour conductance (GH2O)
GH2O <- (E / data$VPD) * 1000
#GH2O <- GH2O /10 # mmol m-2

# calculation of assimilation rate
A <- (((data$Flow * ((data$dCO2ZP * 0.000001) - (data$dCO2MP * 0.000001))) / (Area * 0.1)) - E * (data$ca * 0.000001)) * 1000

# calculation of intercellular CO2 mole fraction (ci)
#ci = ((gCO2 - (E/2)) * ca - A) / (gCO2 + (E/2))
# relationship between conductance for CO2 to conductance of H2O:
# gCO2 = GH2O / 1.56 
# with on-the-fly conversions for units and gCO2
# ci <- (((GCO2 / 1.56) - (E.div / 2)) * data$ca - A) / (GH2O / 1.56) + (E.div / 2)
ci <- ((((GH2O / 1.56) - ((E) / 2))) * (data$ca * 0.0000001) - (A/10000)) / ((GH2O / 1.56) + ((E)/2)) * 10000000


if (inplace == FALSE) {
 # in case of separate result table is asked for
 out <- data.frame(Area = Area,
                   E = E,
                   GH2O = GH2O,
                   A = A,
                   ci = ci,
                   VPD = VPD)
 return(out)
 } else {
  data$Area 	<- Area
  data$E 	<- E
  data$GH2O 	<- GH2O
  data$A 	<- A
  data$ci 	<- ci
  data$VPD 	<- VPD
  return(data)
 }
}

