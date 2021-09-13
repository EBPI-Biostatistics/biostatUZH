#' Misc Tools of the Department of Biostatistics, EBPI, University of Zurich
#' 
#' The package \pkg{biostatUZH} provides functions developed at the
#' Department of Biostatistics at the Epidemiology, Biostatisics and Prevention
#' Institute, University of Zurich, Switzerland.  Covered topics include
#' confidence intervals (CIs) for proportions, CIs for operating characteristics of
#' diagnostic tests, Fagan-Nomogram, bootstrap CIs for the kappa coefficient,
#' intraclass-correlation coefficients including CIs, agreement for continuous
#' measurements (Bland-Altman plot), CI for the area under the curve (AUC),
#' Mantel-Haenszel estimator, Mantel-Cox hazard ratio estimator, CI for the
#' Kaplan-Meier estimate at given time points, CI for quantile of a
#' Kaplan-Meier or cumulative incidence estimate, sample size computations for
#' two-sample Mann-Whitney test, McNemar test, a binary diagnostic test,
#' natural re-parametrization of Weibull
#' output from survreg, hazard ratio and event time ratio interpretations, plot
#' to check the adequacy of the Weibull model.
#' 
#' Most relevant functions: \cr
#' \code{\link{agresti}} \cr                     
#' \code{\link{behrensTest}} \cr                 
#' \code{\link{blandAltman}} \cr                 
#' \code{\link{clopperPearson}} \cr              
#' \code{\link{combineSubgroups}} \cr            
#' \code{\link{computeICCrater}} \cr             
#' \code{\link{confIntAUC}} \cr                  
#' \code{\link{confIntAUCbinorm}} \cr            
#' \code{\link{confIntCorrelation}} \cr          
#' \code{\link{confIntDiagnostic}} \cr           
#' \code{\link{confIntICC}} \cr                  
#' \code{\link{confIntIndependentAUCDiff}} \cr   
#' \code{\link{confIntIndependentDiagnostic}} \cr
#' \code{\link{confIntIndependentProportion}} \cr
#' \code{\link{confIntKappa}} \cr                
#' \code{\link{confIntKM}} \cr                   
#' \code{\link{confIntMedian}} \cr               
#' \code{\link{confIntOddsRatio}} \cr            
#' \code{\link{confIntPairedAUCDiff}} \cr        
#' \code{\link{confIntPairedDiagnostic}} \cr     
#' \code{\link{confIntPairedProportion}} \cr     
#' \code{\link{confIntProportion}} \cr           
#' \code{\link{confIntRate}} \cr                 
#' \code{\link{confIntRateDiff}} \cr             
#' \code{\link{confIntRateWald}} \cr             
#' \code{\link{confIntRateWilson}} \cr           
#' \code{\link{confIntRiskDiff}} \cr             
#' \code{\link{confIntRiskRatio}} \cr            
#' \code{\link{faganPlot}} \cr                   
#' \code{\link{formatCI}} \cr                    
#' \code{\link{formatPercent}} \cr               
#' \code{\link{formatPval}} \cr                  
#' \code{\link{is.wholenumber}} \cr              
#' \code{\link{jeffreys}} \cr                    
#' \code{\link{mantelCoxHR}} \cr                 
#' \code{\link{mantelHaenszel}} \cr              
#' \code{\link{plotcorr}} \cr                    
#' \code{\link{populationSamplePlot}} \cr        
#' \code{\link{printWaldCI}} \cr                 
#' \code{\link{quantileCumInc}} \cr              
#' \code{\link{quantileKM}} \cr                  
#' \code{\link{sampleSizeMcNemar}} \cr           
#' \code{\link{sampleSizeSurvival}} \cr          
#' \code{\link{sampleSizeWilcoxTwoSample}} \cr   
#' \code{\link{standardErrorAUC}} \cr            
#' \code{\link{standardErrorAUCDiff}} \cr        
#' \code{\link{summaryROC}} \cr                  
#' \code{\link{survEvents}} \cr                  
#' \code{\link{survreg2weibull}} \cr             
#' \code{\link{survReport}} \cr                  
#' \code{\link{tableOR}} \cr                     
#' \code{\link{tableRegression}} \cr             
#' \code{\link{wald}} \cr                        
#' \code{\link{waldRate}} \cr                    
#' \code{\link{weibullDiag}} \cr                 
#' \code{\link{weibullReg}} \cr                  
#' \code{\link{wilson}} \cr                      
#' \code{\link{wilsonRate}} \cr
#'
#' Data sets: \cr
#' \code{\link{fischdaten}} \cr
#' \code{\link{larynx}} \cr
#' \code{\link{wiedat2b}} \cr                    
#'
#' @seealso \url{http://www.biostat.uzh.ch}
#'
#'
#'
#' 
#' @name biostatUZH-package
#' @aliases biostatUZH-package biostatUZH
#' @docType package
#' @keywords package
NULL





#' Data of the 'Fischstudie'
#' 
#' This dataset provides data of 308 patients enrolled in the 'Fischstudie',
#' described in Schumacher and Schulgen (2002).
#' 
#' 
#' @name fischdaten
#' @docType data
#' @format A data frame with 308 observations on the following 16 variables.
#' \describe{
#' \item{patnr}{Patient number.}
#' \item{geschl}{Sex: 0 = female, 1 = male.}
#' \item{alter}{Age in years.}
#' \item{groesse}{Height in cm.}
#' \item{gew0}{Weight at study entry, in kg.}
#' \item{gew15}{Weight at 15 days, in kg.}
#' \item{gew28}{Weight at 28 days, in kg.}
#' \item{chol0}{Cholesterin level at study entry, in mg/dl.}
#' \item{chol15}{Cholesterin level at 15 days, in mg/dl.}
#' \item{chol28}{Cholesterin level at 28 days, in mg/dl.}
#' \item{fisch}{Received fish (1) or not (0).}
#' \item{kalorie}{Received reduced food (1) or not (0).} }
#' @references Held, Leonhard, Rufibach, Kaspar, and Seifert, Burkhardt (2013).
#' Medizinische Statistik. Konzepte, Methoden, Anwendungen. Pearson Verlag,
#' Halbergmoos, Germany.
#' 
#' Schumacher, Martin and Schulgen, Gabi. (2002). Methodik klinischer Studien.
#' Springer, Berlin, Germany.
#' @source This dataset is used in Schumacher and Schulgen (2002) and Held et
#' al. (2013).
#' @keywords datasets
NULL





#' Survival Times of Larynx Cancer Patients
#' 
#' A study of 90 males with laryngeal cancer was performed, comparing survival
#' times. Each patient's age, year of diagnosis, and disease stage was noted.
#' Study published in Kardaun (1983), and available at the website for Klein
#' and Moeschberger (2003).
#' 
#' 
#' @name larynx
#' @docType data
#' @format A data frame with 90 observations on the following 5 variables.
#' \describe{
#' \item{stage}{Disease stage (1-4) from TNM cancer staging classification.}
#' \item{time}{Time from first treatment until death, or end of study.}
#' \item{age}{Age at diagnosis.}
#' \item{year}{Year of diagnosis.}
#' \item{death}{Indicator of death [1, if patient died at time t; 0, otherwise].} }
#' @references Kardaun, O. (1983). Statistical survival analysis of male
#' larynx-cancer patients-a case study.  \emph{Statistica Neerlandica},
#' \bold{37}, 103--125.
#' 
#' Klein, J. and Moeschberger, M. (2003). \emph{Survival analysis: techniques
#' for censored and truncated data}.  Springer.
#' @source
#' \url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/Publicfiles/DataFromSection/DataFromSectionTXT/Data_from_section_1.8.txt}
#' @keywords datasets
#' @examples
#' 
#' library(survival)
#' data(larynx)
#' Surv(larynx$time, larynx$death)
#' 
NULL





#' Data of two pancreatic Ca biomarkers
#' 
#' This dataset provides two paired diagnostic marker measurements for
#' diagnosis of pancreatic cancer.
#' 
#' 
#' @name wiedat2b
#' @docType data
#' @format A data.frame with 141 observations on the 3 variables:
#' \describe{
#' \item{y1:}{Biomarker 1.}
#' \item{y2:}{Biomarker 2.}
#' \item{d:}{Disease: 0 = non-cancer, 1 = cancer.} }
#' @references Wieand S, Gail MH, James BR, and James KL. \emph{A family of
#' nonparametric statistics for comparing diagnostic markers with paired or
#' unpaired data.} Biometrika 76(3):585-92. 1989.
#' 
#' Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for
#' classification and prediction}. Oxford University Press.
#' @keywords data
NULL



