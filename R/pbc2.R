#' pbc2 dataset
#'
#' pbc2 data from Mayo clinic
#'
#' @format Longitudinal dataset with 1945 rows and 19 columns for 312 patients
#' \describe{
#'   \item{id}{Patient identifier}
#'   \item{time}{Time measurement}
#'   \item{ascites}{Presence of ascites (Yes/No)}
#'   \item{hepatomegaly}{Presence of hepatomegaly (Yes/No)}
#'   \item{spiders}{Blood vessel malformations in the skin (Yes/No)}
#'   \item{edema}{Edema levels (No edema/edema no diuretics/edema despite diuretics)}
#'   \item{serBilir}{Level of serum bilirubin}
#'   \item{serChol}{Level of serum cholesterol}
#'   \item{albumin}{Level of albumin}
#'   \item{alkaline}{Level of alkaline phosphatase}
#'   \item{SGOT}{Level of aspartate aminotransferase}
#'   \item{platelets}{Platelet count}
#'   \item{prothrombin}{Prothrombin time}
#'   \item{histologic}{Histologic stage of disease}
#'   \item{drug}{Drug treatment (D-penicillmain/Placebo)}
#'   \item{age}{Age at enrollment}
#'   \item{sex}{Sex of patient}
#'   \item{years}{Time-to-event in years}
#'   \item{event}{Event indicator: 0 (alive), 1 (transplanted) and 2 (dead)}
#' }
#' @source pbc2 joineRML
#' @docType data
#' @name pbc2
#'
#' @examples
#' data(pbc2)
NULL
