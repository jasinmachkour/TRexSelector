#' Ground truth for simulated GWAS
#'
#' A table containing the disease SNPs (disease allele, index, position, heterozygote risks,
#' homozygote risks, and rs-ID-number).
#'
#' @format A dataframe containing six columns:
#' \describe{
#'   \item{allele}{Disease allele (zero or one).}
#'   \item{idx}{Index of disease SNP.}
#'   \item{position}{Position of disease SNP.}
#'   \item{risk_hete}{Heterozygote risk.}
#'   \item{risk_homo}{Homozygote risk.}
#'   \item{rs}{rs-ID-number.}
#' }
"ground_truth"
