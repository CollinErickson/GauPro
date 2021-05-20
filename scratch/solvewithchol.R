if (F) {
  backsolve(kchol, backsolve(kchol, Z, transpose = T))
}
#
# If kchol <- chol(k),
# then solvewithchol(kchol, Z) is same as solve(k, Z)
solvewithchol <- function(kchol, Z) {
  backsolve(kchol, backsolve(kchol, Z, transpose = T))
}

calc_inv_from_chol <- FALSE
