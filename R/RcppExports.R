# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

wls_exp <- function(alm0, almc, alpha, m, no, ni, x, r, xv, v, intr, ju, vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr) {
    .Call(`_multiview_wls_exp`, alm0, almc, alpha, m, no, ni, x, r, xv, v, intr, ju, vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr)
}

spwls_exp <- function(alm0, almc, alpha, m, no, ni, x, xm, xs, r, xv, v, intr, ju, vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr) {
    .Call(`_multiview_spwls_exp`, alm0, almc, alpha, m, no, ni, x, xm, xs, r, xv, v, intr, ju, vp, cl, nx, thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr)
}

