## calibrate.R: censoring intercept c_int giving each target censoring rate for each
## censoring scale p.c (the SD-type spread of log C), under data_aft_int()'s model:
##   log X = 0.1 + 0.1 Z + 0.5 log Exp(1);  Y = 1 + (A - X) + Z + e
##   log C = c_int + 0.1 Y + 0.1 Z + p.c log Exp(1)
## Uses one large sample (n = 1e6, fixed seed). Output pasted into main.R (C_INT table).
P_C    <- c(0.75, 1, 1.5, 2, 3)
TARGET <- c(0.50, 0.65, 0.75, 0.85)
set.seed(20261113); n <- 1e6
Z <- rnorm(n); A <- rnorm(n, 2, 1); X <- exp(0.1 + 0.1 * Z + 0.5 * log(rweibull(n, 1, 1)))
y <- 1 + (A - X) + Z + rnorm(n); eg <- log(rweibull(n, 1, 1))
grid <- expand.grid(p.c = P_C, target = TARGET)
grid$c_int <- mapply(function(pc, tg) uniroot(function(ci) mean(X > exp(ci + 0.1 * y + 0.1 * Z + pc * eg)) - tg,
                                              c(-40, 40), tol = 1e-6)$root, grid$p.c, grid$target)
grid$c_int <- round(grid$c_int, 4)
print(grid)
cat("C_INT <- c(", paste(sprintf('"%s_%s" = %s', grid$p.c, grid$target, grid$c_int), collapse = ", "), ")\n")
