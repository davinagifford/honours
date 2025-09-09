



# Packages
library(ncdf4)
library(dplyr)
library(lubridate)

# ---- file ----
fn <- "D:/HONOURS/DavinaG_2025_Honours/data/EAC_filled-daily-distance-depth-gridded-product_20120401-20220727.nc"

nc <- nc_open(fn)

# Helper to pick a variable name from candidates
pick <- function(cands) {
  hit <- cands[cands %in% names(nc$var)]
  if (!length(hit)) stop(sprintf("None of %s found in file", paste(cands, collapse=", ")))
  hit[1]
}
pick_dim <- function(cands) {
  hit <- cands[cands %in% names(nc$dim)]
  if (!length(hit)) stop(sprintf("None of %s found in file dims", paste(cands, collapse=", ")))
  hit[1]
}

# Variables / dims (common names in this CSIRO product)
vvar <- pick(c("VCUR"))
xnm  <- pick_dim(c("DISTANCE","distance","x"))
znm  <- pick_dim(c("DEPTH","depth","z"))
tnm  <- pick_dim(c("TIME","time"))

V   <- ncvar_get(nc, vvar)           # expect [distance, depth, time]
x   <- nc$dim[[xnm]]$vals            # km in this product
z   <- nc$dim[[znm]]$vals            # m (positive down)
t   <- nc$dim[[tnm]]$vals
tu  <- nc$dim[[tnm]]$units
nc_close(nc)

# ---- time parsing (CF units: "days since YYYY-MM-DD ...", etc.) ----
origin_str <- sub(".*since\\s+", "", tu)
origin_str <- sub(" ?UTC.*$", "", origin_str)
origin <- ymd_hms(origin_str, quiet = TRUE)
if (is.na(origin)) origin <- ymd(origin_str, quiet = TRUE)

sec_per <- if (grepl("^days", tu)) 86400 else if (grepl("^hours", tu)) 3600 else 1
time <- origin + t * sec_per

# ---- ensure array is [x, z, t] ----
dims <- dim(V)
if (length(dims) != 3) stop("VCUR is not 3D as expected")
ix <- which(dims == length(x))[1]
iz <- which(dims == length(z))[1]
it <- setdiff(1:3, c(ix, iz))[1]
V <- aperm(V, perm = c(ix, iz, it))

# ---- cell sizes and integration domain ----
# Convert distance from km to m
x_m <- x * 1000
# centred thickness from 1D coords
cell_w <- function(v){
  n <- length(v); if(n<2) stop("Need >=2 coords")
  d <- diff(v)
  w <- c(d[1], (d[-1] + d[-length(d)])/2, d[length(d)])
  if (length(w) != n) w <- c(d[1], rep(d, length.out=n-1))
  w
}
dx <- cell_w(x_m)
dz <- cell_w(z)

# Optional: limit depth to 0–1500 m (matches many EAC transport defs)
z_mask <- (z >= 0 & z <= 1500)
x_mask <- rep(TRUE, length(x))

dxu <- dx[x_mask]
dzu <- dz[z_mask]
Vu  <- V[x_mask, z_mask, , drop = FALSE]

# Cell areas [x,z]
A <- outer(dxu, dzu)  # m^2

# ---- integrate to transport ----
# In file: VCUR is northward positive. For EAC poleward (southward) positive, flip sign.
flip_sign <- TRUE          # southward -> positive
poleward_only <- FALSE     # set TRUE to exclude equatorward (northward) flow

int_one <- function(k){
  Vk <- Vu[,,k]
  if (flip_sign) Vk <- -Vk
  if (poleward_only) Vk[Vk < 0] <- 0
  sum(Vk * A, na.rm = TRUE)  # m^3/s
}

T_m3s <- vapply(seq_len(dim(Vu)[3]), int_one, numeric(1))
T_Sv  <- T_m3s / 1e6
daily <- tibble(time = time, transport_Sv = T_Sv) |>
  arrange(time)

monthly <- daily |>
  mutate(month = floor_date(time, "month")) |>
  summarise(transport_Sv = mean(transport_Sv, na.rm = TRUE), .by = month) |>
  arrange(month)

# Peek
print(head(daily))
print(head(monthly))
