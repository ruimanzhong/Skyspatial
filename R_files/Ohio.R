library(sf)
rm(list=ls())
dohio = read.csv("dataohiocomplete.csv")
map <- read_sf("ohioshp/fe_2007_39_county.shp")


d <- aggregate(
  x = dohio$y,
  by = list(county = dohio$NAME, year = dohio$year),
  FUN = sum
)
names(d) <- c("county", "year", "Y")
head(d)

dohio <- dohio[order(
  dohio$county,
  dohio$year,
  dohio$gender,
  dohio$race
), ]

dohio[1:20, ]

library(SpatialEpi)
n.strata <- 4
E <- expected(
  population = dohio$n,
  cases = dohio$y,
  n.strata = n.strata
)

nyears <- length(unique(dohio$year))
countiesE <- rep(unique(dohio$NAME),
                 each = nyears)

ncounties <- length(unique(dohio$NAME))
yearsE <- rep(unique(dohio$year),
              times = ncounties)

dE <- data.frame(county = countiesE, year = yearsE, E = E)

d <- merge(d, dE, by = c("county", "year"))

d$SIR <- d$Y / d$E

dw <- reshape(d,
              timevar = "year",
              idvar = "county",
              direction = "wide"
)

map <- merge(map, dw, by.x = "NAME", by.y = "county")
map_sf <- gather(map, year, SIR, paste0("SIR.", 1968:1988))
map_sf$year <- as.integer(substring(map_sf$year, 5, 8))












