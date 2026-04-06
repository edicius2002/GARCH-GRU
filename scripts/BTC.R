Sys.setenv(TZ='America/New_York')

library(dplyr)
library(rugarch)
library(TTR)
library(parallel)
library(zoo)

dane <- read.csv("../datasets/BTC.csv", sep = ",", dec=".", header = TRUE,
                 stringsAsFactors=FALSE)

dane$Date <- as.Date(dane$Date, format = "%m/%d/%Y")
dane <- dane[order(dane$Date), ]
dane<-dane[97:(nrow(dane)-187),]
dane <- na.omit(dane)

returns <- diff(log(dane$Close)) * 100
returns <- na.omit(returns)

dane <- dane[-1, ]
dane$Realized_Ret <- returns

ohlc <- dane %>% select(Open, High, Low, Close)
vGKYZ <- as.numeric(volatility(ohlc, calc = "gk.yz", N = 1))

a <- sqrt(mean(dane$Realized_Ret^2, na.rm = TRUE))
b <- sqrt(mean(vGKYZ^2, na.rm = TRUE))
SvOLHC <- (a/b) * vGKYZ

dane$SvOLHC <- SvOLHC
dane$SvOLHC_n10 <- zoo::rollmean(dane$SvOLHC, 10, fill = NA, align = "right")

n_steps <- 504
forecast.length <- length(dane$Realized_Ret) - n_steps
n.cores <- 4

ret <- as.numeric(dane$Realized_Ret)

cl <- makePSOCKcluster(n.cores)

clusterEvalQ(cl, {
  library(rugarch)
})

run_garch <- function(model, dist, data, cl) {
  spec <- ugarchspec(
    variance.model = list(model = model, garchOrder = c(1,1)),
    mean.model     = list(armaOrder = c(1,0), include.mean = TRUE),
    distribution.model = dist
  )
  
  roll_fit <- ugarchroll(
    spec,
    data = data,
    n.ahead = 1,
    forecast.length = forecast.length,
    refit.every = 1,
    refit.window = "moving",
    window.size = n_steps,
    calculate.VaR = FALSE,
    solver = "hybrid",
    cluster = cl
  )
  
  out <- data.frame(
    mu    = roll_fit@forecast$density$Mu,
    sigma = roll_fit@forecast$density$Sigma,
    skew  = roll_fit@forecast$density$Skew,
    shape = roll_fit@forecast$density$Shape
  )
  
  return(out)
}

# "sGARCH", "eGARCH", "gjrGARCH", "apARCH"
models <- c("sGARCH")
# "norm", "std", "sstd"
dists  <- c("norm")

results <- list()

times <- data.frame(
  model = character(),
  dist  = character(),
  time_sec = numeric(),
  stringsAsFactors = FALSE
)

for (m in models) {
  for (d in dists) {
    
    cat("Ejecutando", m, d, "...\n")
    
    t <- system.time({
      res <- run_garch(m, d, ret, cl)
    })
    
    colnames(res) <- paste0(m, "_", d, "_", c("mu","sigma","skew","shape"))
    
    results[[paste0(m,"_",d)]] <- res

    times <- rbind(times, data.frame(
      model = m,
      dist = d,
      time_sec = t["elapsed"]
    ))
    
    cat("Tiempo:", t["elapsed"], "segundos\n\n")
  }
}

stopCluster(cl)

if (length(results) == 0) {
  stop("No hay resultados generados")
}

nn <- nrow(results[[1]])
dane_out <- dane %>% tail(nn)

for (nm in names(results)) {
  dane_out <- cbind(dane_out, results[[nm]])
}

# "Date", "Realized_Ret", "SvOLHC", "SvOLHC_n10"
base_cols <- c("Date", "Realized_Ret", "SvOLHC", "SvOLHC_n10")

selected_cols <- c()

for (m in models) {
  for (d in dists) {
    selected_cols <- c(selected_cols,
                       paste0(m, "_", d, "_mu"),
                       paste0(m, "_", d, "_sigma"))
    
    if (d != "norm") {
      selected_cols <- c(selected_cols,
                         paste0(m, "_", d, "_skew"),
                         paste0(m, "_", d, "_shape"))
    }
  }
}

dane_out <- dane_out %>%
  select(any_of(base_cols), any_of(selected_cols)) %>%
  arrange(Date)

write.csv(dane_out, "../results/BTC_test.csv", row.names = FALSE)

cat("Archivo exportado con", nrow(dane_out), "filas y", ncol(dane_out), "columnas.\n")