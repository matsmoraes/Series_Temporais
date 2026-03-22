# =============================================================================
# Custo dos Servidores Contratados — Modelagem de Séries Temporais
# Matheus de Moraes Neves | Eduardo Henrique de Souza Dias
# =============================================================================


# --- Setup: pacotes ---

library(ggplot2)
library(forecast)
library(tseries)
library(lubridate)
library(tidyverse)
library(deflateBR)
library(urca)
library(knitr)


# --- Leitura dos dados ---

dados = read_csv("dados_prefeitura_limpos_final.csv")
dados$data_ref = as.Date(dados$data_ref)


# --- Deflacionamento (IPCA) e gráfico nominal vs real ---

dados = dados %>%
  mutate(
    Custo_Real = deflate(
      nominal_values = Geral.Custo.dos.Servidores.Contratados,
      nominal_dates  = data_ref,
      real_date      = "10/2025",
      index          = "ipca"
    )
  )

ggplot(dados) +
  geom_line(aes(x = data_ref,
                y = Geral.Custo.dos.Servidores.Contratados,
                color = "Nominal"),
            linetype = "dashed", linewidth = 0.7) +
  geom_line(aes(x = data_ref, y = Custo_Real, color = "Real"), linewidth = 1) +
  scale_color_manual(values = c("Nominal" = "gray60", "Real" = "steelblue")) +
  scale_y_continuous(labels = scales::label_number(prefix = "R$ ",
                                                   big.mark = ".",
                                                   decimal.mark = ",")) +
  theme_minimal() +
  labs(x = "Ano", y = "Valor", color = "Série")


# --- Série temporal e gráfico de evolução ---

serie = ts(dados$Custo_Real,
           start     = c(2021, 1),
           frequency = 12)

ggplot(dados, aes(x = data_ref, y = Custo_Real)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(aes(color = format(data_ref, "%m") == "01"), size = 2) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  scale_y_continuous(labels = scales::label_number(prefix = "R$ ",
                                                   big.mark = ".",
                                                   decimal.mark = ",")) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  theme_minimal() +
  labs(x = "Tempo", y = "Custo (R$)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- Decomposição STL ---

decomposicao = stl(serie, s.window = "periodic")

autoplot(decomposicao) +
  theme_minimal() +
  labs(x = "Ano")


# --- Boxplot mensal ---

ggplot(dados, aes(format(data_ref, "%m"), Custo_Real)) +
  geom_boxplot(fill = "steelblue", alpha = 0.4) +
  scale_y_continuous(labels = scales::label_number(prefix = "R$ ",
                                                   big.mark = ".",
                                                   decimal.mark = ",")) +
  theme_minimal() +
  labs(x = "Mês (01 = Janeiro)", y = "Distribuição de Custos")


# --- Transformação Box-Cox (lambda) ---

lambda_valor = BoxCox.lambda(serie)
cat("Valor de Lambda estimado:", round(lambda_valor, 4))


# --- Teste KPSS e número de diferenciações ---

serie_transf = BoxCox(serie, lambda = lambda_valor)

kpss_result = ur.kpss(serie_transf)
summary(kpss_result)


# --- Número de diferenciações sugeridas ---

ndiffs_val  = ndiffs(serie_transf)
nsdiffs_val = nsdiffs(serie_transf)
cat("Diferenciações simples (d) sugeridas:  ", ndiffs_val,  "\n")
cat("Diferenciações sazonais (D) sugeridas:", nsdiffs_val, "\n")


# --- ACF/PACF: série transformada sem diferenciação ---

ggtsdisplay(serie_transf,
            main    = "Série Transformada — Sem Diferenciação",
            lag.max = 36)


# --- ACF/PACF: série diferenciada (estacionária) ---

serie_estac = diff(diff(serie_transf, lag = 12), differences = 1)

ggtsdisplay(serie_estac,
            main    = "Série Diferenciada — Estacionária",
            lag.max = 36)


# --- auto.arima (busca exaustiva) ---

modelo_auto = auto.arima(serie,
                         lambda        = lambda_valor,
                         stepwise      = FALSE,
                         approximation = FALSE)
summary(modelo_auto)


# --- Modelos candidatos e tabela AICc/BIC ---

m1 = Arima(serie, order = c(0,1,1), seasonal = c(1,0,0), lambda = lambda_valor)
m2 = Arima(serie, order = c(1,1,1), seasonal = c(1,0,0), lambda = lambda_valor)
m3 = Arima(serie, order = c(0,1,1), seasonal = c(0,1,1), lambda = lambda_valor)
m4 = Arima(serie, order = c(1,1,0), seasonal = c(1,0,0), lambda = lambda_valor)
m5 = Arima(serie, order = c(0,1,2), seasonal = c(1,0,0), lambda = lambda_valor)

tab_modelos = data.frame(
  Modelo = c("SARIMA(0,1,1)(1,0,0)[12]",
             "SARIMA(1,1,1)(1,0,0)[12]",
             "SARIMA(0,1,1)(0,1,1)[12]  *",
             "SARIMA(1,1,0)(1,0,0)[12]",
             "SARIMA(0,1,2)(1,0,0)[12]"),
  AICc   = c(m1$aicc, m2$aicc, m3$aicc, m4$aicc, m5$aicc),
  BIC    = c(m1$bic,  m2$bic,  m3$bic,  m4$bic,  m5$bic)
)

tab_modelos = tab_modelos %>%
  mutate(across(c(AICc, BIC), \(x) round(x, 2)))

kable(tab_modelos, align = "lcc")


# --- Modelo final selecionado ---

modelos_lista = list(m1, m2, m3, m4, m5)
aicc_vals     = sapply(modelos_lista, function(m) m$aicc)
modelo_final  = modelos_lista[[which.min(aicc_vals)]]

summary(modelo_final)


# --- Diagnóstico dos resíduos (checkresiduals) ---

checkresiduals(modelo_final)


# --- Ljung-Box ---

lb_test = Box.test(residuals(modelo_final), lag = 12, type = "Ljung-Box", fitdf = 2)
cat("Ljung-Box: Q* =", round(lb_test$statistic, 4),
    "| p-valor =", round(lb_test$p.value, 4), "\n")


# --- Shapiro-Wilk ---

sw_test = shapiro.test(residuals(modelo_final))
cat("Shapiro-Wilk: W =", round(sw_test$statistic, 5),
    "| p-valor =", round(sw_test$p.value, 5), "\n")


# --- Divisão treino/teste: SARIMA, ETS e Prophet ---

library(prophet)

h_teste  = 12
n_total  = length(serie)
n_treino = n_total - h_teste

treino_ts = subset(serie, end   = n_treino)
teste_ts  = subset(serie, start = n_treino + 1)

# --- SARIMA ---
modelo_sarima_treino = Arima(treino_ts,
                             order    = c(0, 1, 1),
                             seasonal = c(0, 1, 1),
                             lambda   = lambda_valor)
prev_sarima = forecast(modelo_sarima_treino, h = h_teste)

# --- ETS ---
modelo_ets_treino = ets(treino_ts)
prev_ets = forecast(modelo_ets_treino, h = h_teste)

# --- Prophet ---
df_prophet_treino = data.frame(
  ds = dados$data_ref[1:n_treino],
  y  = dados$Custo_Real[1:n_treino]
)
df_prophet_teste = data.frame(
  ds = dados$data_ref[(n_treino + 1):n_total]
)

modelo_prophet_treino = prophet(
  yearly.seasonality = TRUE,
  weekly.seasonality = FALSE,
  daily.seasonality  = FALSE
)
modelo_prophet_treino = add_country_holidays(modelo_prophet_treino,
                                             country_name = "BR")
modelo_prophet_treino = fit.prophet(modelo_prophet_treino, df_prophet_treino)

prev_prophet_teste = predict(modelo_prophet_treino, df_prophet_teste)


# --- Tabela de comparação RMSE/MAPE ---

rmse_fn = function(real, prev) sqrt(mean((real - prev)^2))
mape_fn = function(real, prev) mean(abs((real - prev) / real)) * 100

y_real    = as.numeric(teste_ts)
y_sarima  = as.numeric(prev_sarima$mean)
y_ets     = as.numeric(prev_ets$mean)
y_prophet = prev_prophet_teste$yhat

tab_comp = data.frame(
  Modelo = c("SARIMA(0,1,1)(0,1,1)[12]",
             "ETS (seleção automática)",
             "Prophet"),
  RMSE   = c(rmse_fn(y_real, y_sarima),
             rmse_fn(y_real, y_ets),
             rmse_fn(y_real, y_prophet)),
  MAPE   = c(mape_fn(y_real, y_sarima),
             mape_fn(y_real, y_ets),
             mape_fn(y_real, y_prophet))
)

tab_comp = tab_comp %>%
  mutate(
    RMSE = scales::number(RMSE, big.mark = ".", decimal.mark = ",", accuracy = 1),
    MAPE = paste0(round(MAPE, 2), "%")
  )

kable(tab_comp, align = "lcc", col.names = c("Modelo", "RMSE (R$)", "MAPE"))


# --- Gráfico comparação modelos no teste ---

datas_teste = dados$data_ref[(n_treino + 1):n_total]

df_comp = data.frame(
  Data    = datas_teste,
  Real    = y_real,
  SARIMA  = y_sarima,
  ETS     = y_ets,
  Prophet = y_prophet
) %>%
  pivot_longer(cols      = c(Real, SARIMA, ETS, Prophet),
               names_to  = "Serie",
               values_to = "Valor")

ggplot(df_comp, aes(x = Data, y = Valor, color = Serie, linetype = Serie)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("Real"    = "black",
                                "SARIMA"  = "steelblue",
                                "ETS"     = "tomato",
                                "Prophet" = "green4")) +
  scale_linetype_manual(values = c("Real"    = "solid",
                                   "SARIMA"  = "dashed",
                                   "ETS"     = "dotdash",
                                   "Prophet" = "dotted")) +
  scale_y_continuous(labels = scales::label_number(prefix = "R$ ",
                                                   big.mark = ".",
                                                   decimal.mark = ",")) +
  theme_minimal() +
  labs(x = "Mês", y = "Custo Real (R$)", color = NULL, linetype = NULL)


# --- Tabela de coeficientes estimados ---

coef_df = data.frame(
  Coeficiente   = c("ma1  ($\\hat{\\theta}_1$)", "sma1 ($\\hat{\\Theta}_1$)"),
  Estimativa    = c(-0.5885,  0.6412),
  `Erro Padrão` = c( 0.1142,  0.3495),
  check.names   = FALSE
)

kable(coef_df, align = "lcc", row.names = FALSE, escape = FALSE)


# --- Previsão 12 meses: gráfico ---

previsao = forecast(modelo_final, h = 12)

autoplot(previsao) +
  ggtitle("Previsão do Custo Real — Próximos 12 Meses") +
  xlab("Ano") + ylab("Custo Real (R$)") +
  theme_minimal() +
  scale_y_continuous(
    labels = scales::label_number(prefix = "R$ ",
                                  big.mark = ".",
                                  decimal.mark = ",")
  )


# --- Previsão 12 meses: tabela ---

datas_prev = seq(
  from       = max(dados$data_ref) %m+% months(1),
  by         = "month",
  length.out = 12
)

tab_prev = data.frame(
  `Período`  = format(datas_prev, "%b/%Y"),
  `Previsão` = scales::number(as.numeric(previsao$mean),
                              prefix = "R$ ", big.mark = ".",
                              decimal.mark = ",", accuracy = 1),
  `LI 80%`   = scales::number(as.numeric(previsao$lower[, 1]),
                              prefix = "R$ ", big.mark = ".",
                              decimal.mark = ",", accuracy = 1),
  `LS 80%`   = scales::number(as.numeric(previsao$upper[, 1]),
                              prefix = "R$ ", big.mark = ".",
                              decimal.mark = ",", accuracy = 1),
  `LI 95%`   = scales::number(as.numeric(previsao$lower[, 2]),
                              prefix = "R$ ", big.mark = ".",
                              decimal.mark = ",", accuracy = 1),
  `LS 95%`   = scales::number(as.numeric(previsao$upper[, 2]),
                              prefix = "R$ ", big.mark = ".",
                              decimal.mark = ",", accuracy = 1),
  check.names = FALSE
)

kable(tab_prev, align = "lrrrrr")

