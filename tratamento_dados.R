library(readr)
library(dplyr)
library(lubridate)
library(stringr)


# --- Leitura ---

dados = read_csv(
  "dados_brutos.csv",
  locale = locale(grouping_mark = ".", decimal_mark = ",")
)

# Corrige nomes das colunas (read_csv altera caracteres especiais)
names(dados) = names(read.csv("matriz_transposta.csv"))


# --- Função auxiliar: converte texto de moeda para numérico ---

parse_moeda = function(x) {
  x |>
    str_remove_all("\\.") |>
    str_replace_all(",", ".") |>
    as.numeric()
}


# --- Transformação de tipos ---

dados = dados |>
  mutate(
    data_ref = as.Date(data_ref, "%d/%m/%Y"),
    across(where(is.character), parse_moeda)
  )


# --- Substituição de NAs por 0 (apenas em anos com ao menos um valor observado) ---

dados = dados |>
  mutate(ano = year(data_ref)) |>
  group_by(ano) |>
  mutate(
    across(
      where(is.numeric),
      ~ if (all(is.na(.x))) .x else ifelse(is.na(.x), 0, .x)
    )
  ) |>
  ungroup() |>
  select(-ano)


# --- Exportação ---

write.csv(dados, "dados_prefeitura_limpos_final.csv", row.names = FALSE)