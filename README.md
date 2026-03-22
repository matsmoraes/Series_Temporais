# Custo dos Servidores Contratados: Modelagem de Séries Temporais

Trabalho de modelagem e previsão da série temporal referente ao custo mensal com servidores contratados temporariamente pela Prefeitura Municipal de Uberlândia (MG).

**Autores:** Matheus de Moraes Neves e Eduardo Henrique de Souza Dias

## Sobre o projeto

Os dados foram extraídos do [Portal da Transparência de Uberlândia](https://www.uberlandia.mg.gov.br/portal-da-transparencia/dados-abertos/catalogo-de-dados-abertos/) e cobrem o período de janeiro de 2021 a outubro de 2025, totalizando 58 observações mensais. A variável analisada é o custo bruto com servidores contratados, deflacionada pelo IPCA com base em outubro de 2025, de modo a isolar o crescimento real dos gastos dos efeitos nominais da inflação.

A série apresenta três características marcantes: uma tendência de crescimento expressivo entre 2021 e 2023 com posterior estabilização, uma sazonalidade anual muito forte (com quedas drásticas em janeiro, mês de encerramento e renovação de contratos temporários) e heterocedasticidade, com variância crescendo junto com o nível da série.

## Tratamento dos dados

O dado bruto vem do Portal da Transparência no formato de relatório: cada aba corresponde a um ano e setor (Geral ou Educação), e dentro de cada aba as variáveis ficam nas linhas e os meses nas colunas, o chamado formato transposto. Chegar no CSV final exigiu duas etapas distintas.

A primeira foi manual, feita no Excel. Como o relatório cresceu ao longo dos anos, cada aba tinha um número diferente de indicadores. Foi necessário criar uma estrutura unificada com todos os indicadores presentes em qualquer ano e reaplicá-la em todas as abas para garantir que tivessem exatamente as mesmas linhas na mesma ordem. Com isso padronizado, as abas foram consolidadas em um único dataframe e transpostas para o formato analítico padrão (uma linha por mês, uma coluna por variável), gerando a `matriz_transposta.csv`.

A segunda etapa foi o script `tratamento_dados.R`, que recebe essa matriz e faz o tratamento de tipos: converte a coluna de data, parseia as colunas monetárias que chegam como texto com formatação brasileira (ponto como separador de milhar, vírgula como decimal) e substitui NAs por zero, respeitando o critério de só fazer essa substituição em anos que tenham ao menos uma observação válida na coluna. O resultado é exportado como `dados_prefeitura_limpos_final.csv`.

## Metodologia

O processo seguiu as etapas clássicas de Box-Jenkins:

1. **Análise exploratória:** decomposição STL, gráficos de evolução temporal e boxplots mensais para caracterizar tendência, sazonalidade e heterocedasticidade.
2. **Estabilização da variância:** transformação Box-Cox com λ estimado em 0,2997.
3. **Verificação de estacionariedade:** teste KPSS, que indicou a necessidade de uma diferenciação simples (d = 1) e uma diferenciação sazonal (D = 1).
4. **Identificação dos parâmetros:** gráficos ACF e PACF da série estacionária, apontando para estrutura de médias móveis com q = 1 e Q = 1.
5. **Seleção do modelo:** busca exaustiva via `auto.arima()` complementada pela estimação manual de modelos candidatos, comparados por AICc e BIC.
6. **Diagnóstico dos resíduos:** testes de Ljung-Box e Shapiro-Wilk, além de análise gráfica.
7. **Validação preditiva:** hold-out de 12 meses (nov/2024 a out/2025), com comparação entre SARIMA, ETS e Prophet por RMSE e MAPE.

## Modelo selecionado

**SARIMA(0,1,1)(0,1,1)[12]** com transformação Box-Cox (λ = 0,2997).

| Coeficiente | Estimativa | Erro Padrão |
|---|---|---|
| ma1 (θ₁) | -0,5885 | 0,1142 |
| sma1 (Θ₁) | 0,6412 | 0,3495 |

O modelo foi selecionado com AICc de 428,05, mais de 120 pontos abaixo dos modelos alternativos sem diferenciação sazonal. O teste de Ljung-Box (Q* = 6,23, p-valor = 0,796) confirmou ausência de autocorrelação serial nos resíduos. A não-normalidade dos resíduos (Shapiro-Wilk, p ≈ 0), reflexo das quedas abruptas em janeiro, implica cautela na interpretação dos intervalos de predição, sem comprometer a validade das estimativas pontuais.

## Comparação de modelos (conjunto de teste)

| Modelo | RMSE | MAPE |
|---|---|---|
| SARIMA(0,1,1)(0,1,1)[12] | R$ 1.462.587 | 77,04% |
| ETS (seleção automática) | R$ 1.689.251 | 87,60% |
| Prophet | R$ 4.137.290 | 272,05% |

O Prophet teve desempenho muito abaixo dos demais por extrapolar a tendência de crescimento observada entre 2021 e 2023, sem reconhecer a estabilização da série a partir de 2024. O ETS ficou mais próximo, mas ainda superestimou os custos de forma consistente. O SARIMA foi o único que capturou tanto a ausência de tendência crescente no período recente quanto o padrão sazonal de janeiro.

## Previsão (nov/2025 a out/2026)

O modelo projeta a manutenção do padrão histórico: queda acentuada em janeiro de 2026 (previsão de R$ 146 mil), retomada gradual ao longo do primeiro semestre e estabilização entre R$ 800 mil e R$ 1 milhão no segundo semestre, bem abaixo dos picos de 2023, consistente com a estabilização observada desde 2024.

## Estrutura do repositório

```
.
├── tratamento_dados.R           # Limpeza e tipagem da matriz transposta
├── script_analise.R             # Modelagem e previsão da série temporal
├── analisefinal_MatheusEdu.pdf  # Relatório completo do trabalho
└── README.md
```

> Os arquivos de dados (`matriz_transposta.csv` e `dados_prefeitura_limpos_final.csv`) não estão versionados no repositório. Os dados brutos estão disponíveis no Portal da Transparência da Prefeitura de Uberlândia.

## Dependências

```r
library(tidyverse)
library(forecast)
library(tseries)
library(urca)
library(lubridate)
library(deflateBR)
library(prophet)
library(knitr)
```

## Referências

- Box, G. E. P.; Jenkins, G. M.; Reinsel, G. C. *Time Series Analysis: Forecasting and Control*. 5. ed. Wiley, 2015.
- Hyndman, R. J.; Athanasopoulos, G. *Forecasting: Principles and Practice*. 3. ed. OTexts, 2021. Disponível em: https://otexts.com/fpp3/
- Taylor, S. J.; Letham, B. Forecasting at Scale. *The American Statistician*, v. 72, n. 1, p. 37–45, 2018.
