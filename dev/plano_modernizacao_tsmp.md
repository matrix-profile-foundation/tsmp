# Plano de Modernização: tsmp (Camada Analítica em R Puro)
**Package:** `tsmp` (Time Series Matrix Profile Framework)  
**Data:** 10 de Outubro de 2026  
**Documento Relacionado:** [matrixprofiler: Relatório de Gaps Funcionais](file:///home/franz/workspaces/matrixprofiler/dev/relatorio_gaps_funcionais.md)

---

## 1. Visão Arquitetural e Objetivos

O objetivo desta modernização é transformar o `tsmp` em um **pacote de mineração e análise de séries temporais puramente em R** (`NeedsCompilation: no`), delegando todo o processamento matemático e computacional intensivo para o motor base [`matrixprofiler`](file:///home/franz/workspaces/matrixprofiler).

### Divisão de Responsabilidades:
```
┌──────────────────────────────────────────────────────────────┐
│                            tsmp                              │
│       Camada Analítica, Mineração e Visualização (Puro R)    │
│  - Classes S3 (`MatrixProfile`, `MultiMatrixProfile`, etc.)  │
│  - Descoberta de Motifs, Discords, Snippets, Chains          │
│  - Segmentação Semântica (FLUSS, FLOSS)                      │
│  - Pan Matrix Profile (PMP), MPdist, Regras de Análise       │
│  - Visualizações e Gráficos (`plot`, ggplot2)                │
└──────────────────────────────┬───────────────────────────────┘
                               │ Imports: matrixprofiler
┌──────────────────────────────▼───────────────────────────────┐
│                       matrixprofiler                         │
│           Motor Computacional de Ultra-Performance          │
│  - Algoritmos nativos: MPX, STOMP, STAMP, SCRIMP, MASS, RFCP │
│  - Multithreading C++ via `RcppThread` (POSIX fork-safe)     │
│  - Aceleração por GPU via CUDA (`has_cuda()`)                │
│  - Vetores rápidos de estatística móvel (`windowfunc`)       │
└──────────────────────────────────────────────────────────────┘
```

---

## 2. Diagnóstico das Branches (`develop` vs `master`)

### 2.1 Estado do `develop` (Trabalho prévio de desacoplamento)
- **Eliminação de Código C++:** A pasta `src/` foi completamente removida (`mpx.cpp`, `math.cpp`, `windowfunc.cpp`, etc.).
- **Remoção de rotinas redundantes:** Códigos de algoritmos de baixo nível em R (`mpx.R`, `scrimp.R`, `stamp.R`, `stomp.R`, `mass*.R`) foram removidos.
- **Configuração do Pacote:** `DESCRIPTION` alterado para `NeedsCompilation: no`, removendo `LinkingTo: Rcpp, RcppParallel` e dependências pesadas legadas (`audio`, `doSNOW`, `foreach`, `parallel`, `RJSONIO`).
- **Integração Inicial:** Protótipo de delegação em `compute()` e `dist_profile()` para o `matrixprofiler`.

### 2.2 Patches do CRAN aplicados no `master` (v0.4.16) a serem preservados
- **Compatibilidade R >= 4.0:** Uso de `inherits(x, ...)` no lugar de `class(x) == ...`.
- **Proteção de Opções Globais:** `on.exit(options(...))` em `serialize.R` para restauração segura de `digits`.
- **Exportação S3 Explícita:** Tags `@noRd` e métodos de subsetting (`[` e `head`/`tail`) corrigidos em `subset.R`.
- **Metadados e Exemplos:** Substituição de `\dontrun` por `\donttest` em páginas de documentação para aprovação no CRAN.

---

## 3. Plano de Execução em Fases

### Fase 1: Sincronização e Preparação da Base Limpa
- [ ] Mesclar ou aplicar cherry-pick das correções CRAN de `master` na branch de trabalho (`develop`).
- [ ] Confirmar limpeza total de artefatos de compilação C++ (`src/`, flags de compilação, Makevars).
- [ ] Validar integridade do `DESCRIPTION` (`NeedsCompilation: no`, dependências enxutas).

### Fase 2: Integração com o Motor `matrixprofiler`
- [ ] **Alinhamento do `compute()`:**
  - Padronizar despacho para `matrixprofiler::mpx()`, `stomp()`, `stamp()`, `scrimp()`.
  - Tratamento unificado de parâmetros (`window_size`, `n_workers`, `s_size`).
- [ ] **Integração do `dist_profile()`:**
  - Delegar cálculo de perfis de distância para `matrixprofiler::mass()`.
- [ ] **Harmonização de Convenção de Índices (`NA_integer_` vs `-1L`):**
  - Garantir que o `tsmp` processe uniformemente `NA_integer_` para índices sem vizinho mais próximo, eliminando riscos de indexação acidental `data[-1]`.
- [ ] **Suporte a Séries com Valores Ausentes (`NA`):**
  - Habilitar passagem limpa de séries com `NA` utilizando as estratégias nativas do `matrixprofiler`.

### Fase 3: Modernização da Camada Analítica
- [ ] **Motifs e Discords (`find_motif()`, `find_discord()`):**
  - Consumir matrizes e perfis gerados pelo motor `matrixprofiler`.
  - Otimizar buscas e extração de vizinhos em R sem gargalos.
- [ ] **Segmentação Semântica (FLUSS / FLOSS):**
  - Validar cálculo do Arc Count (CAC) e pontos de transição semântica sobre os perfis do novo motor.
- [ ] **Snippets, Chains e PMP (Pan Matrix Profile):**
  - Revisar e garantir compatibilidade dos algoritmos de sumarização e cadeias de subsequências.
- [ ] **Camada de Visualização (`plot()` e ggplot2):**
  - Assegurar renderização correta de perfis, anotações de motifs e arcos para objetos S3 gerados.

### Fase 4: Suíte de Testes e Conformidade CRAN
- [ ] Migrar testes unitários para a 3ª edição do `testthat` (`testthat::test_package()`).
- [ ] Criar testes cruzados de paridade entre as saídas analíticas do `tsmp` e as saídas nativas do `matrixprofiler`.
- [ ] Executar bateria completa de `R CMD check --as-cran` garantindo `0 errors | 0 warnings | 0 notes`.

---

## 4. Matriz de Integração de Algoritmos

| Função Analítica `tsmp` | Algoritmo Base `matrixprofiler` | Modo Suportado | Multithread | NA Support |
| :--- | :--- | :---: | :---: | :---: |
| `compute(..., algo = "mpx")` | `matrixprofiler::mpx` | Self / AB | ✅ Sim | ✅ Sim |
| `compute(..., algo = "stomp")` | `matrixprofiler::stomp` | Self / AB | ✅ Sim | ⚠️ Validar skip |
| `compute(..., algo = "stamp")` | `matrixprofiler::stamp` | Self / AB | ✅ Sim | ⚠️ Parcial |
| `compute(..., algo = "scrimp")` | `matrixprofiler::scrimp` | Self / AB | ✅ Sim | ⚠️ Parcial |
| `dist_profile()` | `matrixprofiler::mass` | Query Match | ✅ Sim (v3 norm) | ❌ Não |
| `pmp()` | `matrixprofiler::mpx` / `stomp` | Multi-escala | ✅ Sim | Conforme motor |
| `find_motif()`, `find_discord()` | *(Analítico interno `tsmp`)* | Pós-processamento | N/A | N/A |
| `fluss()`, `floss()` | *(Analítico interno `tsmp`)* | Pós-processamento | N/A | N/A |
