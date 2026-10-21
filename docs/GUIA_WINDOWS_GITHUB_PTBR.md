# Guia Windows / RStudio / GitHub — artigo REEC

## 1. Extrair o pacote

Extraia `REEC_Cachoeira_Forest_Trajectories_GitHub.zip` em `C:\Projetos\`.
A estrutura ficará em `C:\Projetos\REEC_Cachoeira_Forest_Trajectories\`.

Não é necessário Docker, PostgreSQL ou ferramentas pagas.

## 2. Executar no RStudio

Abra `REEC_Cachoeira_Forest_Trajectories.Rproj`. Em seguida execute:

```r
source("scripts/01_reproduce_article.R", encoding="UTF-8")
source("scripts/02_validate_against_article.R", encoding="UTF-8")
```

Aguarde a instalação inicial dos pacotes R, se o RStudio solicitar autorização. Se tudo estiver correto, o console mostrará `ALL REEC NUMERIC VALIDATIONS PASSED`. Os resultados serão salvos na pasta `outputs`.

O arquivo `outputs/sessionInfo.txt` registra as versões de R e de pacotes utilizadas. Compare-o com a declaração do manuscrito.

## 3. Criar um repositório público exclusivo do artigo

Entre em https://github.com/new, autentique-se na sua conta `ernestoreis` e preencha:

- Repository name: `REEC_Cachoeira_Forest_Trajectories`
- Description: `Reproducible data and R analysis for forest-cover trajectories in the Cachoeira River basin (1985–2022)`
- Visibility: Public
- **Não** marque a criação automática de README, .gitignore ou licença (o pacote já tem os dois primeiros).

Selecione **Create repository**. A página exibirá o endereço do repositório recém-criado.

## 4. Publicar com PowerShell

Abra o PowerShell e execute (ajuste apenas a pasta se extraiu em outro local):

```powershell
cd "C:\Projetos\REEC_Cachoeira_Forest_Trajectories"
git --version
git init
git branch -M main
git add .
git commit -m "Initial reproducibility package for REEC article"
git remote add origin https://github.com/ernestoreis/REEC_Cachoeira_Forest_Trajectories.git
git push -u origin main
```

Se `git commit` solicitar nome e e-mail, configure-os para a identidade que deseja associar ao projeto:

```powershell
git config --global user.name "Ernesto Santana dos Reis"
git config --global user.email "SEU_EMAIL_NO_GITHUB"
```

Em seguida repita o `git commit` e o `git push`. O GitHub pode solicitar autenticação pelo navegador ou por token; não insira senhas dentro dos scripts ou arquivos do repositório.

## 5. Verificar a publicação

Acesse https://github.com/ernestoreis/REEC_Cachoeira_Forest_Trajectories e confira as pastas `data/raw`, `scripts`, `reference` e `figures/published`.

Atenção: o URL só funcionará depois que o repositório for criado e o primeiro `git push` for concluído.

## 6. Atualizar a declaração Data availability no artigo

Após confirmar que o repositório está público, substitua a declaração anterior por:

> MapBiomas Collection 8 provides the original annual land-use and land-cover data. The four aggregated annual percentage datasets used for statistical analysis (1985–2022), the R scripts for reproducing the time-series analyses, principal component analysis, and Figures 3–5, and the supporting numerical results are available at https://github.com/ernestoreis/REEC_Cachoeira_Forest_Trajectories. The source raster and GIS files required to regenerate the cartographic Figures 1–2 are not included in this repository.

Não afirme que a reconstrução integral da cartografia está incluída até acrescentar os insumos espaciais e rotinas GIS.

## 7. Atualizações posteriores

Para acrescentar um arquivo ou corrigir a documentação:

```powershell
cd "C:\Projetos\REEC_Cachoeira_Forest_Trajectories"
git status
git add .
git commit -m "Update REEC reproduction materials"
git push
```

Evite publicar `.RData`, credenciais, dados pessoais, ambientes temporários ou versões intermediárias desnecessárias.
