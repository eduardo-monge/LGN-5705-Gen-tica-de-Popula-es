💻 Prática Dia 1: Organização de Projeto e Exploração de Dados

Nesta primeira prática, vamos configurar nosso ambiente de trabalho e aprender comandos essenciais para manipular e explorar arquivos de sequenciamento.


## 0. Criando o ambiente de CONDA com todos os programas que vamos a usar

``` bash
1. Descarregar o snp_calling_env.yaml e deixar dentro da pasta de UBUNTU. 

2. Rodar os ambientes
conda env create -f snp_calling_env.yaml
conda activate snp_calling
```
---

## 1. Configuração da Estrutura de Diretórios

A organização é o primeiro passo para a reprodutibilidade. Vamos criar uma hierarquia lógica para o nosso projeto de bioinformática. Vamos a seguir a lógica recomendada por [Noble (2009) "A Quick Guide to Organizing Computational Biology Project"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

<img width="4134" height="2308" alt="Noble (2009)" src="https://github.com/user-attachments/assets/84830a55-885a-4b5f-8e4e-014de6ded742" />

```bash
# Navegar para o diretório pessoal
cd ~

# Criar a pasta raiz do curso (caso ainda não tenha criado)
mkdir -p ~/workshop_bioinfo
cd ~/workshop_bioinfo

#Se quiser fazer tudo numa linha só
mkdir -p ~/workshop_bioinfo && cd ~/workshop_bioinfo

# -----------------------------------------
#A opção -p significa “pais”.
#Ela instrui o sistema a:
#1️⃣ Criar todos os diretórios pais, caso eles não existam
#2️⃣ Evitar gerar um erro se o diretório já existir
# -----------------------------------------

# Criar subpastas estruturadas
mkdir -p data/raw data/trimmed data/mapped
mkdir -p results/qc results/snp_calling results/vcf_filtered
mkdir -p scripts #Opcional para o tutorial, recomendado para seus projetos pessoais. 

# Verificar a estrutura criada
ls -R

# -----------------------------------------
#A opção -R significa “recursive”.
#Ele instrui o comando ls a listar o conteúdo dos diretórios e todos os seus subdiretórios recursivamente.
# -----------------------------------------
```
## 2. Download de Dados de Sequenciamento (FASTQ)
Vamos baixar um conjunto de dados reduzido (subset) para praticar os comandos de manipulação. Normalmente utilizamos o comando `wget` (ou `curl`) para realizar o download de conjuntos de dados e genomas de bancos públicos diretamente da internet para os nossos servidores. Ele é incrivelmente eficiente para links estáticos. 

No entanto, quando nossos dados estão hospedados em ecossistemas de armazenamento em nuvem corporativos, como o Google Drive, o `wget` tradicional falha. Isso ocorre porque o Google não fornece um link direto simples para o arquivo bruto; a plataforma utiliza protocolos dinâmicos de segurança, exigindo a resolução de *cookies*, tokens de sessão e, frequentemente, confirmações manuais de antivírus que bloqueiam requisições simples de terminal.

Como alternativa, utilizamos o comando **`gdown`**. O  [gdown](https://github.com/wkentaro/gdown) é um pacote de código aberto desenvolvido em Python, projetado especificamente para interagir com a API do Google Drive emulando o comportamento de um navegador web. Quando você fornece o link da pasta, o algoritmo do `gdown` "conversa" com os servidores do Google, intercepta automaticamente os *cookies* de confirmação de segurança em segundo plano e resolve os redirecionamentos de página. O resultado é a capacidade de espelhar o conteúdo exato da sua pasta do Drive diretamente para o diretório do seu servidor de alto desempenho, contornando os bloqueios de segurança sem a necessidade de configurações complexas.


```bash
# Entrar na pasta de dados brutos
cd data/raw

# Baixar arquivos FASTQ de exemplo desde Google Drive
gdown --folder https://drive.google.com/drive/folders/1HMaZB53tq5jiQbg0tE2Alk8P0WF1Afmt?usp=drive_link

# Listar os arquivos e verificar o tamanho (legível para humanos)
ls -lh
```

💡 **NOTA DO DIA A DIA: Como recebo meus dados?**
> Na vida real, quando o sequenciamento das suas amostras finaliza, você não recebe os dados por pendrive. A empresa geralmente disponibiliza os arquivos em servidores FTP ou na nuvem oficial da [Illumina BaseSpace](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub.html). 
> 
> Nesses casos, você usaria o terminal para baixar os dados diretamente para o seu servidor. 
> 
> **Exemplos de comandos reais (Não execute agora):**
> ```bash
> # Opção A: Usando a ferramenta oficial da Illumina (BaseSpace CLI)
> # (Baixa o projeto inteiro de uma vez)
> bs download project --id 12345678 --output ./data/raw/
> 
> # Opção B: Usando wget com um link de servidor FTP seguro do laboratório
> # (O '-c' permite continuar o download de onde parou caso a internet caia)
> wget -c "ftp://usuario:senha@ftp.centro-sequenciamento.com/projeto/*.fastq.gz"
> ```


## 3. Exploração de Dados via Terminal
Arquivos de bioinformática costumam ser muito grandes. Aprender a "espiar" o conteúdo sem abrir o arquivo inteiro é fundamental.

A. Anatomia do FASTQ (Comando head e tail)
Cada read no formato FASTQ possui 4 linhas. Vamos ver o primeiro read do arquivo:
```bash
zcat AM_04.fq.gz | head -n 4
```
❓ O que significa cada parte do cabeçalho? 

<details>
<summary>🧬 O Cabeçalho FASTQ</b> </summary>

<br>

Este cabeçalho é o identificador único para cada uma única molécula de DNA genômico. Ele contém as coordenadas físicas exatas de onde essa molécula foi processada no equipamento (geralmente uma plataforma Illumina):

* **`@`**: O caractere obrigatório que indica o início de um bloco de sequência no formato FASTQ.
* **`152`**: O Identificador do Sequenciador (*Instrument ID*) ou o número da corrida (*Run ID*). É o "nome" da máquina que realizou o sequenciamento.
* **`1`**: O número da Pista (*Lane*). A lâmina de vidro do sequenciador (*Flowcell*) é dividida em canaletas (frequentemente 8). Esta amostra correu na canaleta 1.
* **`1101`**: O "Azulejo" (*Tile*). Cada canaleta é subdividida em milhares de quadradinhos minúsculos chamados *tiles*, onde o DNA se fixa e é amplificado (formando os *clusters*).
* **`29858`**: A coordenada espacial X (posição horizontal) exata do *cluster* de DNA dentro daquele azulejo específico.
* **`1000`**: A coordenada espacial Y (posição vertical) exata do *cluster* dentro do azulejo.
* **`/1`**: O tipo de leitura (*Read Type*). Em metodologias de sequenciamento modernas (*Paired-End*), lemos o fragmento de DNA pelas duas pontas. O `/1` indica que esta é a leitura *Forward* (Read 1). O arquivo parceiro dessa mesma amostra terá coordenadas idênticas, mas terminará com `/2` (a leitura *Reverse*).

> 💡 **Nota Tecnológica:** O sequenciador necessita registrar meticulosamente as posições X e Y (`29858_1000`) de cada molécula pois ele "lê" o DNA capturando fotografias em altíssima resolução da lâmina a cada ciclo químico de adição de bases. Ele utiliza essas coordenadas matemáticas para garantir que o ponto fluorescente piscando naquela exata posição seja a mesma molécula de DNA em crescimento!

</details>

B. Contagem de Sequências
Podemos usar o contador de linhas (wc -l) para estimar o número de reads.
```bash
# Contar total de linhas
zcat AM_04.fq.gz | wc -l
```

C. Busca com Grep
O comando grep é poderoso para encontrar padrões. Vamos verificar os cabeçalhos das sequências.
```bash
# Listar os primeiros 5 cabeçalhos (começam com @)
grep "^@" AM_04.fq.gz | head -n 5
```

## 4. Gestão de Metadados (Prática com Nano)
Agora, vamos criar nossa tabela de metadados para vincular as amostras às suas informações biológicas.

```bash
# Voltar para a pasta data
cd ..

# Criar o arquivo de metadados
nano metadata.tsv
```
Dentro do editor Nano, digite o seguinte (use a tecla TAB para separar as colunas):

```bash
AM_04	amazonas
AM_05	amazonas
G109	guatemala
G111	guatemala
rifania10	sao_paulo
rifania8	sao_paulo
luz_20m	minas
luz_25m	minas
veracruz-mex-L1p2	mexico
veracruz-mex-L3p2	mexico
```
(Para salvar: Ctrl + O, depois Enter. Para sair: Ctrl + X)



💡 **Dica: Usando o Excel**
> 
> Trabalhar com editores de texto dentro do terminal, como o `nano`, pode ser um desafio no início, especialmente para quem ainda não está acostumados com ambientes de programação. Não se preocupe! 
> 
> Uma alternativa muito prática é organizar sua tabela de metadados no Excel e, em seguida, simplesmente copiar as células e colar dentro do terminal. Se você estiver trabalhando com o WSL (*Windows Subsystem for Linux*), essa integração de "copiar e colar" entre o Windows e o Linux é direta, rápida e facilita muito a vida.
> 
> Só lembra que no terminal do Linux, os atalhos clássicos `Ctrl+C` e `Ctrl+V` **não funcionam** (inclusive, o `Ctrl+C` serve para cancelar e abortar comandos em andamento). Para colar o texto que você copiou do Excel, você precisará **clicar com o botão direito do mouse** dentro da tela do terminal.
> 

**⚠️ Regras de Ouro para Nomes (Arquivos e Metadados):**
> Para evitar que os programas de bioinformática apresentem erros misteriosos no futuro, siga sempre estas três regras na hora de preencher sua tabela:
> 
> 1. **Sem caracteres especiais:** Evite símbolos como `@`, `#`, `$`, `%`, `&`, etc.
> 2. **Sem acentuação ou cedilha:** Não utilize acentos (´, ^, ~, `) ou a letra 'ç' (ex: escreva *populacao* em vez de *população*).
> 3. **Sem espaços em branco:** O Linux entende o espaço como o fim de um comando. Para separar palavras, utilize sempre o sublinhado/underscore (`_`) (ex: `amostra_populacao_01`).
