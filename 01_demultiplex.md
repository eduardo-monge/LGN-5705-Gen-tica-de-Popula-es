# 💻 1: Demultiplexação e controle de qualidade QC

Nesta prática, vamos aprender a analisar as sequências para verificar a qualidade. 

---
## 1. Demultiplexação com STACKS (`process_radtags`)
Quando sequenciamos bibliotecas de representação reduzida (como RADseq ou GBS), várias amostras são misturadas no mesmo tubo (multiplexadas). O primeiro passo é usar os "códigos de barras" (*barcodes*) para separar quem é quem. Geralmente, as empresas já fornecem os dados separados, mas também podemos fazê-lo nós mesmos. 


Existem diferentes programas e métodos, dependendo da biblioteca genômica que usamos. No caso do GBS/RAD-seq, um programa que pode ser utilizado é o módulo `process_radtags` do [STACKS](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).

### 0. Preparar os dados para o programa STACKS (El archivo de Barcodes)
Antes de rodar qualquer comando do STACKS para demultiplexar, o programa precisa de um "mapa" para saber como separar essa "sopa" de sequências. Para que o processo seja um sucesso, precisamos de três elementos fundamentais:

1. **A enzima de restrição:** Saber exatamente qual enzima (ou enzimas) foi usada na bancada do laboratório para cortar o DNA (ex: *EcoRI*, *MspI*, *PstI*).
2. **Os arquivos FASTQ brutos:** Que já baixamos na nossa pasta `/home/genpop/dados_brutos`.
3. **O arquivo de Barcodes:** Um "mapa" de texto que vincula cada código de barras com o nome real da amostra.

### 📝 Como deve ser a planilha de amostras?

O STACKS é extremamente rigoroso com o formato deste arquivo. Ele deve ser um arquivo de texto simples (`.txt`), onde as colunas são separadas **exclusivamente pela tecla TAB** (tabulação). Não deve ter cabeçalhos (títulos de colunas) nem espaços em branco soltos.

A estrutura mais básica tem duas colunas: `[Sequência_do_Barcode]<tab>[Nome_da_Amostra]`

<img width="223" height="222" alt="image" src="https://github.com/user-attachments/assets/427204af-e5fb-4f8d-afae-a6930986df85" />

Os nomes das amostras podem ser os relevantes para o seu projeto. A dica é usar sempre nomes simples, alfanuméricos. 

## Criar o arquivo de barcodes

Primeiro, vamos para a nossa pasta principal e criamos uma pasta para a demultiplexação:

```bash
# Ir para a minha pasta principal
cd

# Criar a pasta "demultiplex" e entrar nela
mkdir -p demultiplex && cd demultiplex

# Criar e abrir o arquivo com o editor nano
nano barcodes.txt
```

Dentro do `nano`, cole os códigos de barras e os nomes das amostras. Você pode copiar de um arquivo Excel ou TXT do seu computador, como fizemos no dia 1. Cada linha deve ter o **barcode**, um **TAB** e o **nome da amostra**:

```
TGACGCCA	AM_04
CAGATA	AM_05
CTCGCGG	G109
AACTGG	G111
ACGCGCG	luz_20m
GTCGCCT	luz_25m
GGACAG	rifania10
ATCTGT	rifania8
TCAGAGAT	veracruz-mex-L1p2
CGTTCA	veracruz-mex-L3p2
```

Para salvar: `Ctrl + O` e depois `Enter`. Para sair: `Ctrl + X`.

### Conferir o arquivo

O `process_radtags` exige que o barcode e o nome estejam separados por um **TAB**, e não por espaços. Como TABs e espaços são invisíveis na tela, vamos usar o `cat -A`, que mostra os caracteres escondidos:

```bash
cat -A barcodes.txt
wc -l barcodes.txt
```

Cada linha deve aparecer assim:

```
TGACGCCA^IAM_04$
```

- `^I` é o **TAB**: tem que aparecer entre o barcode e o nome.
- `$` é o **fim da linha**: não pode ter nada entre o nome e o `$`.
- Se aparecer `^M$` no final, o arquivo veio com quebra de linha do Windows. Corrija com:

```bash
sed -i 's/\r$//' barcodes.txt
```

O `wc -l` deve mostrar **10** linhas, uma por amostra.

> [!WARNING]
> Confira se cada código de barras corresponde realmente à amostra certa. Um barcode trocado não gera nenhum erro: o programa roda normalmente, mas as leituras vão para a amostra errada. Sempre verifique e mantenha o seu mapa da placa em ordem.


### 1. Executar `process_radtags`
Agora estamos prontos para executar o programa. A documentação completa com todos os parâmetros e casos está disponível na página do [STACKS](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).

```bash
# 1. Navegar para a pasta onde criaremos os resultados demultiplexados
cd ~/workshop_bioinfo
mkdir -p data/demultiplex
cd data/demultiplex

# 2. Executar o STACKS (Substitua as informações conforme sua enzima e arquivos)
# -1 e -2: Arquivos brutos (Forward e Reverse)
# -b: Arquivo de texto com os barcodes
# -o: Pasta de saída
# -e: Enzima de restrição usada
# -r, -c, -q: Limpar reads órfãos, corrigir barcodes e checar qualidade básica

process_radtags -f ../raw/mz_gbs_acrocomia.fastq.gz \
                -b ../barcodes.txt \
                -o . \
                --adapter-1 AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG \
                -e pstI \
                --adapter-mm 2 \
                -r \
                -c \
                -q

# -----------------------------------------
O que este código faz:
# -f (File): Indica ao algoritmo que estamos inserindo um arquivo Single-End e aponta o caminho absoluto ou relativo para as leituras brutas comprimidas.
# -b (Barcodes): Caminho para o arquivo de texto tabular contendo os identificadores (índices) para a separação dos indivíduos.
# -o (Output): Diretório de saída onde os arquivos limpos e separados serão salvos (o . instrui o programa a salvar no diretório atual).
# --adapter-1: Fornece a sequência exata do adaptador Illumina. O programa fará uma varredura rigorosa para encontrar e cortar (trim) essa sequência indesejada.
# -e mseI (Enzyme): Informa a enzima de restrição utilizada na digestão do genoma durante a preparação da biblioteca GBS (neste estudo de caso, a MseI). O algoritmo verificará a integridade deste sítio de corte.
# --adapter-mm 2: Permite até 2 mismatches (erros de pareamento) ao procurar a sequência do adaptador, tornando a limpeza mais rigorosa diante de possíveis erros de sequenciamento.
# -r (Rescue): Ativa o algoritmo de resgate. Ele tenta recuperar barcodes e sítios de enzima que falharam na identificação primária devido a erros mínimos do sequenciador.
# -c (Clean): Habilita o filtro de limpeza básica, descartando sumariamente leituras que possuam bases não identificadas (marcadas como N pelo sequenciador).
# -q (Quality): Habilita o filtro avançado de qualidade. Descartará leituras cujo Phred score caia abaixo do limite mínimo aceitável em uma janela deslizante matemática.
```

💡 **Single-end sequencing versus paired-end sequencing**
> Em nosso exemplo, usamos uma biblioteca sequenciada single-end. Caso você tenha uma biblioteca paired-end, o comando a seguir ajuda a indicar onde se encontram as leituras R1 e R2. Os demais parâmetros são semelhantes.
>
```bash
process_radtags -1 ../data/raw/sample_R1.fastq -2 ../../data/raw/sample_R2.fastq
 ```

## 2. Controle de Qualidade Inicial (QC)
Agora que temos nossos arquivos FASTQ separados por amostra, precisamos fazer a avaliação dos arquivos para garantir dados de alta qualidade para análises posteriores. 

Lembre-se de que os arquivos são arquivos FASTQ, que são arquivos FASTA modificados para incluir também informações sobre a qualidade do sequenciamento para cada base (PHRED score). A qualidade Q reflete a probabilidade p de que a base esteja incorreta. Os valores de Q podem variar de 0 a 93 (mas a pontuação máxima que você normalmente verá é 40). Para economizar espaço em um arquivo FASTQ, cada pontuação de qualidade é transformada em um único caractere ASCII.
<img width="877" height="267" alt="5d0a8f1b-4a1d-4777-b8c9-a4325d242560" src="https://github.com/user-attachments/assets/8d1eb519-f2a2-481c-90bc-fbfa91fa0b23" />

### 0. Visualizando os arquivos FASTQ
Antes de gerar os arquivos de qualidade, vamos usar comandos básicos do UNIX para explorar os dados e nos familiarizar com o formato FASTQ.

A. Primeiro use o os comandos `cat`/`head` para inspecionar um arquivo fastq.
❓Quais são os caracteres iniciais do cabeçalho para todas as leituras?

```bash
zcat ../demultiplex/Acrocomia_pop1_01.fastq | head

# -----------------------------------------
#O símbolo de | é o "Pipe".
#Ele instrui o terminal a pegar o resultado (a saída) do comando que está à sua esquerda e "injetá-lo" diretamente como entrada para o comando que está à sua direita, sem precisar salvar um arquivo temporário no meio do caminho. Isso permite realizar mais de uma operação ao mesmo tempo. 
# -----------------------------------------

```
B. Assim como o `head` lê o topo, o `tail` lê o final do arquivo. As últimas sequências sequenciadas em uma corrida costumam ter qualidade pior.

```bash
tail  ../demultiplex/Acrocomia_pop1_01.fastq
tail -n 15 ../demultiplex/Acrocomia_pop1_01.fastq

# -----------------------------------------
#Tanto head quanto tail retornam 10 linhas. Mas é possível usar o -n para indicar o número de linhas que desejo avaliar.
# -----------------------------------------
```

C. Agora vamos usar o comando `grep` para procurar palavras ou expressões específicas dentro do nosso arquivo.
❓Quantas leituras temos para cada indivíduo?

```bash
grep -c "^@SR" ../demultiplex/Acrocomia_pop1_01.fastq

# -----------------------------------------
#O `-c` significa `count` (contar) e o `^` significa 'no início da linha'
# -----------------------------------------
```

### 1. Avaliando a qualidade dos dados com o FastQC
O [FastQC](https://github.com/s-andrews/FastQC?tab=readme-ov-file) é uma ferramenta amplamente utilizada para avaliar a qualidade dos dados de sequenciamento. Ele gera um relatório HTML detalhado que inclui métricas como:

1. **Qualidade por base:** pontuações médias de qualidade ao longo do comprimento das leituras.
2. **Conteúdo GC:** A porcentagem de bases G e C.
3. **Níveis de duplicatas:** Identifica possíveis duplicatas de PCR.
4. **Conteúdo do adaptador:** Verifica se há adaptadores de sequenciamento residuais.

```bash
# Ir para a pasta de resultados de QC
cd ~/workshop_bioinfo
mkdir -p results/qc/raw_qc
cd results/qc/raw_qc

# Rodar o FastQC para todos os arquivos FASTQ demultiplexados (*.fastq)
fastqc -t 4 ../../../data/demultiplex/*.fq.gz -o .

# -----------------------------------------
#Os parâmetros importantes a serem considerados:
# *.fastq o símbolo de asterísco (*) diz ao Linux para selecionar *todos* os arquivos que terminam com `.fastq` dentro daquela pasta. 
# -o é a pasta de saída dos relatórios HTML. O ponto (.) um atalho que significa "esta pasta atual". Portanto, os resultados serão salvos exatamente onde você está agora no terminal.
# -t 4 são os threads. Indica o número de núcleos do processador serem utilizados simultaneamente. 
# -----------------------------------------

```
Abra o arquivo HTML no navegador para examinar as métricas de qualidade.

❓Avalie a qualidade dos arquivos fastq. Qual é o comprimento médio dos reads 1? Qual é a qualidade dos dados? Há alguma queda na qualidade ao longo da sequência? Qual arquivo tem melhor e pior qualidade? Há alguma sequência de adaptador detectada?


### 2. Gerando um único relatório de todas as leituras com o MultiQC
O FastQC é uma ferramenta muito útil e poderosa. Uma limitação é que ele gera um relatório para cada leitura. Quando temos muitas leituras e muitos indivíduos (como geralmente acontece), pode ser difícil abrir e estudar cada relatório individualmente. Para criar um único relatório consensual que seja mais fácil de avaliar, podemos usar o [MultiQC](https://github.com/MultiQC/MultiQC).

O uso é muito simples, basta estar dentro da pasta onde estão todos os arquivos HTML.

```bash
# O ponto (.) indica para o MultiQC procurar na pasta atual
multiqc .
```
Agora abra o relatório final no navegador e vamos avaliar todas as amostras juntas. 

❓Que características da biblioteca GBS poderiam ter causado esses resultados? 
