# 💻 Prática Dia 2: Demultiplexação, QC e Trimming

Nesta prática, vamos aprender a analisar as sequências para verificar a qualidade e fazer o trimming necessário. 

---
## 1. Demultiplexação com STACKS (`process_radtags`)
Quando sequenciamos bibliotecas de representação reduzida (como RADseq ou GBS), várias amostras são misturadas no mesmo tubo (multiplexadas). O primeiro passo é usar os "códigos de barras" (*barcodes*) para separar quem é quem. Geralmente, as empresas já fornecem os dados separados, mas também podemos fazê-lo nós mesmos. 

Existem diferentes programas e métodos, dependendo da biblioteca genômica que usamos. No caso do GBS/RAD-seq, um programa que pode ser utilizado é o módulo `process_radtags` do [STACKS](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).

### 0. Preparar os dados para o programa STACKS (El archivo de Barcodes)
Antes de rodar qualquer comando do STACKS para demultiplexar, o programa precisa de um "mapa" para saber como separar essa "sopa" de sequências. Para que o processo seja um sucesso, precisamos de três elementos fundamentais:

1. **A enzima de restrição:** Saber exatamente qual enzima (ou enzimas) foi usada na bancada do laboratório para cortar o DNA (ex: *EcoRI*, *MspI*, *PstI*).
2. **Os arquivos FASTQ brutos:** Que já baixamos na nossa pasta `data/raw`.
3. **O arquivo de Barcodes:** Um "mapa" de texto que vincula cada código de barras com o nome real da amostra.

### 📝 Como deve ser a planilha de amostras?

O STACKS é extremamente rigoroso com o formato deste arquivo. Ele deve ser um arquivo de texto simples (`.txt`), onde as colunas são separadas **exclusivamente pela tecla TAB** (tabulação). Não deve ter cabeçalhos (títulos de colunas) nem espaços em branco soltos.

A estrutura mais básica tem duas colunas: `[Sequência_do_Barcode]<tab>[Nome_da_Amostra]`

<img width="223" height="222" alt="image" src="https://github.com/user-attachments/assets/427204af-e5fb-4f8d-afae-a6930986df85" />

Os nomes das amostras podem ser os relevantes para o seu projeto. Lembre-se apenas das dicas para nomear amostras aprendidas no dia 1. 
Vamos criar nosso arquivo de *barcodes* diretamente no terminal:

```bash
# Navegar para a pasta principal de dados
cd ~/workshop_bioinfo/data

# Criar e abrir o arquivo com o editor nano
nano barcodes.txt

#Dentro do nano, você pode copiar e colar os códigos de barras e amostras de um arquivo Excel ou TXT do seu computador. Lembre-se de como fizemos no dia 1.

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
(Para salvar: pressione Ctrl + O e depois Enter. Para sair: Ctrl + X)

⚠️ **ATENÇÃO**
>Lembre-se de verificar se cada código de barras corresponde realmente à amostra. Sempre verifique e mantenha seu mapa da placa em ordem. 


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


## 3. Limpeza de Dados (Trimming)
Como vimos nos relatórios, temos algumas leituras de baixa qualidade, bem como a presença de adaptadores Illumina e cola poli-G. Para eliminar essas partes, vamos usar duas ferramentas diferentes muito utilizadas na bioinformática. 

### A. Corte e filtragem com Fastp
O [Fastp](https://github.com/OpenGene/fastp) é uma ferramenta eficiente para controle de qualidade e pré-processamento. Ele pode:
1. Cortar bases de baixa qualidade.
2. Remover adaptadores.
3. Filtrar reads curtas.
4. Gerar relatórios de controle de qualidade.

O `fastp` faz tudo de uma vez: onsegue detectar automaticamente os adaptadores Illumina, as caudas poli-X, corta adaptadores, filtra qualidade e já gera um relatório visual. 

```bash
# Ir para a pasta de trimmed
cd ~/workshop_bioinfo/data
mkdir trimmed
cd trimmed

# Rodar fastp para um para uma amostra
fastp -i ../demultiplex/AM_04.fq.gz \
      -o  AM_04.fq.gz_trimmed.fq.gz \
      -q 20 \
      --trim_poly_g \
      --trim_poly_x \
      -h relatorio_fastp.html

# -----------------------------------------
#Os parâmetros importantes que estamos usando: 
# -i: nome do input
# -o: nome do output
# --trim_poly_g: elimina a cauda polyG
# -q ou --qualified_quality_phred: o valor de qualidade que uma base é qualificada.
# --html: nome do report final no formato html
# -----------------------------------------
```
❓Verifique o relatório HTML gerado pelo Fastp e responda às seguintes perguntas: Qual foi o número de reads antes e depois da filtragem? Qual foi a porcentagem de adaptadores detectados e filtrados? Qual foi a porcentagem de reads com baixa qualidade, reads contendo muitas bases N e reads muito curtas?

🤖 Automatizando o Processo: Rodando o Fastp em todas as amostras

Na vida real, você terá dezenas ou centenas de amostras. Rodar o comando manualmente para cada uma seria impossível. Para resolver isso, usaremos um **Loop `for`** do Bash. Esse comando vai percorrer todos os arquivos da pasta e aplicar o `fastp` em cada um deles automaticamente.

```bash
# 1. Certifique-se de que está na pasta correta
cd ~/workshop_bioinfo/data/demultiplex

# 2. Criar a pasta de saída para os arquivos limpos
mkdir -p ../trimmed_fastp

# 3. Rodar o Loop para arquivos Single-End
for file in *.fq.gz; do
    # Criar um nome para o arquivo de saída adicionando '_trimmed'
    OUTNAME=$(basename $file .fq.gz)_trimmed
    
    echo "Processando a amostra: $file"

    fastp -i $file \
          -o ../trimmed/$OUTNAME.fq.gz \
          -q 20 \
          --trim_poly_g \
          --trim_poly_x \  
          -h ../trimmed/$(basename $file .fq.gz)_report.html
done


# -----------------------------------------
#Os parâmetros importantes que estamos usando: 
# for file in *.fastq; do: Diz ao Linux: "Para cada arquivo que termina com .fastq, guarde o nome dele temporariamente na variável chamada file e comece o processo".
# basename $file .fastq: Esse comando é um truque para pegar apenas o nome da amostra e remover a extensão .fastq. Isso nos permite criar nomes de saída organizados.
# $file: O cifrão indica que estamos chamando o conteúdo da variável (o nome do arquivo atual).
# done: Avisa ao Linux que o comando acabou para aquela amostra e que ele pode passar para a próxima da lista.
# -----------------------------------------
```
Após rodar o loop do `fastp`, você terá um novo arquivo `.fastq` para cada amostra dentro da pasta `trimmed`. Agora, precisamos confirmar se a limpeza foi eficiente. Para isso, geraremos novos relatórios **FastQC** e **MultiQC** para comparar com os dados brutos.
1. Crie uma nova pasta chamada `trimmed_fastp` dentro do seu diretório de resultados de QC. Isso é importante para não misturar os relatórios dos dados brutos com os relatórios dos dados limpos. (Dica: Use o comando `mkdir`).
2. Execute o **FastQC** apenas nos arquivos que estão dentro da pasta `trimmed`. Lembre-se de usar o parâmetro para salvar os resultados dentro da pasta que você criou no Passo 1. (Dica: Use o caractere curinga `*` para selecionar todos os arquivos de uma vez e o parâmetro `-o` para a saída).
3. Entre na pasta `trimmed_fastp` e execute o **MultiQC** para agrupar todos os novos relatórios em um único arquivo interativo. (Dica: O comando é executado dentro da pasta onde os arquivos do FastQC foram gerados).

❓ O que devemos observar agora?
Ao abrir o novo arquivo multiqc_report.html, verifique os seguintes pontos:
1. **Qualidade por Base:** Os "bigodes" dos gráficos de qualidade devem estar agora predominantemente na zona verde.
2. **Presença de Adaptadores:** O gráfico de Adapter Content deve estar completamente limpo (linha plana em 0%).
3. **Comprimento das Sequências:** Se você usou filtros de tamanho, verá uma variação no comprimento das reads (elas não terão mais todas o mesmo tamanho exato).
4. **Caudas de Poli-G:** Se as amostras tinham esse problema, o conteúdo de Guanina (G) no final das leituras deve ter caído drasticamente.

### B. Corte e filtragem com TrimGalore (Alternativa)
O `fastp` é uma das ferramentas mais utilizadas hoje por ser muito rápida, mas não é a única. Às vezes, dependendo do conjunto de dados, ele pode não resolver nosso problema ou remover coisas demais. Outra ferramenta clássica e muito respeitada na comunidade é o [Trim Galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md). Ele funciona como um "facilitador" (wrapper) que roda dois programas famosos ao mesmo tempo: o `cutadapt` (para limpar) e o `fastqc` (para gerar o relatório logo em seguida).

A. Criar as pastas novas
Como vamos testar um programa novo, não queremos misturar os resultados do Trim Galore com os do fastp. Vamos criar uma pasta só para ele.

```bash
# Vá para a pasta onde estamos guardando os dados limpos
cd ~/workshop_bioinfo/data/trimmed

# Crie uma nova pasta específica para o Trim Galore
mkdir -p trimmed_trimgalore
```
B. Rodando para UMA amostra (Single-End)
Vamos testar o comando em uma única amostra primeiro para entender os parâmetros. Como nossas sequências são Single-End (uma única leitura por fragmento), o comando é bem simples:

```bash
# Volte para a pasta com os arquivos demultiplexados
cd ~/workshop_bioinfo/data/processed/demultiplex

# Rodar o Trim Galore em apenas um arquivo
trim_galore \
    --nextseq 20 \
    --illumina \
    --length 20 \
    --fastqc \
    --output_dir ../trimmed/trim_galore \
    AM_04.fq.gz 

# -----------------------------------------
#Os parâmetros importantes que estamos usando: 
--nextseq 20: Identifica e elimina  poly-G e elimina reads com qualidade Q < 20
--illumina: Indica o uso de adaptadores padrão da Illumina.
--length 20: Elimina leituras que ficam muito curtas após o corte.
--fastqc: Executa automaticamente o FastQC após o corte.
--output_dir: Diretório onde os arquivos são salvos.
# -----------------------------------------
```
Verifique o relatório gerado automaticamente pelo programa. 
❓Como esse relatório se compara aos resultados do fasp para a mesma amostra? 

C. Automatizando com um Loop (Para todas as amostras)
Agora vamos usar o mesmo poder da automação que aprendemos antes para processar todos os arquivos FASTQ da pasta de uma só vez.

```bash
# Certifique-se de que está na pasta demultiplex
cd ~/workshop_bioinfo/data/processed/demultiplex

# Rodar Trim Galore para todos os arquivos single-end
for file in *.fq.gz ; do
    echo "Limpando a amostra com Trim Galore: $file"
    trim_galore \
        --nextseq 20 \
        --illumina \
        --length 20 \
        --fastqc \
        --output_dir ../trimmed/trim_galore \
        "$file"

done
```
💡 Nota Prática: O Trim Galore já adiciona automaticamente o sufixo _trimmed.fq no final do nome do arquivo, então não precisamos daquele truque do basename que usamos no fastp. Ele faz a organização dos nomes sozinho!

Agora, basta executar o MultiQC para obter um resumo de tudo o que foi feito pelo TrimGalore.
❓Como os resultados do fasp se comparam aos do multiQC? Há grandes mudanças?


---
💡 Dicas: Melhores práticas para a qualidade dos dados
1. Retenha os metadados: Mantenha registros e resumos do FastQC e Fastp para reproducibilidade.
2. Otimize os parâmetros: personalize as configurações do Fastp/TrimGalore (por exemplo, limites de corte) com base nos relatórios do FastQC.
3. Processamento em lote: use loops ou job arrays para grandes conjuntos de dados.
4. Salve os outputs: mantenha um diretório estruturado com dados brutos, aparados e com qualidade verificada.

Seguindo esse fluxo de trabalho, você garante que seus dados de sequenciamento sejam de alta qualidade e estejam prontos para análise posterior. Essas ferramentas e técnicas são fundamentais na genômica, ajudando a maximizar a confiabilidade e a precisão dos seus resultados.
