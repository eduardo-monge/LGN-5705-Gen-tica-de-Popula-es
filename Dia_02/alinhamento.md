# 💻 Prática Dia 2 (Tarde): Alinhamento e Mapeamento de Sequências
Agora que temos nossos dados limpos e filtrados (seja pelo fastp ou Trim Galore), o próximo passo é o Mapeamento. Vamos comparar cada read da nossa amostra com o Genoma de Referência para descobrir de qual parte do cromossomo aquela sequência veio.

## 0. Instalar CONDA

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```
---
## 1. Genoma de referência
Antes de começar, precisamos conhecer nosso genoma de referência e poder indexá-lo. 

### A. Baixe e conheça o genoma de referência 
O primeiro passo é obter o genoma de referência da Acrocomia aculeata (Acessão: GCA_055471735.1) directamente de NCBI. 

```bash
1. Criar pasta para a referência e entrar nela
cd ~/workshop_bioinfo/data
mkdir -p reference
cd reference
```
Ahora vamos a fazer o download do genoma via FTP (File Transfer Protocol), compactado diretamente do repositório do NCBI.

```bash
# 2. Baixar o genoma da Macaúba (Acrocomia aculeata)
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/055/471/735/GCA_055471735.1_USP_Acracu_1.0/GCA_055471735.1_USP_Acracu_1.0_genomic.fna.gz

gdown --folder https://drive.google.com/drive/folders/1vuZXdhxMmvXq2RY8dAHBtcWsCwPvSG_G?usp=drive_link

# 3. Renomear para facilitar o uso (opcional, mas recomendado)
mv GCA_055471735.1_ASM5547173v1_genomic.fna acrocomia_ref.fna.gz
```

Como o genoma é muito grande, na prática vamos usar apenas o cromossomo 1.
```bash
zcat acrocomia_ref.fna.gz | awk '
/^>/ {
    # Usamos "chromosome 1," com a vírgula no final para garantir que o 
    # awk não pegue o cromossomo 10, 11 ou 12 por acidente!
    if ($0 ~ /chromosome 1,/) {
        extrair = 1
        print $0
    } else {
        extrair = 0
    }
    next
}
extrair == 1 {print}
' | gzip > acrocomia_chr1.fna.gz

# -----------------------------------------
O que o script faz:
# 1. zcat: Descomprime o genoma original (.gz) em tempo real.
# 2. awk verifica cada linha. Se for um cabeçalho (/^>/) e contiver "chromosome 1,", ele liga o "interruptor" (extrair = 1) e imprime o cabeçalho na tela.
# 3. Se for qualquer outro cabeçalho (chromosome 2, chromosome 10, scaffold, etc), ele desliga o interruptor (extrair = 0).
# 4. Enquanto o interruptor estiver ligado (extrair == 1), ele imprime as linhas de DNA.
# 5. gzip: Recomprime apenas as linhas do cromossomo 1 e salva no novo arquivo.
# -----------------------------------------

```

Agora vamos explorar o genoma de referência para nos familiarizarmos com ele. Lembre-se de que o genoma é um arquivo fasta. 

❓Como esperam que seja o genoma?
❓Quantas e quais sequências estão presentes no arquivo fasta que usaremos como referência?

```bash
zcat acrocomia_ref.fna.gz | head
zcat acrocomia_ref.fna.gz | grep "^>"
```
Uma boa prática é modificar o nome dos cromossomos para que seja muito mais simples ler nos próximos programas. Para isso, podemos usar o seguinte script para modificar e deixar o genoma pronto para fazer o índice. 

```bash
zcat acrocomia_ref.fna.gz | awk '
BEGIN {chr=1; un=1}
/^>/ {
    if ($0 ~ /chromosome/) {
        print ">CHR" chr
        chr++
    } else {
        print ">UNRE" un
        un++
    }
    next
}
{print}
' | gzip > acrocomia_ref_rename.fna.gz

# -----------------------------------------
O que o script faz:
# 1. zcat: Descomprime o genoma original (.gz) em tempo real para leitura.
# 2. Detecta cabeçalhos FASTA (/^>/).
# 3. Modifica apenas linhas que começam com >.
# 4. Se o cabeçalho contiver "chromosome" ($0 ~ /chromosome/), então atribui:
# >CHR1
# >CHR2
# >CHR3
# usando o contador chr.
# 5. Se não, presume que é um contíguo não resolvido (scaffold) e atribui:
# >UNRE1
# >UNRE2
# usando o contador un.
# 6. gzip: Recomprime o resultado instantaneamente e salva o novo arquivo .gz.
# -----------------------------------------
```

# 2. Indexar a referência usando o BWA
Com o genoma de referência pronto e conhecido, podemos gerar o índice. Lembre-se de que o índice é a estratégia para pesquisar rapidamente dentro do genoma e otimizar os métodos de mapeamento. 

Este é um processo muito simples, utilizando o programa [bwa](https://github.com/lh3/bwa).

```bash
bwa index acrocomia_ref_rename.fna.gz

Se o programa estiver dentro de uma pasta:
./bwa index ../../data/reference/acrocomia_ref_rename.fna.gz
```

🔍 O que observar: 
Use o comando `ls -lh`. Você verá arquivos com extensões `.bwt`, `.ann`, `.pac`, `.sa` e `.amb`. 

Atenção: Se você mudar o arquivo .fasta de pasta, precisa levar todos esses arquivos de índice junto, senão o alinhamento falhará!


# 3. Mapeamento do genoma de referência utilizando o BWA-MEM
Agora vamos correr o [bwa mem](https://github.com/lh3/bwa) para mapear as reads da macaúba a esta referência. O alinhamento é o processo de encontrar a posição exata de cada read no genoma da Macaúba. Como nossas sequências são Single-End, o comando é direto. 

### A. Criando a estrutura de pastas
Primeiro, vamos organizar onde os resultados (arquivos SAM) serão salvos.

```bash
# 1. Criar pasta para o alinhamento e entrar nela
cd ~/workshop_bioinfo/data/
mkdir -p alignment
cd alignment
```

### B. Rodando o Alinhamento
Vamos testar o comando para a primeira amostra para compreender o que está sendo feito. 

```bash
Comando para alinhar UMA amostra
bwa mem -t 4 \
    -R "@RG\tID:Acro_01\tSM:Acro_01\tPL:ILLUMINA" \
    ../../data/reference/acrocomia_ref_renamed.fna.gz \
    ../../data/trimmed/AM_04_trimmed.fq.gz \
    > ../../data/alignment/AM_04.sam

# -----------------------------------------
Os parámetros que estamos usando:
#bwa mem: Algoritmo de alinhamento BWA otimizado para leituras de 70–1000 bp.
#-t 4: Utiliza 4 threads (núcleos de CPU). Se o seu computador tiver mais núcleos, você pode aumentar (-t 8, -t 16)
#-R "@RG\tID:Acro_01\tSM:Acro_01\tPL:ILLUMINA": Grupo de leitura. Isso adiciona metadados ao BAM. Ferramentas como o GATK exigem isso.
#../../reference/acrocomia_ref_renamed.fasta: Seu genoma renomeado e indexado. Lembre-se de que todos os arquivos da indexação devem estar juntos. 
#../trimmed/Acrocomia_pop1_.fastq: As leituras limpas e sem adaptadores. 
#> Acro_01.sam: O símbolo de "maior que" redireciona a saída do programa. Salva o alinhamento no formato SAM.
# -----------------------------------------
```

Use o less para se mover dentro do arquivo SAM.

❓Como está o arquivo? O que significa da parte dele? 

Agora que validamos o comando, vamos rodar para todos os arquivos que limpamos anteriormente. da forma como rodaríamos na vida real. 

⚠️ **Atenção:**
O arquivo SAM é muito pesado e ocupa muito espaço no computador, por isso normalmente não o geramos. Na prática, usamos `view` do [samtools](https://www.htslib.org/doc/samtools-depth.html) para converter automaticamente o SAM em BAM. 

```bash
# Loop para alinear y convertir todas las muestras Single-End
for file in ../trimmed/*_trimmed.fq.gz; do
    # Extraer solo el nombre de la muestra (ej: Acro_01)
    SAMPLE=$(basename $file _trimmed.fq.gz)
    
    echo "Alineando y convirtiendo la muestra: $SAMPLE"

    bwa mem -t 4 \
        -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" \
        ../reference/acrocomia_chr1.fna.gz \
        $file | \
 samtools view -b -> ${SAMPLE}_nofilers.bam
done
```

Vamos agora usar a mesma ferramenta de `view` para visualizar o resultado do bam.

```bash
samtools view ${SAMPLE}_nofilers.bam | head
```

❓O que é cada parte do arquivo bam?

# 4. Visualizações e filtros pós-alinhamento
Agora já temos nossas sequências alinhadas no formato BAM. No entanto, antes de poder usá-las para chamar os SNPs, precisamos verificar a qualidade do mapeamento e realizar os filtros necessários.

Vamos aplicar um **Filtro de Qualidade de Mapeamento (MAPQ)** e, em seguida, organizaremos as informações.

### A. Filtrar e ordenar (Samtools)
Os programas posteriores precisam que as leituras estejam ordenadas da esquerda para a direita ao longo dos cromossomos (coordenadas genômicas). Podemos filtrar as leituras ruins e ordená-las ao mesmo tempo usando nosso velho amigo, o Pipe (`|`).

```bash
# Filtrar reads não mapeados (-F 4), reter reads com qualidade >= 20 e ordenar por coordenadas
samtools view -F 4 -q 20 -b Acro_01.bam |\
samtools sort -o Acro_01_filter.sorted.bam

# Indexar o arquivo ordenado (Gera o índice posicional .bai)
samtools index Acro_01.sorted.bam

#Para todas as amostras
for file in ./*_nofilers.bam ; do
    # Extraer solo el nombre de la muestra (ej: Acro_01)
    SAMPLE=$(basename $file _nofilers.bam)
    
    echo "Filtrando: $SAMPLE"

    samtools view -F 4 -q 20 -b $file | \
    samtools sort -o $SAMPLE.sorted.bam 
done
```

### B. Mitigação de Vieses Técnicos: Duplicatas de PCR (Picard) - NÃO APLICA PARA GBS. GBS SINGLE-END NÃO PRECISA RODAR. 
Durante a confecção das bibliotecas genômicas, a etapa de amplificação por PCR pode introduzir um viés estocástico, gerando cópias artificiais (clones) de um mesmo fragmento de DNA. Se estas duplicatas não forem expurgadas, o algoritmo de Variant Calling interpretará erroneamente essas cópias como evidência de alta profundidade de cobertura para um alelo específico, distorcendo severamente as estimativas de frequências alélicas da população. No caso de GBS Single-end, não precissa rodar. Fica só de referência para outro tipo de bibliotecas. 

Usaremos a ferramenta `MarkDuplicates` de [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) para mitigar este artefato.

```bash
# Executar a remoção de duplicatas com Picard MarkDuplicates
picard MarkDuplicates \
    I=Acro_01.sorted.bam \
    O=Acro_01.sorted.dedup.bam \
    M=Acro_01_dup_metrics.txt \
    REMOVE_DUPLICATES=true 

# A geração de um novo arquivo BAM exige uma nova indexação
samtools index Acro_01.sorted.dedup.bam

# -----------------------------------------
Os parámetros que estamos usando:
#I=Input
#O=Output
#M= Metrics. Gera um relatório quantitativo contendo a proporção exata de sequências que consistiam em artefatos de PCR.
#REMOVE_DUPLICATES=true. Fala para o programa deletar fisicamente as duplicatas dos dados, não apenas marcá-las.
# -----------------------------------------
```

Para todas as amostras
```Bash
for file in *.sorted.bam; do
    # Pegar o nome da amostra (ej: Acro_01)
    SAMPLE=$(basename "$file" .sorted.bam)

    echo "Removiendo duplicados para la muestra: $SAMPLE"

    # Executar a remoção de duplicatas com Picard MarkDuplicates
    picard MarkDuplicates \
        I="$file" \
        O="${SAMPLE}.sorted.dedup.bam" \
        M="${SAMPLE}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true

    echo "Indexando ${SAMPLE}.sorted.dedup.bam..."
    # A geração de um novo arquivo BAM exige uma nova indexação
    samtools index "${SAMPLE}.sorted.dedup.bam"

done
```


Podemos gerar automaticamente um resumo da eliminação de duplicatas em nossos arquivos.
```bash
echo -e "Amostra\tPares_Avaliados\tTaxa_Duplicacao" > resumo_duplicatas.txt

for file in *_dup_metrics.txt; do
    # Extrair o nome base da amostra
    SAMPLE=$(basename "$file" _dup_metrics.txt)
    
    # Extrair a linha de dados (logo abaixo do cabeçalho PERCENT_DUPLICATION)
    DADOS=$(grep -A 1 "PERCENT_DUPLICATION" "$file" | tail -n 1)
    
    # Pegar as colunas específicas (3 = Pares, 9 = Porcentagem de Duplicação)
    PARES=$(echo "$DADOS" | awk '{print $3}')
    TAXA=$(echo "$DADOS" | awk '{print $9}')
    
    # Adicionar os dados ao nosso arquivo de resumo
    echo -e "${SAMPLE}\t${PARES}\t${TAXA}" >> resumo_duplicatas.txt
done
```



# 5. Check as estatísticas básicas do bam
O último passo é revisar as estatísticas finais do mapeamento para verificar se tudo está correto e se aproveitamos bem as sequências. Este é um passo vital, pois nos permite considerar se precisamos ajustar os parâmetros de mapeamento ou se podemos seguir em frente. 

Para esta fase, utilizaremos duas ferramentas do samtools: 
1. `flagstat` para gerar estatísticas rápidas de alinhamentos do arquivo bam.
2. `depth` para calcular a cobertura do nosso arquivo bam.

```bash
# 1. Navegar para o diretório com os BAM finais
cd ~/workshop_bioinfo/data/processed/alignment

# 2. Criar pasta para os stats
mkdir stats

# 3. Criar arquivo de resumo da profundidade média
printf "Sample\tTotal_reads\tMapped_reads\tMapped_percent\tDuplicates\tProperly_paired\n" > stats/flagstat_summary.tsv

# 4. Loop para todas as amostras
for bam in *.sorted.dedup.bam; do

    SAMPLE=$(basename "$bam" .sorted.dedup.bam)
    echo "Processando $SAMPLE"

    # generar flagstat
    samtools flagstat "$bam" > stats/${SAMPLE}_flagstat.txt

    FILE=stats/${SAMPLE}_flagstat.txt

    # extraer métricas
    TOTAL=$(awk '/in total/ {print $1}' "$FILE")
    MAPPED=$(awk '/ mapped \(/ && !/primary mapped/ {print $1}' "$FILE")
    MAP_PERC=$(awk -F'[()%]' '/ mapped \(/ && !/primary mapped/ {print $2}' "$FILE")
    DUP=$(awk '/duplicates$/ {print $1}' "$FILE" | head -n1)
    PROPER=$(awk '/properly paired/ {print $1}' "$FILE")

    # agregar a la tabla
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$TOTAL" "$MAPPED" "$MAP_PERC" "$DUP" "$PROPER" \
    >> stats/flagstat_summary.tsv

done

# -----------------------------------------
Os parâmetros que estamos usando:

# printf = escreve o cabeçalho da tabela que irá resumir as métricas de alinhamento obtidas para cada amostra. 
# O cabeçalho contém as seguintes colunas:
# - Sample — identificação da amostra analisada.
# - Total_reads — número total de leituras presentes no arquivo BAM.
# - Mapped_reads — número de leituras que foram alinhadas ao genoma de referência.
# - Mapped_percent — porcentagem de leituras que foram alinhadas com sucesso.
# - Duplicates — número de leituras marcadas como duplicadas.
# - Properly_paired — número de leituras corretamente pareadas durante o alinhamento.

# O símbolo > indica redirecionamento de saída, o que significa que o conteúdo será escrito no arquivo especificado. 
# Caso o arquivo já exista, ele será sobrescrito, garantindo que um novo resumo seja gerado a cada execução do script.

# samtools flagstat = calcula um conjunto de estatísticas fundamentais do alinhamento presentes em um arquivo BAM. 
# Entre essas estatísticas estão o número total de leituras, leituras alinhadas, duplicadas e leituras corretamente pareadas.

# awk = é utilizado para extrair valores específicos do arquivo de saída gerado por samtools flagstat, permitindo selecionar apenas
# as métricas de interesse que serão adicionadas à tabela final de resumo.
# -----------------------------------------

```
Agora abra os arquivos gerados e veja as estatísticas gerais. 

❓Como foi o nosso mapeamento?

---
# Mapeando no dia a dia
Como vimos, durante esse processo são gerados muitos arquivos intermediários que não usamos para a chamada de SNPs. Isso ocupa espaço em nosso computador/servidor. 

Para mitigar esse problema, na vida real usamos um único script que faz todo o mapeamento, conversão e filtragem dos arquivos e retorna no final apenas um arquivo .bam final pronto para a chamada de SNPs. 

Quando estamos começando com a análise de dados, é recomendável seguir o processo passo a passo para entender o que está sendo feito e ter certeza do processo. Quando já estivermos familiarizados com os programas e parâmetros, podemos otimizar essa etapa. 

💡 Script "Tudo-em-Um"
```bash
# 1. Indexar el genoma de referencia (solo una vez)
echo "Indexando el genoma de referencia..."
bwa index ../reference/acrocomia_chr1.fna.gz


# 2. Loop de procesamiento para todas las muestras
for file in ../trimmed/*_trimmed.fq.gz; do

    # Extraer nombre de muestra
    SAMPLE=$(basename "$file" _trimmed.fq.gz)

    echo "Procesando muestra: $SAMPLE"

    # PASO A: Alineamiento + filtrado + ordenamiento
    bwa mem -t 4 \
        -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" \
        ../reference/acrocomia_chr1.fna.gz \
        "$file" | \
    samtools view -F 4 -q 20 -b - | \
    samtools sort -o ${SAMPLE}.tmp.sorted.bam

    # PASO B: Remoción de duplicados
    picard MarkDuplicates \
        I=${SAMPLE}.tmp.sorted.bam \
        O=${SAMPLE}.sorted.dedup.bam \
        M=${SAMPLE}_dup_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT

    # PASO C: Indexar BAM final
    samtools index ${SAMPLE}.sorted.dedup.bam

    # PASO D: eliminar archivo temporal
    rm ${SAMPLE}.tmp.sorted.bam

    echo "Muestra $SAMPLE finalizada."

done
```

