# 💻 Prática Dia 3: Chamada de SNPs

Nesta prática, aprenderemos duas maneiras diferentes de chamar variantes, com foco na chamada de SNPs. Faremos isso considerando dois cenários diferentes: com e sem genoma de referência. 

---
## 0. Lembrando os tipos de variantes 
As variantes ou variaçoes genômicas são as posições em que as reads não são idênticas à referência. Assim identificaremos variações genética entre indivíduos de uma espécie ou entre espécies, o que fornecerá informações importantes sobre diversidade genética, saúde, evolução e muito mais.

Essas variantes podem ocorrer em uma única base ou em cromossomos inteiros. Entre os tipos mais importantes estão os SNPs, os indels e as variantes estruturais (SVs).

1. **SNPs (Single Nucleotide Polymorphisms):** É a substitução de uma única base no DNA por outra base. É o tipo de variante mais simple e frequente no genoma. Elas são marcadores altamente informativos, ya que podem ocorrer em regiões codificantes ou não codificantes, influenciando uma série de processos, genes e proteínas.
2. **INDELs (Inserções e Deleções):** Indels são variantes em que ocorrem inserções ou deleções de pequenas sequências de nucleotídeos em comparação com a referência. Elas geralmente envolvem poucos pares de bases (<50pb), mas podem causar um impacto funcional significativo — especialmente em regiões codificantes de genes.
3. **Variantes Estruturais (SVs):** São alterações genômicas em larga escala, geralmente envolvendo segmentos com dezenas a milhões de bases. Ao contrário dos anteriores, as SVs podem modificar significativamente a arquitetura genômica. Entre os principais tipos de SVs estão:
    + _Deleções_: remoção de segmentos inteiros do DNA
    + _Duplicações_: cópia adicional de uma região
    + _Inversões_: quando um trecho do DNA é invertido
    + _Translocações_: troca de segmentos entre cromossomos diferentes

![Types-of-Genetic-Variation-Single-nucleotide-variants-SNVs-and-indels-are-changes-that_W640](https://github.com/user-attachments/assets/a8a404bb-07a4-432e-bd0d-2a869c573161)

Nós nos concentraremos nos SNPs

## 1. Chamada com genoma de referência 
Com a redução dos custos de sequenciamento, cada vez é mais comum ver genomas de referência para muitas espécies. Está se tornando rotineiro em muitos projetos gerar o próprio genoma de referência e usá-lo para estudos populacionais posteriores. Por isso, começamos hoje aprendendo como fazer a chamada de variantes (SNPs) usando as sequências mapeadas contra o genoma de referência. 

### Programas para chamada de SNPs
Existem diversas ferramentas bioinformáticas desenvolvidas para identificar variantes genômicas a partir de dados de sequenciamento. Cada ferramenta tem suas abordagens, pontos fortes e requisitos, e a escolha ideal depende do tipo de projeto, da cobertura dos dados, da qualidade do genoma de referência e dos objetivos da análise.

Entre as ferramentas mais utilizadas, destacam-se:
1. **[GATK (Genome Analysis Toolkit)](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller):** é considerado o padrão-ouro na genética humana e de populações. Ele utiliza um algoritmo sofisticado chamado HaplotypeCaller, que não olha apenas para uma base isolada, mas reconstrói localmente as sequências (haplótipos) para determinar com precisão se um SNP ou um pequeno InDel é biológico ou um artefato de alinhamento. No entanto, exige etapas complexas de pré-processamento (como realinhamento e recalibração de base) e depende de genomas de alta qualidade, o que pode ser um desafio em espécies não-modelo com genomas menos completos.
2. **[FreeBayes](https://github.com/freebayes/freebayes):** é uma alternativa similar a GATK, ja que é baseado em haplótipos. É frequentemente a ferramenta de escolha para genomas complexos devido à sua flexibilidade com ploidias variadas e sua capacidade nativa de analisar múltiplos indivíduos simultaneamente, utilizando o poder da probabilidade bayesiana para inferir a presença de polimorfismos na população como um todo. No entanto, ele também depende de genomas de alta qualidade. 
3. **[bcftools](https://samtools.github.io/bcftools/howtos/variant-calling.html):** é uma ferramenta rápido e de fácil uso, o bcftools é altamente recomendado para projetos com espécies não-modelo, genomas ainda em construção ou muito fragmentados. Ele oferece um fluxo simples, com poucos requisitos e excelente compatibilidade com outras ferramentas como samtools e vcftools. 

Como nosso foco são principalmente espécies vegetais comerciais, a maioria das quais já possui um genoma de referência de qualidade, aprenderemos a usar o GATK como ferramenta principal. 

### A. Pré-processamento antes de chamar SNPs
O GATK é extremamente exigente quanto aos metadados da referência. Antes de executá-lo, precisamos gerar um dicionário de sequências `(.dict)` e um índice FASTA `(.fai)`. 

```bash
# 1. Navegar até a pasta da referência
cd ~/workshop_bioinfo/data/reference

# 2. Criar o índice FASTA com Samtools
gunzip acrocomia_chr1.fna.gz
samtools faidx acrocomia_chr1.fna

# 3. Criar o Dicionário de Sequências com Picard
picard CreateSequenceDictionary \
    R=acrocomia_chr1.fna \
    O=acrocomia_chr1.dict

```

💡*Opcional, mas recomendado*

Recalibra os arquivos bam em relação a um “conjunto verdadeiro” de SNPs

Os vcfs do conjunto verdadeiro são variantes nas quais temos alta confiança e tendem a ser gerados a partir de bibliotecas de alta cobertura sem PCR. Se tiver dispon[ivel, é recomendado fazer a recalibração. Na maioria das plantas não modelo, não temos esse tipo de dados, portanto, vamos pular essa etapa. 

### B. Chamada de variantes
A chamada de variantes não é um passo simples, pois o GATK faz a chamada de variantes por indivíduo. Posteriormente, é necessário consolidar em um único arquivo e chamar as variantes comuns entre todos os indivíduos. 

O primeiro passo é chamar as variantes para cada indivíduo. Nesta etapa, o GATK analisará cada arquivo .dedup.bam e gerará um arquivo .g.vcf.gz. O sufixo .g indica que o arquivo contém informações de confiança para todas as bases do genoma, não apenas para os polimorfismos. É um arquivo intermediário que precisamos usar antes de criar nosso arquivo VCF final.

```bash
#Instalar GATK
conda install bioconda::gatk4

# Criar diretório para os gVCFs
mkdir -p ~/workshop_bioinfo/results/gatk_gvcfs
cd ~/workshop_bioinfo/data/alignment

# Loop para processar todas as amostras de Acrocomia simultaneamente
for bam in *.sorted.dedup.bam; do
    SAMPLE=$(basename $bam .sorted.dedup.bam)
    echo "Gerando gVCF para a amostra: $SAMPLE"

    gatk HaplotypeCaller \
        -R ~/workshop_bioinfo/data/reference/acrocomia_chr1.fna \
        -I $bam \
        -O ~/workshop_bioinfo/results/gatk_gvcfs/${SAMPLE}.g.vcf.gz \
        -native-pair-hmm-threads 4 \
        -ERC GVCF
done

Esta é a parte mais demorada de todo o processo. É hora de tomar um café ☕. 
```
💡 Verifique seu arquivo `gvcf` para garantir que ele tenha um arquivo de índice `.idx`. Se o `haplotypecaller` travar, ele produzirá um arquivo `gvcf` truncado que acabará travando a etapa `genotypegvcf`. 

O próximo passo é importar nossos arquivos `gvcf` para um arquivo `genomicsDB`. E uma representação compactada em banco de dados de todos os dados lidos em nossas amostras. O arquivo `GenomicsDB` contém todas as informações dos seus arquivos `GVCF`, mas não pode ser adicionado e não pode ser transformado de volta em um `gvcf`. Isso significa que, se você obtiver mais amostras, não poderá simplesmente adicioná-las ao seu arquivo `genomicdDB`, terá que voltar aos arquivos `gvcf`. 

💡 A ferrameta `GenomicsDBImport` substituiu o antigo `CombineGVCFs` porque utiliza um banco de dados estruturado, sendo infinitamente mais rápido e escalável para genomas grandes e populações numerosas. A grande particularidade computacional do `GenomicsDBImport` é que ele exige a definição de intervalos (cromossomos ou scaffolds) para particionar o trabalho, e utiliza um arquivo de mapeamento de amostras (Sample Map). Primeiro, devemos gerar estes dois arquivos

O algoritmo exige um arquivo tabular indicando o nome biológico da amostra e o caminho físico do seu respectivo gVCF (Sample Map).

```bash
# Navegar até o diretório dos gVCFs
cd ~/workshop_bioinfo/results/gatk_gvcfs

# Gerar o Sample Map automaticamente com um loop
> sample_map.txt
for gvcf in *.g.vcf.gz; do
    SAMPLE=$(basename "$gvcf" .g.vcf.gz)
    echo -e "${SAMPLE}\t${gvcf}" >> sample_map.txt
done
```
A segunda particularidade é a necessidade estrita de paralelização espacial: ele precisa saber exatamente em quais regiões (intervalos) do genoma da *Acrocomia* ele deve atuar. Como nosso genoma está nos cromossomos, podemos indicar uma lista com o nome dos cromossomos. Isso a partir do nome do índice `(.fai)`. 

```bash
# Extrair a primeira coluna (nomes das sequências) do arquivo .fai para criar a lista
awk '{print $1}' ~/workshop_bioinfo/data/reference/acrocomia_chr1.fna.fai > intervalos_acrocomia.list

# Visualizar os primeiros intervalos gerados
head intervalos_acrocomia.list
```

Em seguida, chamamos `GenomicsDBImport` para criar o banco de dados

```bash
# Importar os dados para o GenomicsDB utilizando a lista de intervalos
gatk GenomicsDBImport \
    --sample-name-map sample_map.txt \
    --genomicsdb-workspace-path banco_dados_acrocomia \
    -L intervalos_acrocomia.list

# -----------------------------------------
# sample-name-map: é o nome de saída da biblioteca que estamos construíndo
#-L: Também pode ser usado para um único cromossomo. Substituímos a lista pelo nome da região de interesse. 
# -----------------------------------------
```

Com o `genomicsDB` criado, estamos finalmente prontos para identificar variantes e gerar um `vcf`.
```bash
# Realizar a chamada de variantes lendo diretamente do banco de dados recém-criado
gatk GenotypeGVCFs \
    -R ~/workshop_bioinfo/data/reference/acrocomia_chr1.fna \
    -V gendb://banco_dados_acrocomia \
    -O ./acrocomia_populacao_bruta.vcf.gz
```
❓Vamos fazer um less neste arquivo. Como ele fica? Que variantes você consegue ver? 


O arquivo gerado ao final deste processo `(acrocomia_populacao_bruta.vcf.gz)` é categorizado como "bruto" (unfiltered). Ele contém todos os candidatos a polimorfismos. A próxima e mandatória etapa analítica é a filtragem criteriosa (Hard Filtering).

Como estamos interessados apenas em SNPs, o primeiro passo é selecionar apenas essas marcas. 
```bash
gatk SelectVariants \
-R ~/workshop_bioinfo/data/reference/acrocomia_chr1.fna \
-V acrocomia_populacao_bruta.vcf.gz \
--select-type SNP \
-O raw_snps.vcf.gz
```
 
### C. Filtragem de SNPs
Como podem ver, o arquivo final é chamado `raw_snps.vcf.gz`. Isso porque é um arquivo final sem filtros de qualidade nem nada. O GATK tem a função `VariantFiltration` que permite fazer esses filtros. 

```bash
gatk VariantFiltration \
-R acrocomia_chr1.fna \
-V raw_snps.vcf.gz \
-O filtered_snps.vcf.gz \
--filter-name "QD_filter" --filter-expression "QD < 2.0" \
--filter-name "FS_filter" --filter-expression "FS > 60.0" \
--filter-name "MQ_filter" --filter-expression "MQ < 40.0"
```

No entanto, existem outras duas ferramentas mais comuns e utilizadas para realizar esta etapa: `bcftools` e `VCFtools`. Iremos utilizá-las na próxima prática para realizar os filtros no nosso VCF e deixar o arquivo final pronto para as nossas análises. 


## 2. Chamada com genoma de referência (Alternativa)
Como vimos nas seções anteriores, o GATK é uma ferramenta computacionalmente robusta e de altíssima precisão metodológica.

Contudo, é fundamental ter em mente um viés histórico: o GATK foi originalmente desenvolvido e calibrado para analisar dados genômicos humanos. Consequentemente, seus algoritmos probabilísticos (como o *HaplotypeCaller* e o *VQSR*) pressupõem implicitamente a existência de um genoma de referência extremamente bem curado e dados de sequenciamento com a mesma qualidade. 

Na botânica e melhoramento de plantas não sempre é possível. Para pesquisadores que trabalham com espécies de plantas nativas ou menos estudas, as condições ideais exigidas pelo GATK nem sempre podem ser atendidas.

Diante desses desafios inerentes à genômica de espécies não-modelo, destacam-se alternativas de altíssima confiabilidade e flexibilidade na comunidade científica: o **FreeBayes**, **BCFtools** ou **STACK**. 


### 0. Preparação da Lista de Amostras
O FreeBayes opera de maneira otimizada lendo os dados de todos os indivíduos da sua coorte de uma só vez. Isso confere ao algoritmo um enorme poder estatístico: se uma mutação rara aparece com alta qualidade em um indivíduo da população, o programa usará essa evidência para refinar a busca da mesma mutação nos demais indivíduos, mesmo naqueles com menor cobertura de sequenciamento.

Vamos primero instalar o programa

```bash
conda install bioconda::freebayes
```

Antes de poder rodar o algoritmo, ele exige um arquivo de texto simples listando o caminho de todos os arquivos BAM que serão analisados. Vamos gerar essa lista automaticamente.

```bash
# 1. Navegue para a pasta onde estão os BAMs processados e deduplicados
mkdir ~/workshop_bioinfo/results/freebayes
cd ~/workshop_bioinfo/data/alignment

# 2. Indexar os arquivos bam
#Nós já fizemos isso, então não é necessário repetir. Mas é importante lembrar que, antes de tudo, precisamos ter os arquivos indexados para que os programas possam funcionar.

for bam in *.bam; do
    samtools index $bam
done

# 3. Crie uma lista contendo os nomes de todos os arquivos .dedup.bam
ls *.sorted.dedup.bam > lista_amostras_acrocomia.txt

# 4. Verifique se a lista foi criada corretamente
cat lista_amostras_acrocomia.txt

```

Agora estamos prontos para rodar freebayes. Vemos que é muito mais fácil e simples do que o GATK. 

```bash
freebayes \
-f ~/workshop_bioinfo/data/reference/acrocomia_chr1.fna \
-L lista_amostras_acrocomia.txt \
> ~/workshop_bioinfo/results//freebayes_raw.vcf

# -----------------------------------------
Parâmetros utilizados: 
# f: genoma de referência (indexado)
#-L: lista de BAMs
#--standard-filters
#>: arquivo VCF de saída
# -----------------------------------------
```

Agora vamos inspecionar o VCF bruto sem filtros. 

```bash
less -S ~/workshop_bioinfo/results/freebayes_raw.vcf
#O parâmetro -S evita a quebra de linha, permitindo rolar a tela lateralmente com as setas do teclado. 
```
❓Como o resultado do GATK difere do FreeBayes?



## 3. Chamada sem genoma de referência (_de novo_)
Até agora as ferramentas que exploramos (como GATK e FreeBayes) compartilham uma premissa metodológica fundamental: a dependência estrita de um genoma de referência. Esses algoritmos operam mapeando e comparando cada leitura de sequenciamento contra uma sequência genômica preestabelecida para, então, identificar os SNPs.Porém, na prática da biologia da conservação e da ecologia evolutiva, nem sempre temos o genoma disponível. 

Para isso foram desenvolvidos algoritmos capazes de realizar a montagem e a chamada de variantes ***de novo*** (sem qualquer referência genômica prévia). Historicamente, a vasta maioria dessas ferramentas foi concebida para processar bibliotecas de Representação Reduzida do Genoma, como o GBS (*Genotyping-by-Sequencing*) e o RADseq. Existe uma clara tendência na genômica de transição das abordagens de GBS para o Sequenciamento de Genoma Completo (WGS - *Whole Genome Sequencing*), mas no conteto latino-americano, o GBS continua sendo uma realidade metodológica e econômica vital para a viabilização de estudos populacionais.

Dentre os *pipelines* bioinformáticos mais conhecidos para a inferência *de novo* de variantes a partir de bibliotecas de representação reduzida, destacam-se:

1. **[STACKS](https://catchenlab.life.illinois.edu/stacks/)**: Uma das ferramentas mais robustas e amplamente adotadas, baseada na construção de catálogos de *loci* consensuais a partir do agrupamento de sequências idênticas.
2. **[ipyrad](https://ipyrad.readthedocs.io/en/master/index.html):** Uma evolução do antigo *PyRAD*, com excelente capacidade para lidar com InDels e otimizado para análises filogenéticas e de introgressão em dados RADseq.
3. **[TASSEL](https://tassel.bitbucket.io):** Um software clássico, extensivamente utilizado no melhoramento genético de plantas (especialmente para milho e soja), com um *pipeline* muito consolidado para dados de GBS (o *TASSEL-GBS*).

No nosso caso, vamos usar o `denovo_map.pl` [STACKS](https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php). Ele é uma ferramenta muito utilizada para estudos de genética populacional em espécies não modelo. Tem a vantagem de já incluir todas as etapas até o VCF final. 


### 0. Visão geral do `denovo_map.pl`

Quando executamos o `denovo_map.pl`, ele executa autonomamente os seguintes módulos, na exata ordem:
1. **`ustacks` (Unique Stacks - Construção Intra-indivíduo):** O algoritmo processa os dados limpos de cada amostra separadamente. Ele agrupa sequências idênticas em pilhas (*stacks*) e, em seguida, funde pilhas com pequenas diferenças (possíveis alelos) para formar *loci* primários dentro de um único indivíduo.
2. **`cstacks` (Catalog Stacks - O Genoma de Referência Virtual):** O programa seleciona os indivíduos representativos da população e mescla os seus *loci* primários. O produto desta etapa é um "Catálogo", que funcionará como um genoma de referência sintético para o restante da análise.
3. **`sstacks` (Sets of Stacks - Ancoragem):** Cada *locus* identificado individualmente no passo 1 é agora alinhado contra o Catálogo  gerado no passo 2. Isso garante que o Locus 150 do Indivíduo A seja comparado exatamente com o Locus 150 do Indivíduo B.
4. **`tsv2bam` e `gstacks` (Genotipagem e Faseamento):** Estes módulos reorganizam os dados alinhados, otimizam a chamada de variantes (SNP Calling) utilizando modelos probabilísticos avançados e realizam o faseamento de SNPs adjacentes no mesmo *locus*, construindo micro-haplótipos.
5. **`populations` (Genética de Populações):** O módulo final. Ele lê o catálogo genotipado, aplica filtros populacionais (como frequência alélica e presença mínima nas populações) e calcula estatísticas-chave ($F_{ST}$, $\pi$, heterozigosidade). Por fim, ele exporta os dados no formato padrão VCF.

<img width="1165" height="895" alt="stacks_pipeline" src="https://github.com/user-attachments/assets/cf71b50b-6c96-4ee2-943f-534b4f3c4075" />


### 1. Criar os arquivos necessários para executar o programa
Antes de executarmos o *pipeline* integral do STACKS, precisamos do o arquivo de metadados populacionais, comumente referido como *PopMap* (`popmap.txt`). A arquitetura do arquivo é estritamente tabular. Ele deve conter exatamente duas colunas, separadas **por uma tabulação (TAB)**:
1. **ID da Amostra:** O nome exato do indivíduo (deve ser estritamente idêntico ao nome do arquivo `.fq` ou `.fastq.gz` correspondente, mas **sem** a extensão).
2. **Identificador da População:** Uma *string* (texto sem espaços) indicando a origem ou agrupamento do indivíduo.

Exemplo:
```bash
Acro_01    Cerrado
Acro_02    Cerrado
Acro_03    Cerrado
Acro_04    Mata_Atlantica
Acro_05    Mata_Atlantica
```

A metodologia mais eficiente e segura consiste em organizar seus metadados previamente em um software de planilhas (como Microsoft Excel ou Google Sheets). Ao copiar duas colunas de uma planilha e colá-las no terminal, o sistema operacional automaticamente insere o caractere de tabulação (TAB) entre elas.

```bash
# 1. Navegue até a pasta de dados onde o STACKS buscará as informações
mkdir -p ~/workshop_bioinfo/results/stacks_denovo

# 2. Abra o editor de texto nativo do Linux (nano) criando um novo arquivo
nano popmap.txt

#Copiar e colar
AM_04_trimmed   amazonas
AM_05_trimmed   amazonas
G109_trimmed    guatemala
G111_trimmed    guatemala
rifania10_trimmed       sao_paulo
rifania8_trimmed        sao_paulo
luz_20m_trimmed minas
luz_25m_trimmed minas
veracruz-mex-L1p2_trimmed       mexico
veracruz-mex-L3p2_trimmed       mexico


#Para salvar o arquivo no nano, pressione Ctrl+O, confirme o nome apertando Enter, e saia do editor pressionando Ctrl+X.
```

### 2. Executando o Pipeline `denovo_map.pl`
A beleza do `denovo_map.pl` reside na sua capacidade de condensar toda essa complexidade em um único comando de terminal. Para executá-lo, precisaremos dos nossos arquivos FASTQ limpos e desmultiplexados, além dos arquivos gerados anteriormente. 

```bash
# 2. Executar o pipeline orquestrador De Novo
~/workshop_bioinfo/programs/stacks-2.68/scripts/denovo_map.pl \
-T 4 \
-m 3 \
-M 3 \
-n 3 \
-o . \
-O popmap.txt \
--samples ~/workshop_bioinfo/data/trimmed

# -----------------------------------------
Parâmetros utilizados: 
# -T (Threads):Alocação de núcleos de processamento para paralelizar a análise.
# -M: O número máximo de nucleotídeos divergentes permitidos entre duas sequências para que o programa as considere como sendo alelos do mesmo locus em um único indivíduo.
# n: O número máximo de nucleotídeos divergentes permitidos ao fundir loci de indivíduos diferentes para a construção do catálogo.
# -o: arquivo VCF de saída
#-b: o número que damos para o processo
# --p: a pasta onde os arquivos das amostras estão armazenados
# --O: arquivo descritivo das subpopulações
# -----------------------------------------
```

Dentro da pasta de saída você terá:
1. `populations.snps.vcf`
2. `populations.sumstats.tsv`
3. `catalog.fa.gz`

O arquivo principal para análises é:

`populations.snps.vcf`

Abra o arquivo usando o comando less

❓Como o resultado do GATK difere do FreeBayes?


💡 **Nota Metodológica:** 
> Neste exercício prático, optamos por executar o *pipeline* orquestrador completo (`denovo_map.pl`). Como vimos, ele percorre autonomamente todas as etapas computacionais e nos entrega diretamente o arquivo de resultados finais. Esta é, sem dúvida, a abordagem mais ágil e direta.
> 
> No entanto, é  possível executar o STACKS passo a passo. Isso significa invocar manualmente no terminal cada um dos módulos isolados (`ustacks`, depois `cstacks`, `sstacks`, `gstacks` e, finalmente, `populations`). 
> 
> Esta abordagem modular, embora mais trabalhosa, confere ao pesquisador um controle absoluto sobre a análise, permitindo refinar, testar e calibrar os parâmetros matemáticos específicos de cada fase de forma independente para otimizar a descoberta de variantes na sua espécie de estudo.
