# 💻 Prática Dia 3 (tarde): Filtros do arquivo VCF

A última etapa antes de iniciar as análises populacionais é filtrar o VCF final. Isso é feito para garantir que os SNPs sejam de alta qualidade e padronizados entre todas as amostras. 

---
## 0. Lembrando do arquivo VCF
O `Variant Call Format (VCF)` é o o formato universal que armazena toda a variabilidade genética e informação de variantes estructurais. Estruturalmente, um arquivo VCF divide-se em dois grandes blocos hierárquicos:
1. **Cabeçalho (Header)**: Um extenso bloco inicial de metadados, onde cada linha obrigatoriamente começa com ##. Ele atua como um dicionário, definindo de forma rigorosa as abreviações estatísticas, os filtros aplicados, os formatos dos dados e a versão das ferramentas bioinformáticas utilizadas na geração do arquivo.
2. **Matriz de Variantes (Data Table)**: Inicia-se com a linha de cabeçalho das colunas (marcada por um único #CHROM). A partir deste ponto, cada linha subsequente descreve detalhadamente um único polimorfismo (SNP ou InDel) encontrado na população.

O cabeçalho do VCF segue esse padrão:
| Campo | Descrição Acadêmica |
| :--- | :--- |
| **`#CHROM`** | O cromossomo ou *scaffold* de referência onde o polimorfismo está ancorado. |
| **`POS`** | A coordenada genômica exata (posição física em pares de bases) da mutação. |
| **`ID`** | Identificador alfanumérico da variante (ex: em humanos, o número dbSNP). Em espécies não-modelo sequenciadas *de novo*, este campo geralmente recebe um ponto (`.`). |
| **`REF`** | O Alelo de Referência (a base nucleotídica ancestral ou presente no genoma de referência). |
| **`ALT`** | O Alelo Alternativo (a base mutada inferida na nossa população amostrada). |
| **`QUAL`** | Escore de qualidade (*Phred-scaled*) da chamada da variante. Valores maiores indicam uma maior probabilidade estatística de que o polimorfismo seja biologicamente real e não um erro técnico. |
| **`FILTER`** | Status de aprovação na filtragem de qualidade (*Hard Filtering*). Variantes validadas recebem a etiqueta **`PASS`**. Variantes reprovadas exibem o nome do filtro que falhou (ex: `LowQual`, `StrandBias`). |
| **`INFO`** | Metadados e métricas globais consolidadas do *locus* (ex: profundidade total da coorte, frequência alélica da população). |
| **`FORMAT`** | Um dicionário (ex: `GT:AD:DP:GQ:PL`) que define a estrutura exata dos dados genotípicos contidos nas colunas individuais à direita. |
| **`<amostras>`** | Colunas terminais nomeadas com o ID de cada indivíduo (ex: `Acro_01`, `Acro_02`). Contêm o genótipo inferido (`0/0` homozigoto referência, `0/1` heterozigoto, `1/1` homozigoto alternativo) e suas métricas locais, seguindo a regra ditada na coluna `FORMAT`. |
<img width="1504" height="619" alt="vcf_format" src="https://github.com/user-attachments/assets/cbda7aa1-e532-4f98-a0ec-a1029d461076" />


Para manipular os arquivos VCF existem duas ferramentas muito utilizadas: 
1. **[BCFtools](https://www.htslib.org/doc/1.0/bcftools.html)**:  é a ferramenta de eleição para a manipulação técnica e estrutural da matriz de dados. Ele permete conversão de formatos (de VCF textual para BCF binário comprimido), concatenação e fusão de múltiplos cromossomos, indexação rápida e a aplicação de filtros técnicos complexos baseados nas anotações matemáticas da coluna `INFO` (como o *Mapping Quality* ou *Strand Bias*).
2. **[VCFtools](https://vcftools.github.io/man_latest.html)**: Ele é mais antigo, e tem alguns friltros não disponíveis no BCFtools, muitos com sentido biológico importante. Utilizamos o VCFtools para expurgar genótipos com excesso de dados faltantes (*Missing Data*) e remover alelos excessivamente raros (*Minor Allele Frequency* - MAF).

Vamos a instalar eles
```bash
conda install -c bioconda -c conda-forge bcftools openssl -y
conda install -c bioconda -c conda-forge vcftools -y
```


## 1. Filtrando com bcftools
O nosso objetivo é filtrar os VCF gerados com GATK e FreeBayes para reter apenas SNPs:
+ bialélicos
+ alta qualidade
+ MAF ≥ 0.05
+ presentes em ≥60% das amostras
+ remover SNPs com >50% de missing data
+ calcular estatísticas finais do VCF


Utilizaremos o `bcftools view` para reter estritamente os SNPs e forçar a exclusão de qualquer *locus* que possua mais de um alelo alternativo.

```bash
bcftools view \
-v snps \
-m2 -M2 \
-q 0.05:minor \
-i 'F_MISSING<=0.4' \
acrocomia_populacao_bruta.vcf.gz \
-Oz \
-o filtered_snps.vcf.gz

# -----------------------------------------
O que o script faz
# bcftools view:
# -v snps: mantém apenas SNPs
# -m2 -M2: mantém apenas variantes bialélicas (min e max)
# -q 0.05:minor: frequência mínima do alelo menor (MAF ≥ 0.05)
# -i 'F_MISSING<=0.4': máximo de 40% de missing (≥60% presentes)
# raw_variants.vcf.gz: arquivo de entrada
# -Oz: saída comprimida em VCF.gz
# -o filtered_snps.vcf.gz:  arquivo de saída
# -----------------------------------------
```


Agora que temos apenas os SNPs que nos interessam e que passaram pelos filtros de características dos SNPs, vamos filtrar por qualidade.  Para isso, utilizaremos o comando `bcftools filter`, aplicando uma expressão lógica combinada que exige simultaneamente uma alta probabilidade de acerto do sequenciador e um mínimo de leituras (reads) ancoradas naquela posição.

```bash
# Filtrar variantes mantendo apenas aquelas com Qualidade > 20 e Profundidade Total > 5

bcftools filter \
-i 'QUAL>30 && INFO/DP>5' \
-Oz \
-o acrocomia_snps_quality.vcf.gz \
filtered_snps.vcf.gz
```

Após aplicarmos os primeiros filtros temos que entender a taxa de retenção dos dados. Quantos polimorfismos tínhamos na população bruta e quantos sobreviveram depois da filtragem? 

Podemos descobrir isso usando a mesma ferramenta `bcftools view`

```bash
echo "SNPs no VCF original:"
bcftools view -H raw_variants.vcf.gz | wc -l

echo "SNPs após filtragem:"
bcftools view -H snp_boa_qualidade.vcf.gz | wc -l
```

❓ Quantos SNPs havia antes e depois? Como podemos alterar esse número? 


## 2. Filtrando com vcftool
Até o momento, filtramos nossa matriz focando na qualidade das mutações (os *loci*). Como vimos, dependendo de como a biblioteca foi construída e da tecnologia de sequenciamento, é comum que algumas amostras de DNA rendam bibliotecas de baixa qualidade, resultando em indivíduos com poucas leituras (baixa profundidade) e muitos dados faltantes.

Manter indivíduos ruins na análise pode distorcer severamente os cálculos de estrutura populacional e diversidade genética. Utilizaremos o **VCFtools** para gerar um diagnóstico clínico da coorte e removeremos qualquer indivíduo que apresente mais de 30% de dados faltantes.

### A Diagnóstico de Profundidade e Missing Data
Primeiro, vamos calcule as estatísticas de qualidade para cada indivíduo da nossa população.

```bash
# 1. Calcular a proporção de dados faltantes por indivíduo (Missing Data)
vcftools --gzvcf acrocomia_snps_quality.vcf.gz  --missing-indv --out stats_missing

# 2. Calcular a profundidade média de sequenciamento por indivíduo (Depth)
vcftools --gzvcf acrocomia_snps_quality.vcf.gz --depth --out stats_depth
```

Estes comandos não alteram o VCF. Eles geram arquivos de relatório em texto simples (`stats_acrocomia.imiss` e `stats_acrocomia.idepth`), que contêm as métricas exatas de cada planta sequenciada.

❓ Qual é a cobertura média das amostras? Temos alguma amostra com muitos dados ausentes? 


### Remover os indivíduos com muitos dados faltantes
Primeiro, precisamos gerar um arquivo .txt com a lista das pessoas que queremos excluir. 

```bash
# Ler o arquivo .imiss e extrair o nome (coluna 1) das amostras com mais de 30% de falha (coluna 5 > 0.3)
# O comando 'NR>1' ignora o cabeçalho da tabela
awk 'NR>1 && $5 > 0.3 {print $1}' stats_missing.imiss > amostras_para_remover.txt

# Conferir quantas e quais amostras reprovaram no controle de qualidade
echo "Amostras com >30% de missing data:"
cat amostras_para_remover.txt
```

Com a nossa lista de indivíduos de baixa qualidade gerada (`amostras_para_remover.txt`), instruímos o VCFtools a ler o nosso VCF, deletar inteiramente essas amostras e recalcular as frequências.

```bash
# Remover os indivíduos problemáticos e gerar o VCF final limpo
vcftools \
    --gzvcf acrocomia_snps_quality.vcf.gz \
    --remove amostras_para_remover.txt \
    --recode --recode-INFO-all \
    --out final_filtered_snps

#Verificar número final de SNPs
grep -vc "^#" final_filtered_snps.recode.vcf
```

Em seguida, é necessário fazer a indexação do vcf gerado.
```bash
bgzip final_filtered_snps.recode.vcf
bcftools index final_filtered_snps.recode.vcf.gz
```
❓O número final do SNPS mudou?


## 3. Estadísticas finais
Depois de aplicar todos os filtros (bialélicos, MAF, missing data, etc.), é importante verificar se o dataset final está consistente.

```bash
# 1. Gerar as estatísticas compreensivas do VCF final lapidado
bcftools stats final_filtered_snps.recode.vcf > estatisticas_finais_acrocomia.txt

# 2. Visualizar o sumário de números básicos (Summary Numbers)
grep "^SN" estatisticas_finais_acrocomia.txt
```


## 3. Filtrando os outros arquivos. 

Geramos três VCF diferentes, usando programas diferentes. Agora filtramos apenas um deles. Agora, usando os mesmos parâmetros e filtros, filtre os gerados pelo Freebayes e STACK. 


❓Que diferenças você encontra entre os três VCF finais? 


