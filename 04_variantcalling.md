# 💻 4: Chamada de SNPs

Nesta prática, vamos aprender a chamar variantes, com foco em SNPs, em dois cenários diferentes: **com** e **sem** genoma de referência. Para que tudo caiba no tempo da aula, vamos trabalhar só com o **cromossomo 5**.

---

## 0. Lembrando os tipos de variantes

As variantes genômicas são as posições em que os reads não são idênticos à referência. Elas revelam a variação genética entre indivíduos de uma espécie ou entre espécies, e são a base para estudar diversidade genética, estrutura populacional, adaptação e evolução.

Essas variantes podem envolver uma única base ou cromossomos inteiros. Os tipos mais importantes são:

1. **SNPs (*Single Nucleotide Polymorphisms*):** a substituição de uma única base por outra. É o tipo de variante mais simples e mais frequente no genoma. São marcadores muito informativos, que podem ocorrer em regiões codificantes ou não codificantes.
2. **InDels (inserções e deleções):** a inserção ou deleção de pequenas sequências em comparação com a referência, geralmente com menos de 50 pb. Podem ter um impacto funcional importante, especialmente em regiões codificantes, onde podem alterar a fase de leitura.
3. **Variantes estruturais (SVs):** alterações em larga escala, envolvendo de dezenas a milhões de bases, que modificam a arquitetura do genoma. Entre elas estão:
    + *Deleções*: perda de segmentos inteiros de DNA
    + *Duplicações*: cópias adicionais de uma região
    + *Inversões*: um trecho do DNA invertido
    + *Translocações*: troca de segmentos entre cromossomos diferentes

![Tipos de variação genética](https://github.com/user-attachments/assets/a8a404bb-07a4-432e-bd0d-2a869c573161)

Nesta prática, vamos nos concentrar nos SNPs.

---

## 1. Chamada com genoma de referência

Com a redução dos custos de sequenciamento, cada vez mais espécies têm um genoma de referência, e é cada vez mais comum que um projeto gere o seu próprio genoma para depois usá-lo em estudos populacionais. Por isso, começamos chamando os SNPs a partir dos reads que mapeamos contra a referência.

### Programas para chamada de SNPs

Existem diversas ferramentas para identificar variantes. Cada uma tem a sua abordagem, os seus pontos fortes e os seus requisitos, e a escolha depende do tipo de dado, da cobertura, da qualidade do genoma de referência e dos objetivos do projeto. Entre as mais usadas:

1. **[GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)**: é considerado o padrão-ouro na genética humana. O seu algoritmo principal, o HaplotypeCaller, não olha cada base isoladamente: ele reconstrói localmente as sequências (haplótipos) de cada região para decidir se um SNP ou InDel é real ou um artefato de alinhamento. Em compensação, o fluxo de trabalho tem várias etapas (gVCFs por indivíduo, banco de dados, genotipagem conjunta) e foi calibrado com dados humanos.
2. **[FreeBayes](https://github.com/freebayes/freebayes)**: também é baseado em haplótipos. Analisa todos os indivíduos ao mesmo tempo, usando inferência bayesiana para estimar os polimorfismos na população como um todo, e lida bem com diferentes ploidias.
3. **[bcftools](https://samtools.github.io/bcftools/howtos/variant-calling.html)**: é rápido, simples de usar e tem poucos requisitos. Analisa cada posição do genoma a partir da "pilha" de reads que a cobrem. Funciona muito bem com espécies não-modelo e genomas ainda em construção, e integra-se naturalmente com o samtools e o vcftools.

Nesta prática, vamos usar o **bcftools** como ferramenta principal e o **FreeBayes** como alternativa, para comparar os resultados.

### A. Preparação

Para chamar variantes, precisamos de duas coisas que já temos prontas da prática de mapeamento: a referência indexada e os arquivos BAM de cada amostra. Vamos conferir:

```bash
# A referência e os seus índices (precisamos do .fai)
ls -lh ~/referencia/

# Os alinhamentos finais: um .bam e um .bam.bai por amostra
ls -lh ~/alinhamento/*.bam*

# Criar a pasta para os resultados
mkdir -p ~/vcf

# Guardar o caminho da referência em uma variável
REF=~/referencia/chr5.fa
```

> [!NOTE]
> O índice `.fai` da referência, que criamos com o `samtools faidx` na prática de mapeamento, é obrigatório para a chamada de variantes. Se ele não estiver na pasta, crie de novo com `samtools faidx ~/referencia/chr5.fa`.

> [!IMPORTANT]
> Lembre-se: a variável `REF` só existe na sessão atual do terminal. Se você sair do servidor e entrar de novo, defina outra vez com `REF=~/referencia/chr5.fa`.

### B. Chamando as variantes com o bcftools

A chamada com o bcftools tem duas etapas, ligadas por um *pipe* (`|`):

1. **`bcftools mpileup`**: para cada posição do cromossomo, empilha (*pileup*) os reads de todas as amostras que a cobrem e calcula a probabilidade de cada genótipo possível.
2. **`bcftools call`**: a partir dessas probabilidades, decide em quais posições existe uma variante e qual é o genótipo de cada indivíduo.

Vamos rodar primeiro e entender as opções enquanto o comando trabalha:

```bash
bcftools mpileup -a AD,DP,SP -d 1000 -Ou -f $REF ~/alinhamento/*.bam | \
bcftools call -m -v -f GQ,GP -Oz -o ~/vcf/acrocomia_bcftools.vcf.gz
```

**`bcftools mpileup`:**

| Opção | Função |
|---|---|
| `-a AD,DP,SP` | Acrescenta informações ao VCF: profundidade de cada alelo (`AD`), profundidade total por amostra (`DP`) e viés de fita (`SP`) |
| `-d 1000` | Número máximo de reads lidos por amostra em cada posição. O padrão é 250, baixo demais para os locos de GBS |
| `-Ou` | Saída em formato binário sem compressão (`u`), a mais rápida para passar ao próximo programa pelo *pipe* |
| `-f` | O genoma de referência |

**`bcftools call`:**

| Opção | Função |
|---|---|
| `-m` | Usa o modelo multialélico, que permite mais de um alelo alternativo por posição |
| `-v` | Grava só as posições variantes, ignorando as idênticas à referência |
| `-f GQ,GP` | Acrescenta a qualidade do genótipo (`GQ`) e as probabilidades de cada genótipo (`GP`) |
| `-Oz` | Saída em VCF comprimido (`.vcf.gz`) |
| `-o` | Arquivo de saída |

Em resumo: empilhamos os reads em cada posição do cromossomo e chamamos uma variante onde há reads suficientes sustentando a presença de um alelo alternativo.

Quando terminar, vamos criar o índice do VCF, necessário para várias das próximas etapas:

```bash
bcftools index ~/vcf/acrocomia_bcftools.vcf.gz
```

### C. Explorando o VCF

O VCF (*Variant Call Format*) tem um **cabeçalho**, com as linhas que começam com `##`, e uma **tabela**, com uma linha por variante.

```bash
# As últimas linhas do cabeçalho
bcftools view -h ~/vcf/acrocomia_bcftools.vcf.gz | tail -n 15

# Os nomes das amostras
bcftools query -l ~/vcf/acrocomia_bcftools.vcf.gz

# As variantes (setas para navegar, q para sair)
bcftools view -H ~/vcf/acrocomia_bcftools.vcf.gz | less -S
```

Cada linha da tabela tem as seguintes colunas:

| Coluna | Conteúdo |
|---|---|
| CHROM | Cromossomo |
| POS | Posição |
| ID | Identificador da variante (aqui, vazio: `.`) |
| REF | Alelo da referência |
| ALT | Alelo(s) alternativo(s) |
| QUAL | Qualidade da chamada (escala Phred) |
| FILTER | Se a variante passou em filtros aplicados |
| INFO | Informações da variante em toda a população (profundidade total, etc.) |
| FORMAT | A ordem dos campos nas colunas das amostras (por exemplo, `GT:PL:DP:AD:GQ`) |
| Amostras | Uma coluna por indivíduo, com os campos descritos em FORMAT |

No campo `GT` (genótipo), `0` representa o alelo da referência e `1` o alelo alternativo: `0/0` é homozigoto para a referência, `0/1` é heterozigoto, `1/1` é homozigoto para o alternativo e `./.` é um dado faltante.

O `bcftools query` permite extrair só os campos que nos interessam, em formato de tabela:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n' ~/vcf/acrocomia_bcftools.vcf.gz | head -n 20
```

> [!TIP]
> No `bcftools query`, `%CHROM`, `%POS` etc. são as colunas, e o que está entre colchetes (`[\t%GT]`) é repetido para cada amostra.

E para contar as variantes por tipo:

```bash
# Total de variantes
bcftools view -H ~/vcf/acrocomia_bcftools.vcf.gz | wc -l

# Só SNPs
bcftools view -H -v snps ~/vcf/acrocomia_bcftools.vcf.gz | wc -l

# Só InDels
bcftools view -H -v indels ~/vcf/acrocomia_bcftools.vcf.gz | wc -l

# Um resumo completo
bcftools stats ~/vcf/acrocomia_bcftools.vcf.gz > ~/vcf/acrocomia_bcftools.stats.txt
grep "^SN" ~/vcf/acrocomia_bcftools.stats.txt
```

❓ Quantas variantes foram encontradas no cromossomo 5? Quantas são SNPs e quantas são InDels? Encontre uma variante com QUAL baixa e uma com QUAL alta: o que muda nos genótipos e nas profundidades (`DP` e `AD`) das amostras? Há muitos dados faltantes (`./.`)? Por que isso é esperado em dados de GBS?

### D. Selecionando só os SNPs

Como o nosso foco são os SNPs, vamos criar um VCF só com eles:

```bash
bcftools view -v snps -Oz -o ~/vcf/acrocomia_bcftools_snps.vcf.gz ~/vcf/acrocomia_bcftools.vcf.gz
bcftools index ~/vcf/acrocomia_bcftools_snps.vcf.gz
```

Este arquivo ainda é **bruto** (sem filtros): contém todos os candidatos a SNP, inclusive os de baixa qualidade, com baixa profundidade ou com muitos dados faltantes. A filtragem, com o bcftools e o vcftools, será o tema da próxima prática.

---

## 2. Chamada com genoma de referência: alternativa com o FreeBayes

O FreeBayes analisa os reads de todos os indivíduos ao mesmo tempo. Isso dá a ele um poder estatístico maior: se um alelo aparece com boa qualidade em alguns indivíduos, o programa usa essa evidência para avaliar a mesma posição nos demais, inclusive naqueles com menor cobertura.

Primeiro, ele precisa de um arquivo de texto com a lista de alinhamentos a analisar:

```bash
ls ~/alinhamento/*.bam > ~/vcf/lista_bams.txt
cat ~/vcf/lista_bams.txt
```

Agora rodamos o FreeBayes. Como ele grava um VCF sem compressão, passamos a saída para o `bgzip`, que comprime no formato que o bcftools consegue indexar:

```bash
freebayes -f $REF -L ~/vcf/lista_bams.txt --standard-filters | \
bgzip > ~/vcf/acrocomia_freebayes.vcf.gz

bcftools index ~/vcf/acrocomia_freebayes.vcf.gz
```

> [!TIP]
> - `-f`: o genoma de referência (indexado).
> - `-L`: o arquivo com a lista de alinhamentos.
> - `--standard-filters`: aplica um conjunto de filtros recomendados, por exemplo ignorar reads com qualidade de mapeamento menor que 30 e bases com qualidade menor que 20.
> - `bgzip`: compressão em blocos, uma versão do `gzip` que permite indexar o arquivo e acessar regiões específicas sem descomprimir tudo.

Vamos explorar o resultado com os mesmos comandos de antes:

```bash
bcftools view -H ~/vcf/acrocomia_freebayes.vcf.gz | less -S

bcftools view -H ~/vcf/acrocomia_freebayes.vcf.gz | wc -l
bcftools view -H -v snps ~/vcf/acrocomia_freebayes.vcf.gz | wc -l
```

> [!NOTE]
> O FreeBayes classifica as variantes em mais categorias que o bcftools: além de `snp` e `ins`/`del`, há `mnp` (várias bases vizinhas trocadas ao mesmo tempo) e `complex` (combinações de SNPs e InDels). Você vai encontrar essa informação no campo `TYPE` da coluna INFO.

E também vamos separar só os SNPs:

```bash
bcftools view -v snps -Oz -o ~/vcf/acrocomia_freebayes_snps.vcf.gz ~/vcf/acrocomia_freebayes.vcf.gz
bcftools index ~/vcf/acrocomia_freebayes_snps.vcf.gz
```

### Comparando o bcftools e o FreeBayes

O `bcftools isec` compara dois VCFs e separa os SNPs exclusivos de cada um e os compartilhados:

```bash
bcftools isec -p ~/vcf/comparacao \
    ~/vcf/acrocomia_bcftools_snps.vcf.gz \
    ~/vcf/acrocomia_freebayes_snps.vcf.gz

# O README explica o que há em cada arquivo
cat ~/vcf/comparacao/README.txt

# Quantos SNPs há em cada grupo
for arquivo in ~/vcf/comparacao/000*.vcf; do
    echo "$(basename $arquivo): $(grep -vc '^#' $arquivo) SNPs"
done
```

❓ Quantos SNPs cada programa encontrou? Quantos são compartilhados e quantos são exclusivos de cada um? Os SNPs exclusivos têm qualidade e profundidade parecidas com as dos compartilhados? O que isso sugere sobre a escolha do programa de chamada?

---

## 3. Chamada sem genoma de referência (*de novo*)

Até agora, as ferramentas que usamos compartilham uma premissa fundamental: dependem de um genoma de referência, contra o qual os reads são mapeados para então identificar os SNPs. Porém, na prática da biologia da conservação e da ecologia evolutiva, nem sempre temos um genoma disponível.

Para esses casos, existem algoritmos que montam os locos e chamam as variantes ***de novo***, sem nenhuma referência prévia. A maioria dessas ferramentas foi desenvolvida para bibliotecas de representação reduzida, como o GBS (*Genotyping-by-Sequencing*) e o RADseq. Embora haja uma clara tendência de transição para o sequenciamento de genoma completo (WGS), o GBS continua sendo uma realidade metodológica e econômica importante para viabilizar estudos populacionais, especialmente no contexto latino-americano.

Entre os *pipelines* mais conhecidos para a chamada *de novo* de variantes em bibliotecas de representação reduzida, destacam-se:

1. **[Stacks](https://catchenlab.life.illinois.edu/stacks/)**: uma das ferramentas mais robustas e usadas, baseada na construção de um catálogo de locos a partir do agrupamento de sequências semelhantes.
2. **[ipyrad](https://ipyrad.readthedocs.io/en/master/index.html)**: uma evolução do antigo *PyRAD*, que lida bem com InDels e é otimizado para análises filogenéticas e de introgressão.
3. **[TASSEL](https://tassel.bitbucket.io)**: um software clássico do melhoramento genético de plantas (especialmente milho e soja), com um *pipeline* consolidado para GBS (o *TASSEL-GBS*).

Vamos usar o [`denovo_map.pl`](https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) do Stacks, muito usado em genética de populações de espécies não-modelo, e que tem a vantagem de executar todas as etapas até o VCF final com um único comando.

### A. Visão geral do `denovo_map.pl`

Quando executamos o `denovo_map.pl`, ele roda automaticamente os seguintes módulos, nesta ordem:

1. **`ustacks`** (construção dos locos em cada indivíduo): processa os reads de cada amostra separadamente. Agrupa sequências idênticas em pilhas (*stacks*) e depois junta pilhas com poucas diferenças (possíveis alelos) para formar os locos de cada indivíduo.
2. **`cstacks`** (o catálogo): junta os locos de todos os indivíduos em um **catálogo**, que funciona como um genoma de referência "virtual" para o restante da análise.
3. **`sstacks`** (ancoragem): os locos de cada indivíduo são comparados com o catálogo. Isso garante que o loco 150 do indivíduo A seja comparado exatamente com o loco 150 do indivíduo B.
4. **`tsv2bam` e `gstacks`** (genotipagem e faseamento): reorganizam os dados, chamam os SNPs com um modelo probabilístico e fazem o faseamento dos SNPs de um mesmo loco, construindo micro-haplótipos.
5. **`populations`** (genética de populações): aplica filtros populacionais, calcula estatísticas como $F_{ST}$, $\pi$ e heterozigosidade, e exporta os dados em diversos formatos, entre eles o VCF.

<img width="1165" alt="stacks_pipeline" src="https://github.com/user-attachments/assets/cf71b50b-6c96-4ee2-943f-534b4f3c4075" />

### B. O arquivo de populações (*popmap*)

O Stacks precisa de um arquivo que indique a que população pertence cada amostra, o *popmap*. Ele tem exatamente duas colunas, separadas por **TAB**:

1. **ID da amostra:** o nome do arquivo FASTQ **sem** a extensão. Os nossos arquivos se chamam, por exemplo, `AM_04_trimmed.fq.gz`, então o ID é `AM_04_trimmed`.
2. **População:** um texto sem espaços indicando a origem ou o grupo do indivíduo.

```bash
# Criar a pasta para os resultados do Stacks
mkdir -p ~/stacks_denovo

# Criar o popmap com o nano
nano ~/stacks_denovo/popmap.txt
```

Cole o conteúdo abaixo (lembre-se: barcode e população separados por **TAB**):

```
AM_04_trimmed	amazonas
AM_05_trimmed	amazonas
G109_trimmed	guatemala
G111_trimmed	guatemala
rifania10_trimmed	sao_paulo
rifania8_trimmed	sao_paulo
luz_20m_trimmed	minas
luz_25m_trimmed	minas
veracruz-mex-L1p2_trimmed	mexico
veracruz-mex-L3p2_trimmed	mexico
```

Para salvar: `Ctrl + O` e `Enter`. Para sair: `Ctrl + X`. E, como no arquivo de barcodes, confira os TABs:

```bash
cat -A ~/stacks_denovo/popmap.txt
```

Cada linha deve ter um `^I` entre o nome e a população, e terminar só com `$`.

> [!TIP]
> Uma forma prática de montar o *popmap* é organizar os dados em uma planilha (Excel ou Google Sheets) e copiar as duas colunas: ao colar no terminal, o separador entre elas já é um TAB.

### C. Rodando o `denovo_map.pl`

De preferência dentro de um `screen`, porque esta etapa pode demorar:

```bash
denovo_map.pl \
    -T 1 \
    -M 3 \
    -n 3 \
    --samples ~/trimmed \
    --popmap ~/stacks_denovo/popmap.txt \
    -o ~/stacks_denovo \
    -X "populations: --vcf"
```

| Opção | Função |
|---|---|
| `-T` | Número de núcleos (*threads*) |
| `-M` | Número máximo de diferenças entre duas pilhas para que sejam consideradas alelos do mesmo loco **dentro de um indivíduo** (`ustacks`) |
| `-n` | Número máximo de diferenças entre locos de **indivíduos diferentes** para que sejam unidos no catálogo (`cstacks`) |
| `--samples` | Pasta com os arquivos FASTQ limpos |
| `--popmap` | O arquivo de populações |
| `-o` | Pasta de saída |
| `-X` | Passa opções extras para um módulo específico. Aqui, pede ao `populations` que gere um VCF |

> [!NOTE]
> Outro parâmetro importante é o `-m` do `ustacks`: o número mínimo de reads idênticos para formar uma pilha. O padrão é 3, e é o que estamos usando. Para alterá-lo, use `-X "ustacks: -m 5"`, por exemplo.

### D. Resultados

```bash
ls -lh ~/stacks_denovo
```

Entre os arquivos de saída, destacam-se:

1. `populations.snps.vcf`: os SNPs, em formato VCF.
2. `populations.sumstats.tsv`: estatísticas de diversidade para cada SNP e população.
3. `catalog.fa.gz`: as sequências consenso de todos os locos do catálogo.
4. `denovo_map.log` e `populations.log`: os registros de cada etapa.

```bash
less -S ~/stacks_denovo/populations.snps.vcf

# Quantos SNPs?
grep -vc "^#" ~/stacks_denovo/populations.snps.vcf

# Quantos locos há no catálogo?
zcat ~/stacks_denovo/catalog.fa.gz | grep -c "^>"
```

❓ Quantos SNPs o Stacks encontrou? Como esse número se compara aos do bcftools e do FreeBayes? (Dica: nas chamadas com referência usamos só o cromossomo 5. E no *de novo*, que parte do genoma foi usada?)

❓ O que aparece na coluna CHROM do VCF do Stacks, já que não há cromossomos? E na coluna POS?

> [!NOTE]
> **O Stacks passo a passo.** Neste exercício, usamos o `denovo_map.pl`, que roda automaticamente todas as etapas e entrega o resultado final. É a forma mais rápida e direta. Também é possível rodar cada módulo separadamente (`ustacks`, `cstacks`, `sstacks`, `tsv2bam`, `gstacks` e `populations`). Embora mais trabalhosa, essa abordagem dá controle total sobre cada etapa e permite testar e calibrar os parâmetros de cada uma de forma independente, o que é recomendado para escolher os melhores valores de `-M` e `-n` para a sua espécie.
