# 💻5: Filtros do arquivo VCF

A última etapa antes das análises populacionais é filtrar o VCF. O VCF que geramos é **bruto**: contém todos os candidatos a SNP, inclusive erros de sequenciamento, erros de mapeamento e posições com pouca informação. O objetivo da filtragem é manter só os SNPs confiáveis e informativos.

Filtrar parece uma tarefa simples, mas é uma das etapas mais importantes de toda a análise, e merece tempo e cuidado. Não existe um conjunto de filtros "correto" para todos os projetos: os limites dependem do tipo de biblioteca, da profundidade de sequenciamento, do número de indivíduos e das perguntas do estudo. O mais importante é **conhecer os seus dados** e saber **justificar** cada escolha.

Por isso, vamos trabalhar em três passos:

1. **Diagnóstico:** calcular estatísticas dos SNPs e dos indivíduos.
2. **Filtros um por um:** ver quanto cada filtro remove.
3. **Filtragem final:** aplicar todos os filtros juntos.

---

## 0. Lembrando do arquivo VCF

O VCF (*Variant Call Format*) é o formato padrão para armazenar variantes genéticas. Ele tem dois blocos:

1. **Cabeçalho:** as linhas que começam com `##`. Funciona como um dicionário: define as abreviações usadas no arquivo, os filtros aplicados, os formatos dos dados e os programas usados para gerá-lo.
2. **Tabela de variantes:** começa com a linha `#CHROM`, que nomeia as colunas. A partir dela, cada linha descreve uma variante encontrada na população.

| Campo | Descrição |
| :--- | :--- |
| **`#CHROM`** | O cromossomo (ou *scaffold*) onde está a variante. |
| **`POS`** | A posição da variante no cromossomo, em pares de bases. |
| **`ID`** | Identificador da variante (em humanos, por exemplo, o número do dbSNP). Em espécies não-modelo, geralmente fica vazio (`.`). |
| **`REF`** | O alelo presente no genoma de referência. |
| **`ALT`** | O alelo (ou alelos) alternativo encontrado nas amostras. |
| **`QUAL`** | Qualidade da chamada da variante, em escala Phred. Quanto maior, maior a probabilidade de que a variante seja real e não um erro técnico. |
| **`FILTER`** | Se a variante passou em filtros aplicados (`PASS`) ou o nome do filtro em que falhou. Quando nenhum filtro foi aplicado, aparece `.`. |
| **`INFO`** | Informações da variante considerando todas as amostras (por exemplo, a profundidade total). |
| **`FORMAT`** | A ordem dos campos nas colunas das amostras (por exemplo, `GT:PL:DP:AD:GQ`). |
| **`<amostras>`** | Uma coluna por indivíduo, com o genótipo (`0/0` homozigoto para a referência, `0/1` heterozigoto, `1/1` homozigoto para o alternativo, `./.` faltante) e as métricas descritas em `FORMAT`. |

<p align="center">
 <img width="1504" alt="vcf_format" src="https://github.com/user-attachments/assets/cbda7aa1-e532-4f98-a0ec-a1029d461076" />
  <br>
  <sub><b>Figura 1.</b> Formato do arquivo VCF.</sub>
</p>

Nesta prática, vamos usar duas ferramentas:

1. **[VCFtools](https://vcftools.github.io/man_latest.html)**: para calcular as estatísticas e aplicar os filtros.
2. **[bcftools](https://samtools.github.io/bcftools/bcftools.html)**: para tarefas auxiliares, como contar SNPs, indexar o VCF e gerar o resumo final.

---

## 1. Preparação

Vamos criar uma pasta para os filtros dentro da pasta dos VCFs e definir duas variáveis, para não precisar digitar os caminhos completos a cada comando:

```bash
mkdir -p ~/vcf/filtros

# O VCF de entrada: os SNPs brutos do bcftools
VCF=~/vcf/acrocomia_bcftools_snps.vcf.gz

# O prefixo dos arquivos de saída
OUT=~/vcf/filtros/acrocomia
```

> [!IMPORTANT]
> Lembre-se: as variáveis só existem na sessão atual do terminal. Se você sair do servidor e entrar de novo, defina outra vez `VCF` e `OUT` (e as variáveis dos limites, na seção 3).

Quantos SNPs temos antes de qualquer filtro?

```bash
bcftools view -H $VCF | wc -l
```

---

## 2. Diagnóstico: conhecendo os dados

Antes de escolher os filtros, vamos calcular algumas estatísticas. Os comandos abaixo **não alteram o VCF**: apenas geram arquivos de texto com as métricas.

### A. Estatísticas por SNP

```bash
# Qualidade de cada SNP (QUAL)
vcftools --gzvcf $VCF --site-quality --out $OUT

# Profundidade média de cada SNP (média entre os indivíduos)
vcftools --gzvcf $VCF --site-mean-depth --out $OUT

# Proporção de dados faltantes em cada SNP
vcftools --gzvcf $VCF --missing-site --out $OUT

ls -lh ~/vcf/filtros
```

Cada comando gera um arquivo com uma linha por SNP: `.lqual` (qualidade), `.ldepth.mean` (profundidade média) e `.lmiss` (dados faltantes). Como são milhares de linhas, vamos resumir cada um. Para isso, criamos uma pequena função no Bash, que calcula o mínimo, os percentis, a mediana, a média e o máximo de uma coluna (como o `summary()` do R):

```bash
resumo() {
    tail -n +2 "$1" | cut -f "$2" | grep -v "nan" | sort -g | \
    awk '{ v[NR] = $1; soma += $1 }
         END {
             printf "N de SNPs:     %d\n", NR
             printf "Mínimo:        %s\n", v[1]
             printf "Percentil 5%%:  %s\n", v[int(NR * 0.05) + 1]
             printf "Mediana:       %s\n", v[int((NR + 1) / 2)]
             printf "Média:         %.2f\n", soma / NR
             printf "Percentil 95%%: %s\n", v[int(NR * 0.95)]
             printf "Máximo:        %s\n", v[NR]
         }'
}
```

> [!TIP]
> **Como funciona a função `resumo`:** ela recebe um arquivo (`$1`) e o número de uma coluna (`$2`).
> - `tail -n +2` pula a primeira linha (o cabeçalho).
> - `cut -f` pega só a coluna desejada.
> - `grep -v "nan"` remove valores indefinidos.
> - `sort -g` ordena os valores do menor para o maior.
> - O `awk` guarda os valores ordenados e, no final (`END`), imprime o resumo.
>
> Assim como as variáveis, a função só existe na sessão atual do terminal.

Agora, resumimos cada estatística:

```bash
# Qualidade dos SNPs (coluna 3 do .lqual)
resumo ${OUT}.lqual 3

# Profundidade média dos SNPs (coluna 3 do .ldepth.mean)
resumo ${OUT}.ldepth.mean 3

# Proporção de dados faltantes nos SNPs (coluna 6 do .lmiss)
resumo ${OUT}.lmiss 6
```

Antes de responder às perguntas, vamos entender o que cada estatística nos diz:

**Qualidade (QUAL):** a confiança de que o SNP é real. Em escala Phred, QUAL = 30 significa uma probabilidade de erro de 1 em 1000.

**Profundidade (DP):** quantos reads cobrem o SNP, em média, em cada indivíduo. Precisamos de dois limites:
- um **mínimo**, porque com poucos reads não conseguimos distinguir um heterozigoto de um erro de sequenciamento;
- um **máximo**, porque posições com profundidade muito acima da média geralmente são **regiões repetitivas ou parálogas**, colapsadas em uma única cópia no genoma de referência. Os reads de todas as cópias mapeiam no mesmo lugar, e as diferenças entre as cópias aparecem como falsos SNPs, quase sempre heterozigotos em todos os indivíduos. Uma regra prática é usar um máximo de cerca de **2 vezes a profundidade média** (ou o percentil 95).

**Dados faltantes:** a proporção de indivíduos sem genótipo naquele SNP. Em GBS, esse valor é naturalmente mais alto que em genoma completo, porque nem todos os locos são sequenciados em todos os indivíduos (por exemplo, por mutações no sítio de restrição ou pela amostragem aleatória dos fragmentos).

❓ Qual é a qualidade mediana dos SNPs? E a profundidade mediana? Há SNPs com profundidade muito acima da média? Que proporção de dados faltantes é típica nos nossos SNPs?

### B. Estatísticas por indivíduo

```bash
# Proporção de dados faltantes em cada indivíduo
vcftools --gzvcf $VCF --missing-indv --out $OUT

# Profundidade média de cada indivíduo
vcftools --gzvcf $VCF --depth --out $OUT

# Heterozigosidade e coeficiente de endogamia (F) de cada indivíduo
vcftools --gzvcf $VCF --het --out $OUT
```

Como temos só 10 indivíduos, em vez de resumir, podemos ver as tabelas inteiras, ordenadas pela coluna de interesse:

```bash
# Dados faltantes (ordenado pela coluna 5: F_MISS)
sort -k5,5g ${OUT}.imiss | column -t

# Profundidade média (ordenado pela coluna 3: MEAN_DEPTH)
sort -k3,3g ${OUT}.idepth | column -t

# Heterozigosidade (ordenado pela coluna 5: F)
sort -k5,5g ${OUT}.het | column -t
```

> [!TIP]
> No `sort`, `-k5,5` indica que a ordenação deve usar a coluna 5, e o `g` faz a ordenação numérica (inclusive com decimais e notação científica).

Indivíduos com muito mais dados faltantes ou com profundidade muito menor que os demais provavelmente tiveram problemas na extração do DNA, na biblioteca ou no sequenciamento. Mantê-los na análise pode distorcer as estimativas de diversidade e de estrutura populacional.

❓ Qual é a profundidade média dos indivíduos? Há algum indivíduo com muito mais dados faltantes ou com profundidade muito menor que os outros? Você removeria algum? Por quê?

> [!NOTE]
> **Sobre a heterozigosidade.** Um indivíduo com heterozigosidade muito maior que os demais (F muito negativo) pode indicar **contaminação** da amostra com DNA de outro indivíduo. Use esta estatística apenas como referência: as nossas 10 amostras vêm de populações muito distantes (do México ao Brasil), então o F tende a sair alto para todos, porque ele é calculado como se todos pertencessem a uma única população. Esse é o **efeito Wahlund**: a estrutura populacional reduz a heterozigosidade observada em relação à esperada, mesmo sem endogamia.

---

## 3. Filtros um por um

### A. Definindo os limites

Com base no diagnóstico, definimos os limites de cada filtro. Vamos guardá-los em variáveis, para usá-los tanto nos testes individuais quanto na filtragem final:

```bash
# Filtros dos SNPs
MAF=0.05            # frequência mínima do alelo menor
MIN_Q=30            # qualidade mínima do SNP (QUAL)
MIN_MEANDP=10       # profundidade média mínima do SNP
MAX_MEANDP=100      # profundidade média máxima do SNP
MAX_MISSING=0.75    # proporção mínima de indivíduos COM genótipo (ver explicação abaixo)

# Filtros dos genótipos
MIN_DP=5            # profundidade mínima de cada genótipo
MIN_GQ=20           # qualidade mínima de cada genótipo

# Filtro dos indivíduos
MAX_MISS_IND=0.5    # proporção máxima de dados faltantes por indivíduo
```

Por que esses valores?

- **`MAF=0.05`:** remove os SNPs em que o alelo menos frequente aparece em menos de 5% das cópias. Alelos muito raros são difíceis de distinguir de erros de sequenciamento, e trazem pouca informação para a maioria das análises populacionais. É o valor mais usado na literatura.
- **`MIN_Q=30`:** probabilidade de erro de 1 em 1000 na chamada do SNP.
- **`MIN_MEANDP` e `MAX_MEANDP`:** definidos a partir do resumo da profundidade dos SNPs (seção 2A).
- **`MAX_MISSING`:** no VCFtools, este filtro funciona "ao contrário": indica a proporção mínima de indivíduos **com** genótipo. `0.75` significa que toleramos até 25% de dados faltantes por SNP.
- **`MIN_DP` e `MIN_GQ`:** um GQ de 20 corresponde a 99% de confiança no genótipo.
- **`MAX_MISS_IND`:** definido a partir da tabela de dados faltantes por indivíduo (seção 2B).

> [!NOTE]
> **Sobre o MAF com 10 indivíduos.** Com 10 indivíduos diploides, temos 20 cópias de cada SNP. Um MAF de 0,05 corresponde a 1 cópia em 20: ou seja, com o nosso número de amostras, este filtro remove muito pouco. Em estudos populacionais, com dezenas ou centenas de indivíduos, o MAF de 0,05 é o padrão recomendado.

### B. Testando cada filtro separadamente

Agora vamos aplicar cada filtro **sozinho** sobre o VCF bruto e ver quantos SNPs sobram. Como só queremos a contagem, não geramos nenhum VCF: o VCFtools informa o resultado no seu relatório, e o `grep` mostra só a linha que interessa.

```bash
echo "Só bialélicos:"
vcftools --gzvcf $VCF --min-alleles 2 --max-alleles 2 --out ${OUT}_teste 2>&1 | grep "After filtering, kept"

echo "MAF:"
vcftools --gzvcf $VCF --maf $MAF --out ${OUT}_teste 2>&1 | grep "After filtering, kept"

echo "Qualidade do SNP:"
vcftools --gzvcf $VCF --minQ $MIN_Q --out ${OUT}_teste 2>&1 | grep "After filtering, kept"

echo "Profundidade média do SNP:"
vcftools --gzvcf $VCF --min-meanDP $MIN_MEANDP --max-meanDP $MAX_MEANDP --out ${OUT}_teste 2>&1 | grep "After filtering, kept"

echo "Dados faltantes por SNP:"
vcftools --gzvcf $VCF --max-missing $MAX_MISSING --out ${OUT}_teste 2>&1 | grep "After filtering, kept"
```

> [!TIP]
> O VCFtools escreve o seu relatório na "saída de erro" (*stderr*), e não na saída normal. O `2>&1` junta as duas saídas, para que o `grep` consiga ler o relatório.

Os filtros de **genótipo** funcionam de outro jeito: eles não removem SNPs, e sim **transformam em dado faltante** os genótipos que não passam no filtro. Por isso, o efeito deles só aparece quando combinados com o filtro de dados faltantes:

```bash
echo "Dados faltantes por SNP, SEM filtro de genótipo:"
vcftools --gzvcf $VCF --max-missing $MAX_MISSING --out ${OUT}_teste 2>&1 | grep "After filtering, kept"

echo "Dados faltantes por SNP, COM filtro de genótipo:"
vcftools --gzvcf $VCF --minDP $MIN_DP --minGQ $MIN_GQ --max-missing $MAX_MISSING --out ${OUT}_teste 2>&1 | grep "After filtering, kept"
```

❓ Qual filtro remove mais SNPs? Por que o filtro de genótipo, combinado com o de dados faltantes, remove mais SNPs do que o de dados faltantes sozinho?

### C. Indivíduos a remover

Vamos gerar a lista de indivíduos com mais dados faltantes do que o limite que definimos:

```bash
awk -v limite=$MAX_MISS_IND 'NR > 1 && $5 > limite {print $1}' ${OUT}.imiss > ${OUT}_remover.txt

echo "Indivíduos a remover:"
cat ${OUT}_remover.txt
```

> [!TIP]
> - `-v limite=$MAX_MISS_IND` passa o valor da variável do Bash para dentro do `awk`.
> - `NR > 1` pula a linha do cabeçalho.
> - `$5 > limite` seleciona as linhas em que a coluna 5 (F_MISS) é maior que o limite, e `print $1` imprime o nome do indivíduo.
>
> Se nenhum indivíduo passar do limite, o arquivo fica vazio, e nenhuma amostra é removida.

---

## 4. Filtragem final: todos os filtros juntos

Agora aplicamos todos os filtros em um único comando e geramos o VCF final:

```bash
vcftools --gzvcf $VCF \
    --remove ${OUT}_remover.txt \
    --remove-indels \
    --min-alleles 2 --max-alleles 2 \
    --maf $MAF \
    --minQ $MIN_Q \
    --min-meanDP $MIN_MEANDP --max-meanDP $MAX_MEANDP \
    --minDP $MIN_DP --minGQ $MIN_GQ \
    --max-missing $MAX_MISSING \
    --recode --recode-INFO-all --stdout 2> ${OUT}_filtrado.log | \
bgzip > ${OUT}_filtrado.vcf.gz

bcftools index ${OUT}_filtrado.vcf.gz
```

| Opção | Tipo | Função |
|---|---|---|
| `--gzvcf` | Entrada | O VCF de entrada, comprimido (`.vcf.gz`) |
| `--remove` | Indivíduos | Remove os indivíduos listados no arquivo |
| `--remove-indels` | SNPs | Remove InDels. O nosso VCF já tem só SNPs, mas o filtro garante isso caso o VCF de entrada tenha outras variantes |
| `--min-alleles 2 --max-alleles 2` | SNPs | Mantém só os SNPs bialélicos |
| `--maf` | SNPs | Frequência mínima do alelo menor |
| `--minQ` | SNPs | Qualidade mínima do SNP |
| `--min-meanDP` / `--max-meanDP` | SNPs | Profundidade média mínima e máxima do SNP |
| `--minDP` | Genótipos | Profundidade mínima de **cada genótipo**. Os que não passam viram dado faltante |
| `--minGQ` | Genótipos | Qualidade mínima de cada genótipo. Os que não passam viram dado faltante |
| `--max-missing` | SNPs | Proporção mínima de indivíduos com genótipo (0 = aceita tudo faltante, 1 = nenhum faltante) |
| `--recode` | Saída | Gera um novo VCF com os dados filtrados |
| `--recode-INFO-all` | Saída | Mantém todas as informações da coluna INFO no novo VCF |
| `--stdout` | Saída | Envia o VCF para a saída padrão, para comprimi-lo com o `bgzip` pelo *pipe* |
| `2>` | Saída | Salva o relatório do VCFtools em um arquivo `.log` |

Quantos SNPs e quantos indivíduos sobraram?

```bash
cat ${OUT}_filtrado.log
```

---

## 5. Estatísticas finais

Vamos comparar o VCF antes e depois da filtragem:

```bash
echo "SNPs antes:"
bcftools view -H $VCF | wc -l

echo "SNPs depois:"
bcftools view -H ${OUT}_filtrado.vcf.gz | wc -l

echo "Amostras no VCF final:"
bcftools query -l ${OUT}_filtrado.vcf.gz

# Resumo completo do VCF final
bcftools stats ${OUT}_filtrado.vcf.gz > ${OUT}_filtrado.stats.txt
grep "^SN" ${OUT}_filtrado.stats.txt
```

E como ficaram os indivíduos depois da filtragem?

```bash
vcftools --gzvcf ${OUT}_filtrado.vcf.gz --missing-indv --out ${OUT}_filtrado
vcftools --gzvcf ${OUT}_filtrado.vcf.gz --depth --out ${OUT}_filtrado

sort -k5,5g ${OUT}_filtrado.imiss | column -t
sort -k3,3g ${OUT}_filtrado.idepth | column -t
```

❓ Quantos SNPs havia antes e quantos sobraram depois? Que porcentagem foi mantida?

❓ Os dados faltantes por indivíduo diminuíram? Por quê?

❓ Se você precisasse de mais SNPs, que filtro afrouxaria? E se precisasse de SNPs mais confiáveis, qual endureceria? Quais seriam as consequências de cada escolha para as análises populacionais?

O arquivo `acrocomia_filtrado.vcf.gz` é o ponto de partida para as análises de genética de populações: diversidade, estrutura, diferenciação e muito mais. 🌴
