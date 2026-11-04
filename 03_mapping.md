# 💻 3. Alinhamento e mapeamento de sequências

Agora que temos nossos dados limpos e filtrados, o próximo passo é o **mapeamento**: comparar cada read de cada amostra com o genoma de referência para descobrir de qual região do cromossomo aquela sequência veio.

---

## 1. Genoma de referência

### A. Baixar o cromossomo 5

O genoma de referência de *Acrocomia aculeata* está disponível no NCBI (montagem `GCA_058219215.1`). O genoma completo tem cerca de 2 Gb distribuídos em 15 cromossomos, e poderia ser baixado inteiro pelo FTP do NCBI:

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/058/219/215/GCA_058219215.1_ASM5821921v1/GCA_058219215.1_ASM5821921v1_genomic.fna.gz
```

> [!NOTE]
> **Não vamos rodar o download do genoma completo.** Para que a prática caiba no nosso tempo, vamos trabalhar só com o **cromossomo 5**. No final deste roteiro há uma seção explicando como trabalhar com o genoma inteiro.

No NCBI, cada cromossomo tem o seu próprio número de acesso. O do cromossomo 5 é `CM177308.1`. Vamos baixar só ele:

```bash
# Criar a pasta para a referência
mkdir -p ~/referencia

# Baixar o cromossomo 5 do NCBI
wget -O ~/referencia/chr5_ncbi.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CM177308.1&rettype=fasta"
```

> [!TIP]
> Se o download falhar ou o arquivo vier muito pequeno, espere alguns segundos e rode o comando de novo. O NCBI limita o número de downloads simultâneos vindos de um mesmo servidor.

### B. Conhecer o genoma de referência

O genoma é um arquivo FASTA. Vamos explorá-lo com comandos do Linux.

❓ Antes de olhar: como vocês esperam que seja este arquivo?

```bash
# As primeiras linhas: o cabeçalho e o começo da sequência
head -n 3 ~/referencia/chr5_ncbi.fa

# Quantas sequências há no arquivo? (cada sequência começa com ">")
grep -c "^>" ~/referencia/chr5_ncbi.fa

# Qual é o cabeçalho completo?
grep "^>" ~/referencia/chr5_ncbi.fa

# Tamanho do arquivo
ls -lh ~/referencia/chr5_ncbi.fa
```

Agora vamos calcular algumas características do cromossomo. Em todos os comandos abaixo, o `grep -v "^>"` remove o cabeçalho, para contarmos só as bases:

```bash
# Comprimento total do cromossomo (em pb)
grep -v "^>" ~/referencia/chr5_ncbi.fa | tr -d '\n' | wc -c

# Quantas bases G e C?
grep -v "^>" ~/referencia/chr5_ncbi.fa | tr -cd 'GCgc' | wc -c

# Quantas bases N (regiões não resolvidas da montagem)?
grep -v "^>" ~/referencia/chr5_ncbi.fa | tr -cd 'Nn' | wc -c

# Quantas bases em letra minúscula?
grep -v "^>" ~/referencia/chr5_ncbi.fa | tr -cd 'acgtn' | wc -c
```

> [!TIP]
> - `tr -d '\n'` apaga as quebras de linha, juntando a sequência em uma só linha.
> - `tr -cd 'GCgc'` apaga tudo **menos** (`-c`, de complemento) as letras indicadas.
> - `wc -c` conta os caracteres.

❓ Qual é o comprimento do cromossomo 5? Qual é a porcentagem de GC? Há regiões com N? Há bases em letra minúscula? (Dica: pesquise o que é *soft-masking* de regiões repetitivas.)

### C. Renomear o cromossomo

O cabeçalho que vem do NCBI é longo. Os programas usam a primeira palavra do cabeçalho como nome da sequência (aqui, `CM177308.1`), e mais adiante vamos precisar digitar esse nome várias vezes. Uma boa prática é trocar por algo simples, como `Chr5`:

```bash
# Trocar o cabeçalho por ">Chr5" e salvar em um arquivo novo
sed 's/^>.*/>Chr5/' ~/referencia/chr5_ncbi.fa > ~/referencia/chr5.fa

# Conferir
head -n 1 ~/referencia/chr5.fa

# Apagar o arquivo original, que não vamos mais usar
rm ~/referencia/chr5_ncbi.fa
```

> [!TIP]
> No `sed`, `s/padrão/substituição/` troca o texto: `^>.*` significa "uma linha que começa com `>` seguida de qualquer coisa", e a linha inteira é substituída por `>Chr5`.

### D. Indexar a referência

O índice é uma estrutura que permite procurar sequências dentro do genoma muito rapidamente, sem percorrê-lo do começo ao fim a cada read. Cada programa usa o seu próprio tipo de índice:

```bash
# Índice do samtools (.fai): usado pelo samtools, bcftools e freebayes
samtools faidx ~/referencia/chr5.fa

# Índice do bwa: usado no alinhamento (leva alguns minutos)
bwa index ~/referencia/chr5.fa

# Dicionário de sequências (.dict): exigido pelo GATK e pelo Picard
picard CreateSequenceDictionary R=~/referencia/chr5.fa O=~/referencia/chr5.dict

# Ver os arquivos criados
ls -lh ~/referencia
```

Você verá arquivos com as extensões `.fai`, `.dict`, `.amb`, `.ann`, `.bwt`, `.pac` e `.sa`.

> [!WARNING]
> Os índices precisam ficar sempre **na mesma pasta** que o arquivo `.fa`. Se você mover o FASTA, leve todos os índices junto, senão os programas não vão encontrá-los.

❓ Olhe o conteúdo do arquivo `.fai` com `cat ~/referencia/chr5.fa.fai`. O que significa cada coluna? O comprimento bate com o que calculamos na parte B?

---

## 2. Mapeamento com o BWA-MEM: uma amostra, passo a passo

Vamos usar o [bwa mem](https://github.com/lh3/bwa) para alinhar os reads ao cromossomo 5. Primeiro, faremos tudo com **uma amostra** (`AM_04`), para entender cada etapa e cada formato de arquivo.

```bash
# Criar as pastas para os alinhamentos e as estatísticas
mkdir -p ~/alinhamento/stats

# Guardar o caminho da referência em uma variável, para não ter que digitá-lo toda vez
REF=~/referencia/chr5.fa
echo $REF
```

> [!IMPORTANT]
> Variáveis como `REF` só existem na sessão atual do terminal. Se você sair do servidor e entrar de novo, rode outra vez `REF=~/referencia/chr5.fa` antes de continuar.

### A. Alinhar e gerar o arquivo SAM

```bash
bwa mem -t 1 \
    -R "@RG\tID:AM_04\tSM:AM_04\tPL:ILLUMINA" \
    $REF \
    ~/trimmed/AM_04_trimmed.fq.gz \
    > ~/alinhamento/AM_04.sam
```

> [!TIP]
> - `-t 1`: número de núcleos (*threads*).
> - `-R`: o *read group*, que adiciona metadados ao alinhamento: identificador (`ID`), nome da amostra (`SM`) e plataforma (`PL`). Programas como o GATK e o freebayes usam o `SM` para saber a que amostra pertence cada read.
> - `$REF`: a referência indexada. O bwa encontra os arquivos de índice sozinho, desde que estejam na mesma pasta.
> - `>`: redireciona a saída do programa, que normalmente apareceria na tela, para um arquivo.

### B. Explorar o arquivo SAM

O SAM (*Sequence Alignment/Map*) é um arquivo de texto. Vamos abri-lo com o `less` (use as setas para navegar e `q` para sair). O `-S` evita que as linhas longas quebrem na tela:

```bash
less -S ~/alinhamento/AM_04.sam
```

As linhas que começam com `@` formam o **cabeçalho**. Depois, cada linha é um read, com 11 colunas obrigatórias:

| Coluna | Nome | Conteúdo |
|---|---|---|
| 1 | QNAME | Nome do read |
| 2 | FLAG | Código numérico que descreve o alinhamento (mapeado ou não, fita, etc.) |
| 3 | RNAME | Cromossomo onde o read alinhou |
| 4 | POS | Posição de início do alinhamento |
| 5 | MAPQ | Qualidade do mapeamento: a confiança de que o read está no lugar certo |
| 6 | CIGAR | Como o read alinha: `M` (alinhado), `I` (inserção), `D` (deleção), `S` (*soft-clip*) |
| 7–9 | RNEXT, PNEXT, TLEN | Informações do par (vazias em *single-end*) |
| 10 | SEQ | Sequência do read |
| 11 | QUAL | Qualidade de cada base |

❓ O que aparece no cabeçalho (linhas `@SQ`, `@RG` e `@PG`)? Encontre um read com MAPQ igual a 0. O que isso significa? Encontre um read com `*` na coluna 3. O que aconteceu com ele?

> [!TIP]
> Para decifrar o valor da coluna FLAG, use o [Explain SAM Flags](https://broadinstitute.github.io/picard/explain-flags.html).

### C. Converter para BAM e ordenar

O SAM ocupa muito espaço. Na prática, convertemos para **BAM**, a versão binária e comprimida do SAM. Ao mesmo tempo, **ordenamos** os reads pela posição no cromossomo, como exigem os programas das próximas etapas:

```bash
samtools sort -o ~/alinhamento/AM_04.sorted.bam ~/alinhamento/AM_04.sam
samtools index ~/alinhamento/AM_04.sorted.bam
```

O BAM não é texto, então para olhar o conteúdo usamos o `samtools view`:

```bash
samtools view ~/alinhamento/AM_04.sorted.bam | head
```

❓ O que mudou em relação ao SAM? Por que os reads agora aparecem em outra ordem?

### D. Estatísticas do mapeamento

Antes de filtrar, vamos ver como foi o alinhamento:

```bash
samtools flagstat ~/alinhamento/AM_04.sorted.bam
```

❓ Quantos reads há no total? Que porcentagem mapeou no cromossomo 5? Esse valor é alto ou baixo, considerando que a nossa referência é só **um** dos 15 cromossomos?

### E. Filtrar e converter para CRAM

Agora vamos remover os reads não mapeados (`-F 4`) e os reads com qualidade de mapeamento menor que 20 (`-q 20`), e salvar o resultado no formato **CRAM**.

O CRAM comprime ainda mais que o BAM porque **não guarda a sequência completa de cada read**: guarda apenas as diferenças em relação à referência. Por isso, ele precisa da referência tanto para ser criado quanto para ser lido (`-T $REF`).

```bash
samtools view -F 4 -q 20 -C -T $REF -o ~/alinhamento/AM_04.cram ~/alinhamento/AM_04.sorted.bam
samtools index ~/alinhamento/AM_04.cram
```

> [!TIP]
> - `-F 4`: exclui os reads com o FLAG 4 ("não mapeado").
> - `-q 20`: mantém só os reads com MAPQ ≥ 20.
> - `-C`: grava a saída em formato CRAM.
> - `-T`: a referência usada para comprimir.

Quantos reads o filtro removeu?

```bash
samtools view -c ~/alinhamento/AM_04.sorted.bam
samtools view -c -T $REF ~/alinhamento/AM_04.cram
```

### F. Comparar SAM, BAM e CRAM

Para comparar só o efeito do formato, vamos criar também um CRAM **sem filtro**, com exatamente os mesmos reads do SAM e do BAM:

```bash
samtools view -C -T $REF -o ~/alinhamento/AM_04.sorted.cram ~/alinhamento/AM_04.sorted.bam

ls -lh ~/alinhamento/AM_04*
```

| Arquivo | Conteúdo |
|---|---|
| `AM_04.sam` | Todos os reads, em texto |
| `AM_04.sorted.bam` | Todos os reads, binário comprimido |
| `AM_04.sorted.cram` | Todos os reads, comprimido em relação à referência |
| `AM_04.cram` | Só os reads que passaram no filtro: **arquivo final** |

❓ Quantas vezes o BAM é menor que o SAM? E o CRAM, em relação ao BAM? Pensando em um projeto com centenas de amostras de genoma completo, por que isso importa?

Feita a comparação, vamos manter só o arquivo final e apagar os intermediários:

```bash
rm ~/alinhamento/AM_04.sam ~/alinhamento/AM_04.sorted.bam* ~/alinhamento/AM_04.sorted.cram
ls -lh ~/alinhamento
```

---

## 3. Mapeamento de todas as amostras

Agora que entendemos cada etapa, vamos juntar tudo em um único loop, como faríamos na vida real. Para cada amostra, o loop:

1. Alinha com o `bwa mem` e ordena com o `samtools sort`, em um BAM temporário (sem passar pelo SAM).
2. Calcula as estatísticas de mapeamento com o `flagstat`.
3. Filtra e converte para CRAM.
4. Indexa o CRAM.
5. Apaga o BAM temporário.

```bash
REF=~/referencia/chr5.fa

for arquivo in ~/trimmed/*_trimmed.fq.gz; do
    amostra=$(basename "$arquivo" _trimmed.fq.gz)
    echo "Processando a amostra: $amostra"

    # 1. Alinhar e ordenar
    bwa mem -t 1 \
        -R "@RG\tID:${amostra}\tSM:${amostra}\tPL:ILLUMINA" \
        $REF "$arquivo" | \
    samtools sort -o ~/alinhamento/${amostra}.tmp.bam

    # 2. Estatísticas de mapeamento (antes do filtro)
    samtools flagstat ~/alinhamento/${amostra}.tmp.bam > ~/alinhamento/stats/${amostra}_flagstat.txt

    # 3. Filtrar e converter para CRAM
    samtools view -F 4 -q 20 -C -T $REF \
        -o ~/alinhamento/${amostra}.cram \
        ~/alinhamento/${amostra}.tmp.bam

    # 4. Indexar o CRAM
    samtools index ~/alinhamento/${amostra}.cram

    # 5. Apagar o BAM temporário
    rm ~/alinhamento/${amostra}.tmp.bam

    echo "Amostra $amostra finalizada."
done
```

> [!TIP]
> Aqui o `bwa mem` não grava nenhum SAM: a saída vai direto, pelo *pipe* (`|`), para o `samtools sort`. Isso economiza espaço e tempo.

```bash
ls -lh ~/alinhamento
```

Você deve ter um `.cram` e um `.crai` (o índice) para cada uma das 10 amostras.

---

## 4. Estatísticas do mapeamento

### A. Taxa de mapeamento

Vamos juntar os resultados do `flagstat` de todas as amostras em uma única tabela:

```bash
printf "Amostra\tTotal_reads\tReads_mapeados\tPorcentagem_mapeada\n" > ~/alinhamento/stats/resumo_mapeamento.tsv

for arquivo in ~/alinhamento/stats/*_flagstat.txt; do
    amostra=$(basename "$arquivo" _flagstat.txt)
    total=$(awk '/ primary$/ {print $1}' "$arquivo")
    mapeados=$(awk '/primary mapped/ {print $1}' "$arquivo")
    porcentagem=$(awk -F'[(%]' '/primary mapped/ {print $2}' "$arquivo")
    printf "%s\t%s\t%s\t%s\n" "$amostra" "$total" "$mapeados" "$porcentagem" >> ~/alinhamento/stats/resumo_mapeamento.tsv
done

column -t ~/alinhamento/stats/resumo_mapeamento.tsv
```

> [!TIP]
> - Usamos as linhas *primary* do `flagstat` porque elas contam cada read uma única vez. A linha *in total* também soma os alinhamentos secundários e suplementares, quando um mesmo read alinha em mais de um lugar.
> - O `awk` extrai valores específicos de cada linha, e o `printf` escreve uma linha da tabela por amostra.
> - O `column -t` alinha as colunas para facilitar a leitura na tela.

### B. Cobertura e profundidade

O `samtools coverage` resume, para cada amostra, quantas bases do cromossomo foram cobertas por pelo menos um read e qual é a profundidade média. Ele roda nos arquivos CRAM finais, ou seja, nos reads que vão de fato para a chamada de SNPs:

```bash
printf "Amostra\tCromossomo\tInicio\tFim\tN_reads\tBases_cobertas\tCobertura_pct\tProfundidade_media\tQual_base_media\tMAPQ_medio\n" > ~/alinhamento/stats/resumo_cobertura.tsv

for cram in ~/alinhamento/*.cram; do
    amostra=$(basename "$cram" .cram)
    samtools coverage -H --reference $REF "$cram" | \
        awk -v a="$amostra" '{print a "\t" $0}' >> ~/alinhamento/stats/resumo_cobertura.tsv
done

column -t ~/alinhamento/stats/resumo_cobertura.tsv
```

❓ Que porcentagem do cromossomo 5 está coberta por reads em cada amostra? Por que esse valor é tão baixo no GBS? A profundidade média de todo o cromossomo representa bem a profundidade nos locos que foram sequenciados? (Dica: lembre-se de como funciona uma biblioteca de representação reduzida.)

❓ Como foi o nosso mapeamento? Há amostras muito diferentes das outras?

> [!NOTE]
> **Para bibliotecas de genoma completo (WGS e lcWGS)**, a avaliação do mapeamento precisa ir além destas estatísticas. Também é importante olhar, por exemplo: a taxa de duplicação, a porcentagem de reads corretamente pareados e o tamanho do inserto (em dados *paired-end*), e a distribuição da profundidade ao longo do genoma inteiro, e não só a média. Ferramentas como `samtools stats`, [Qualimap](http://qualimap.conesalab.org/) e [mosdepth](https://github.com/brentp/mosdepth) ajudam nessa avaliação.

---

## 5. Duplicatas de PCR: quando remover?

<!-- Inserir imagem aqui -->

Durante a preparação das bibliotecas, a amplificação por PCR pode gerar cópias de um mesmo fragmento de DNA. Se essas duplicatas não forem removidas, os programas de chamada de variantes as interpretam como evidências independentes, inflando a profundidade de um alelo e distorcendo os genótipos e as frequências alélicas.

> [!IMPORTANT]
> **No GBS *single-end*, não removemos duplicatas.** Em uma biblioteca GBS, todos os reads de um mesmo locus começam exatamente no mesmo sítio de corte da enzima. Como o *single-end* não tem a informação do par para diferenciar os fragmentos, reads de moléculas diferentes ficam idênticos em posição, e o programa os confundiria com duplicatas de PCR, descartando dados reais.
>
> **Em bibliotecas de genoma completo (WGS e lcWGS), marcar e remover duplicatas é uma etapa essencial.** Nelas, o DNA é fragmentado aleatoriamente, então dois reads com exatamente a mesma posição de início e fim muito provavelmente vêm da mesma molécula original. No lcWGS, com poucos reads por posição, uma única duplicata pesa ainda mais nas estimativas.

Fica aqui, como referência, o comando para remover duplicatas com o `MarkDuplicates` do [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) em uma amostra de WGS (**não rode com os nossos dados**):

```bash
picard MarkDuplicates \
    I=amostra.sorted.bam \
    O=amostra.sorted.dedup.bam \
    M=amostra_dup_metrics.txt \
    REMOVE_DUPLICATES=true

samtools index amostra.sorted.dedup.bam
```

> [!TIP]
> - `I` / `O`: BAM de entrada (ordenado) e de saída.
> - `M`: relatório com a porcentagem de reads duplicados.
> - `REMOVE_DUPLICATES=true`: apaga as duplicatas. Sem essa opção, elas são só **marcadas** (FLAG 1024), e a maioria dos programas as ignora mesmo assim.

---

## 📌 Para os seus projetos: trabalhando com o genoma completo

Nesta prática usamos só o cromossomo 5. Em um projeto real, você vai baixar o genoma inteiro e renomear todos os cromossomos.

**1. Baixar o genoma completo pelo FTP do NCBI:**

```bash
mkdir -p ~/referencia_completa

wget -c -O ~/referencia_completa/acrocomia_ncbi.fna.gz \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/058/219/215/GCA_058219215.1_ASM5821921v1/GCA_058219215.1_ASM5821921v1_genomic.fna.gz
```

> [!TIP]
> O `-c` permite **continuar** um download interrompido, em vez de começar do zero. Muito útil para arquivos grandes.

**2. Ver os cabeçalhos de todas as sequências:**

```bash
zcat ~/referencia_completa/acrocomia_ncbi.fna.gz | grep "^>"
```

Além dos 15 cromossomos, uma montagem pode trazer sequências não posicionadas (*scaffolds* ou *contigs*) e, às vezes, os genomas de organelas.

**3. Renomear os cromossomos.** O script abaixo procura o número do cromossomo **dentro de cada cabeçalho** (por exemplo, "chromosome 5" vira `Chr5`). As sequências sem número de cromossomo mantêm o número de acesso original. O script também grava uma tabela com a correspondência entre os nomes antigos e os novos, importante para a reprodutibilidade:

```bash
zcat ~/referencia_completa/acrocomia_ncbi.fna.gz | awk '
/^>/ {
    nome_antigo = substr($1, 2)
    if (match($0, /chromosome [0-9]+/)) {
        nome_novo = "Chr" substr($0, RSTART + 11, RLENGTH - 11)
    } else {
        nome_novo = nome_antigo
    }
    print nome_antigo "\t" nome_novo > "/dev/stderr"
    print ">" nome_novo
    next
}
{ print }
' > ~/referencia_completa/acrocomia.fa 2> ~/referencia_completa/nomes_cromossomos.tsv

# Conferir
grep "^>" ~/referencia_completa/acrocomia.fa
cat ~/referencia_completa/nomes_cromossomos.tsv
```

> [!TIP]
> - `match($0, /chromosome [0-9]+/)` procura o texto "chromosome" seguido de um número. Quando encontra, `RSTART` guarda onde o texto começa e `RLENGTH` o seu tamanho.
> - `substr($0, RSTART + 11, ...)` pula as 11 letras de "chromosome " e pega só o número.
> - Diferentemente de numerar os cromossomos pela ordem em que aparecem no arquivo, esse jeito funciona mesmo se o FASTA não estiver em ordem.

**4. Indexar**, exatamente como fizemos com o cromossomo 5. Com um genoma de 2 Gb, o `bwa index` leva bem mais tempo e precisa de alguns GB de memória, por isso costuma ser rodado como um *job* na fila do servidor (Slurm):

```bash
samtools faidx ~/referencia_completa/acrocomia.fa
bwa index ~/referencia_completa/acrocomia.fa
picard CreateSequenceDictionary R=~/referencia_completa/acrocomia.fa O=~/referencia_completa/acrocomia.dict
```
