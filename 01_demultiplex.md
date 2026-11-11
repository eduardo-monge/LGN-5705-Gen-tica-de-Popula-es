# 💻 1: Demultiplexação e controle de qualidade QC

Nesta prática, vamos aprender a analisar as sequências para verificar a qualidade. 

---
## 1. Demultiplexação com STACKS (`process_radtags`)
Quando sequenciamos bibliotecas de representação reduzida (como RADseq ou GBS), várias amostras são misturadas no mesmo tubo (multiplexadas). O primeiro passo é usar os "códigos de barras" (*barcodes*) para separar quem é quem. Geralmente, as empresas já fornecem os dados separados, mas também podemos fazê-lo nós mesmos. 
<p align="center">
<img width="2498" height="1679" alt="gr1_lrg" src="https://github.com/user-attachments/assets/4c7c1cf9-2000-47ca-8b08-2457b9c36b30" />
  <br>
  <sub><b>Figura 1.</b> Steps in GBS library construction. Fonte: Elshire et. al 2011. DOI:10.1371/journal.pone.0019379 .</sub>
</p>



Existem diferentes programas e métodos, dependendo da biblioteca genômica que usamos. No caso do GBS/RAD-seq, um programa que pode ser utilizado é o módulo `process_radtags` do [STACKS](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).

### 0. Preparar os dados para o programa STACKS (El archivo de Barcodes)
Antes de rodar qualquer comando do STACKS para demultiplexar, o programa precisa de um "mapa" para saber como separar essa "sopa" de sequências. Para que o processo seja um sucesso, precisamos de três elementos fundamentais:

1. **A enzima de restrição:** Saber exatamente qual enzima (ou enzimas) foi usada na bancada do laboratório para cortar o DNA (ex: *EcoRI*, *MspI*, *PstI*).
2. **Os arquivos FASTQ brutos:** Que já baixamos na nossa pasta `/home/genpop/dados_brutos`.
3. **O arquivo de Barcodes:** Um "mapa" de texto que vincula cada código de barras com o nome real da amostra.

### 📝 Como deve ser a planilha de amostras?

O STACKS é extremamente rigoroso com o formato deste arquivo. Ele deve ser um arquivo de texto simples (`.txt`), onde as colunas são separadas **exclusivamente pela tecla TAB** (tabulação). Não deve ter cabeçalhos (títulos de colunas) nem espaços em branco soltos.

A estrutura mais básica tem duas colunas: `[Sequência_do_Barcode]<tab>[Nome_da_Amostra]`
<p align="center">
  <img width="223" height="222" alt="image" src="https://github.com/user-attachments/assets/427204af-e5fb-4f8d-afae-a6930986df85" />
  <br>
  <sub><b>Figura 2.</b> Formato do arquivo de barcodes para STACKS.</sub>
</p>

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
Já estamos na pasta `demultiplex`, onde criamos o arquivo `barcodes.txt`. Se você saiu dela, volte com:

```bash
cd ~/demultiplex
```

Agora vamos rodar o `process_radtags`. Os dados brutos estão em uma pasta compartilhada, usada por toda a turma, então indicamos o caminho completo até eles:

```bash
process_radtags -f /home/genpop/dados_brutos/gbs_acrocomia.fastq.gz \
                -b ~/demultiplex/barcodes.txt \
                -o ~/demultiplex \
                -e pstI \
                --adapter-1 AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG \
                --adapter-mm 2 \
                -r -c -q
```

O que faz cada opção:

| Opção | Função |
|---|---|
| `-f` | Arquivo FASTQ bruto (*single-end*), ainda sem demultiplexar |
| `-b` | Arquivo com os barcodes e os nomes das amostras |
| `-o` | Pasta de saída |
| `-e` | Enzima de restrição usada na biblioteca: o programa confere se o sítio de corte aparece logo depois do barcode |
| `--adapter-1` | Sequência do adaptador: reads que contêm o adaptador são descartados |
| `--adapter-mm` | Número de diferenças (*mismatches*) permitidas ao procurar o adaptador |
| `-r` | Recupera barcodes e sítios de restrição com erros de sequenciamento |
| `-c` | Remove reads que contêm bases não identificadas (`N`) |
| `-q` | Descarta reads de baixa qualidade |

### Conferir os resultados

```bash
ls -lh ~/demultiplex
```

Deve aparecer um arquivo `.fq.gz` para cada uma das 10 amostras (por exemplo, `AM_04.fq.gz`), mais o relatório `process_radtags.log`. Vamos olhar o resumo:

```bash
head -n 20 ~/demultiplex/process_radtags.log
```

> [!NOTE]
> Como usamos só 10 das 96 amostras do lane, a maior parte dos reads vai aparecer como *barcode não encontrado*. Isso é esperado: esses reads pertencem às outras 86 amostras.

💡 **Single-end sequencing versus paired-end sequencing**
> Em nosso exemplo, usamos uma biblioteca sequenciada single-end. Caso você tenha uma biblioteca paired-end, o comando a seguir ajuda a indicar onde se encontram as leituras R1 e R2. Os demais parâmetros são semelhantes.
>
```bash
process_radtags -1 ../data/raw/sample_R1.fastq -2 ../../data/raw/sample_R2.fastq
 ```

## 2. Controle de Qualidade Inicial (QC)
Agora que temos nossos arquivos FASTQ separados por amostra, precisamos fazer a avaliação dos arquivos para garantir dados de alta qualidade para análises posteriores. 

Lembre-se de que os arquivos são arquivos FASTQ, que são arquivos FASTA modificados para incluir também informações sobre a qualidade do sequenciamento para cada base (PHRED score). 
<p align="center">
 <img width="553" height="311" alt="fastqc" src="https://github.com/user-attachments/assets/880b65b5-b995-475e-8a04-229565039229" />
  <br>
  <sub><b>Figura 3.</b> Formato do arquivo FASTQC.</sub>
</p>

A qualidade Q reflete a probabilidade p de que a base esteja incorreta. Os valores de Q podem variar de 0 a 93 (mas a pontuação máxima que você normalmente verá é 40). Para economizar espaço em um arquivo FASTQ, cada pontuação de qualidade é transformada em um único caractere ASCII.
<p align="center">
 <img width="877" height="267" alt="5d0a8f1b-4a1d-4777-b8c9-a4325d242560" src="https://github.com/user-attachments/assets/8d1eb519-f2a2-481c-90bc-fbfa91fa0b23" />
  <br>
  <sub><b>Figura 4.</b> Interpretação dos caracteres ASCII em valores de qualidade das sequencias (Q).</sub>
</p>

### 0. Visualizando os arquivos FASTQ

Antes de gerar os relatórios de qualidade, vamos usar comandos básicos do UNIX para explorar os dados e nos familiarizar com o formato FASTQ. Os arquivos demultiplexados estão na pasta `~/demultiplex`, com a extensão `.fq.gz` (o `.gz` indica que estão comprimidos).

**A.** Como os arquivos estão comprimidos, usamos o `zcat` para descomprimi-los "no voo", sem criar um arquivo novo, e o `head` para ver só o começo:

```bash
zcat ~/demultiplex/AM_04.fq.gz | head
```

> [!TIP]
> O símbolo `|` é o *pipe*. Ele pega a saída do comando à esquerda e a entrega como entrada para o comando à direita, sem salvar nenhum arquivo intermediário. Assim, encadeamos várias operações em uma única linha.

❓ Cada read ocupa quantas linhas? O que tem em cada uma delas? Qual é o caractere inicial do cabeçalho de todos os reads?

**B.** Assim como o `head` lê o começo, o `tail` lê o final do arquivo. Isso é útil para conferir se um arquivo está completo: um arquivo cortado no meio de uma transferência termina com um read incompleto.

```bash
zcat ~/demultiplex/AM_04.fq.gz | tail
zcat ~/demultiplex/AM_04.fq.gz | tail -n 16
```

> [!TIP]
> Por padrão, `head` e `tail` mostram 10 linhas. A opção `-n` define quantas linhas queremos ver.

**C.** Quantos reads tem cada amostra? Como cada read ocupa exatamente 4 linhas, basta contar as linhas do arquivo e dividir por 4:

```bash
echo $(( $(zcat ~/demultiplex/AM_04.fq.gz | wc -l) / 4 ))
```

E para todas as amostras de uma vez, usando um *loop*:

```bash
for arquivo in ~/demultiplex/*.fq.gz; do
    echo "$(basename $arquivo): $(( $(zcat $arquivo | wc -l) / 4 )) reads"
done
```

❓ Quantos reads tem cada indivíduo? Os números batem com o que aparece no `process_radtags.log`? Há amostras com muito menos reads que as outras?

❓ Por que não podemos simplesmente contar quantas linhas começam com `@`? (Dica: olhe com atenção a quarta linha de cada read.)

### 1. Avaliando a qualidade dos dados com o FastQC

O [FastQC](https://github.com/s-andrews/FastQC) é uma ferramenta amplamente usada para avaliar a qualidade de dados de sequenciamento. Ele gera um relatório HTML com métricas como:

1. **Qualidade por base:** a qualidade média em cada posição ao longo dos reads.
2. **Conteúdo GC:** a porcentagem de bases G e C.
3. **Níveis de duplicação:** identifica possíveis duplicatas de PCR.
4. **Conteúdo de adaptador:** verifica se sobraram sequências de adaptador.

```bash
# Criar a pasta para os relatórios de qualidade dos dados brutos
mkdir -p ~/qc/brutos

# Rodar o FastQC em todos os arquivos demultiplexados
fastqc -t 1 ~/demultiplex/*.fq.gz --outdir ~/qc/brutos
```

> [!TIP]
> - O `*` seleciona **todos** os arquivos que terminam em `.fq.gz` dentro da pasta.
> - `-o` indica a pasta de saída dos relatórios.
> - `-t` define quantos núcleos (*threads*) usar. Neste servidor, cada um de nós tem 1 núcleo disponível.

Para abrir os relatórios HTML, é preciso copiá-los para o seu computador (veja a seção [Copiando arquivos do servidor para o seu computador](#copiando-arquivos-do-servidor-para-o-seu-computador)).

❓ Avalie a qualidade dos arquivos. Qual é o comprimento dos reads? A qualidade cai ao longo da sequência? Qual amostra tem a melhor e qual tem a pior qualidade? Foi detectada alguma sequência de adaptador?

### 2. Um relatório único de todas as amostras com o MultiQC

O FastQC gera um relatório separado para cada arquivo. Com muitas amostras (o que é o normal), abrir e comparar um por um fica inviável. O [MultiQC](https://github.com/MultiQC/MultiQC) junta todos em um único relatório:

```bash
multiqc ~/qc/brutos -o ~/qc/brutos
```

Agora copie o arquivo `multiqc_report.html` para o seu computador, abra no navegador e vamos avaliar todas as amostras juntas.

❓ Que características de uma biblioteca GBS poderiam explicar esses resultados?

---

## Copiando arquivos do servidor para o seu computador

O servidor não tem tela nem navegador, então, para ver os relatórios HTML, precisamos copiá-los para o nosso computador. Para isso usamos o `scp` (*secure copy*), que funciona como o `cp`, mas entre computadores diferentes, usando a mesma conexão do SSH.

> [!IMPORTANT]
> O `scp` é rodado **no terminal do seu computador**, e não no servidor. Abra um terminal novo:
> - **Windows:** PowerShell (o `scp` já vem instalado no Windows 10 e 11)
> - **Mac e Linux:** Terminal

Primeiro, vá para a pasta onde você quer salvar os arquivos, por exemplo Downloads:

```bash
cd Downloads
```

O formato do comando é `scp -P PORTA usuario@servidor:caminho/no/servidor destino/no/computador`. Para copiar o relatório do MultiQC:

```bash
scp -P 59022 USUARIO@dayhoff1.esalq.usp.br:~/qc/brutos/multiqc_report.html .
```

Para copiar vários arquivos de uma vez, por exemplo todos os relatórios do FastQC, use o `*` **entre aspas**, para que ele seja interpretado no servidor e não no seu computador:

```bash
scp -P 59022 "USUARIO@dayhoff1.esalq.usp.br:~/qc/brutos/multiqc_report.html" .
```

E para copiar uma pasta inteira, com tudo o que tem dentro, use `-r`:

```bash
scp -P 59022 -r USUARIO@dayhoff1.esalq.usp.br:~/qc/brutos .
```

> [!TIP]
> - O `.` no final significa "a pasta atual", que no exemplo é Downloads.
> - No `scp` a porta se indica com **`-P` maiúsculo**, enquanto no `ssh` é com `-p` minúsculo. É um erro muito comum!
> - O `scp` vai pedir a sua senha, a mesma do SSH.

