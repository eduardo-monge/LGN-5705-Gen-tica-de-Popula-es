# 💻 2: Trimagem (*trimming*) das sequencias

Nesta prática, vamos aprender a analisar a fazer a trimagem das sequencias para deixar só os dados de boa qualidade, sem adaptador e sem cauda poli-G. 

---
## 1. Limpeza dos dados (*trimming*) com fastp
O [fastp](https://github.com/OpenGene/fastp) é uma ferramenta rápida de controle de qualidade e pré-processamento. Em um único comando, ele:

1. Corta bases de baixa qualidade.
2. Remove adaptadores.
3. Remove caudas poli-G e poli-X.
4. Filtra reads curtos demais.
5. Gera um relatório de qualidade.

Na vida real, você terá dezenas ou centenas de amostras, e rodar o comando manualmente para cada uma seria inviável. Por isso, vamos usar desde o início um **loop `for`**, que percorre todos os arquivos e aplica o fastp em cada um deles automaticamente.

```bash
# Criar as pastas de saída: uma para os reads limpos e outra para os relatórios do fastp
mkdir -p ~/trimmed ~/qc/fastp

# Rodar o fastp em todas as amostras
for arquivo in ~/demultiplex/*.fq.gz; do
    amostra=$(basename "$arquivo" .fq.gz)
    echo "Processando a amostra: $amostra"

    fastp -i "$arquivo" \
          -o ~/trimmed/${amostra}_trimmed.fq.gz \
          -w 1 \
          -q 20 \
          --cut_tail \
          --cut_mean_quality 20 \
          -l 30 \
          --trim_poly_g \
          --trim_poly_x \
          -h ~/qc/fastp/${amostra}_fastp.html \
          -j ~/qc/fastp/${amostra}_fastp.json
done
```

> [!TIP]
> **Como funciona o loop:**
> - `for arquivo in ~/demultiplex/*.fq.gz; do`: "para cada arquivo que termina em `.fq.gz`, guarde o nome dele na variável `arquivo` e execute os comandos abaixo".
> - `basename "$arquivo" .fq.gz`: tira o caminho e a extensão, deixando só o nome da amostra (por exemplo, `AM_04`). Assim criamos nomes de saída organizados.
> - `$arquivo` e `${amostra}`: o cifrão chama o conteúdo da variável.
> - `done`: marca o fim dos comandos. O loop passa para o próximo arquivo da lista.
>
> **Os parâmetros do fastp:**
> - `-i` / `-o`: arquivo de entrada e de saída.
> - `-w 1`: número de núcleos (*threads*).
> - `-q 20`: qualidade mínima para uma base ser considerada "boa". Reads com mais de 40% de bases abaixo de Q20 são descartados.
> - `--cut_tail` e `--cut_mean_quality 20`: percorre o read a partir do final e corta enquanto a qualidade média da janela estiver abaixo de Q20.
> - `-l 30`: descarta reads com menos de 30 pb depois dos cortes.
> - `--trim_poly_g` / `--trim_poly_x`: remove caudas de G (ou de qualquer base repetida) no final dos reads.
> - `-h` / `-j`: relatório em HTML (para nós) e em JSON (para o MultiQC).

> [!WARNING]
> Nos comandos com várias linhas, a barra invertida `\` precisa ser o **último caractere** da linha. Um espaço depois dela quebra o comando!

### Conferindo a limpeza com o FastQC e o MultiQC

Agora precisamos confirmar se a limpeza funcionou. Vamos gerar novos relatórios e salvá-los em uma pasta separada, para não misturar com os dos dados brutos:

```bash
mkdir -p ~/qc/trimmed

# FastQC nos reads limpos
fastqc -t 1 ~/trimmed/*_trimmed.fq.gz --outdir ~/qc/trimmed

# MultiQC juntando os relatórios do FastQC e do fastp
multiqc ~/qc/trimmed ~/qc/fastp -o ~/qc/trimmed
```

Copie para o seu computador os arquivos `~/qc/trimmed/multiqc_report.html` e `~/qc/comparacao/multiqc_report.html` (veja a seção [Copiando arquivos do servidor para o seu computador](#copiando-arquivos-do-servidor-para-o-seu-computador)).

❓ Abra um dos relatórios HTML do fastp (em `~/qc/fastp`). Quantos reads havia antes e depois da filtragem? Qual foi a porcentagem de reads descartados por baixa qualidade, por excesso de bases N e por serem curtos demais?

❓ No relatório de comparação, o que mudou entre os dados brutos e os limpos? Observe:

1. **Qualidade por base:** a qualidade deve estar agora predominantemente na zona verde, inclusive no final dos reads.
2. **Comprimento das sequências:** depois dos cortes, os reads não têm mais todos o mesmo tamanho.
3. **Caudas poli-G:** se existiam, o conteúdo de G no final dos reads deve ter caído.
4. **Conteúdo de adaptador:** quase não deve haver adaptadores. Por quê? (Dica: lembre-se das opções que usamos no `process_radtags`.)

---

> [!TIP]
> **Boas práticas no controle de qualidade:**
> 1. **Guarde os relatórios:** mantenha os resultados do FastQC e do fastp, para que a análise seja reprodutível.
> 2. **Ajuste os parâmetros:** escolha os limites de corte com base no que os relatórios mostram, e não só nos valores padrão.
> 3. **Automatize:** use loops (ou *job arrays*, em servidores com fila) para muitas amostras.
> 4. **Organize as pastas:** separe dados brutos, dados limpos e relatórios de qualidade.
