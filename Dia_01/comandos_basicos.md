# 💻 Prática: Introdução ao Linux para Bioinformática

Nesta prática, vamos aprender os comandos essenciais do terminal Linux para navegar, manipular arquivos e explorar dados — habilidades fundamentais para qualquer pipeline de bioinformática.

---

## 1. Navegação e Estrutura de Diretórios

Antes de qualquer análise, é preciso saber onde você está e como se mover pelo sistema de arquivos.

### `ls` — Listar o conteúdo de um diretório

```bash
# Listar arquivos e pastas no diretório atual
ls

# Listar com detalhes (permissões, tamanho, data)
ls -l

# Listar com tamanho legível para humanos (KB, MB, GB)
ls -lh

# Listar recursivamente (incluindo subdiretórios)
ls -R

# Listar arquivos ocultos (que começam com .)
ls -la
```

### `cd` — Mudar de diretório

```bash
# Ir para o diretório home do usuário
cd ~

# Entrar em uma pasta
cd data/raw

# Voltar um nível
cd ..

# Voltar dois níveis
cd ../..

# Ir para o último diretório visitado
cd -
```

### `clear` — Limpar a tela

```bash
# Limpar o terminal (os comandos anteriores continuam no histórico)
clear

# Atalho equivalente
Ctrl + L
```

---

## 2. Criação e Manipulação de Arquivos e Diretórios

### `mkdir` — Criar diretórios

```bash
# Criar um único diretório
mkdir resultados

# Criar uma hierarquia completa de uma vez (opção -p)
mkdir -p data/raw data/trimmed data/mapped
mkdir -p results/qc results/snp_calling

# -----------------------------------------
# A opção -p significa "pais" (parents).
# Ela instrui o sistema a:
# 1️⃣ Criar todos os diretórios pai, caso não existam
# 2️⃣ Evitar erro se o diretório já existir
# -----------------------------------------
```

### `touch` — Criar um arquivo vazio

```bash
# Criar um arquivo vazio
touch metadata.tsv

# Criar vários arquivos de uma vez
touch arquivo1.txt arquivo2.txt arquivo3.txt
```

### `cp` — Copiar arquivos ou diretórios

```bash
# Copiar um arquivo para outro diretório
cp amostra.fastq.gz data/raw/

# Copiar e renomear ao mesmo tempo
cp amostra.fastq.gz data/raw/amostra_backup.fastq.gz

# Copiar um diretório inteiro (opção -r, recursivo)
cp -r scripts/ scripts_backup/
```

### `mv` — Mover ou renomear arquivos

```bash
# Mover um arquivo para outro diretório
mv amostra.fastq.gz data/raw/

# Renomear um arquivo (é um "mover" para o mesmo lugar com outro nome)
mv metadata.txt metadata.tsv

# Mover vários arquivos de uma vez
mv *.fastq.gz data/raw/
```

### `rm` — Remover arquivos ou diretórios

```bash
# Remover um arquivo
rm arquivo_temporario.txt

# Remover um diretório e todo o seu conteúdo (CUIDADO: irreversível!)
rm -r pasta_temporaria/

# Perguntar antes de remover (opção -i, interativo)
rm -i arquivo_importante.txt
```

> ⚠️ **Atenção:** No Linux, **não existe lixeira**. O `rm` apaga definitivamente. Use com cautela, especialmente com `-r` e curingas (`*`).

---

## 3. Leitura e Edição de Arquivos

### `nano` — Editor de texto no terminal

```bash
# Abrir ou criar um arquivo para edição
nano metadata.tsv
```

Dentro do editor, os atalhos mais importantes são:

| Atalho     | Ação                        |
|------------|-----------------------------|
| `Ctrl + O` | Salvar (Write Out)          |
| `Enter`    | Confirmar o nome do arquivo |
| `Ctrl + X` | Sair do editor              |
| `Ctrl + W` | Buscar texto                |
| `Ctrl + K` | Recortar linha              |
| `Ctrl + U` | Colar linha                 |

### `cat` — Exibir o conteúdo completo de um arquivo

```bash
# Mostrar todo o conteúdo de um arquivo
cat metadata.tsv

# Concatenar dois arquivos e mostrar na tela
cat arquivo1.txt arquivo2.txt

# Concatenar dois arquivos e salvar o resultado em um terceiro
cat arquivo1.txt arquivo2.txt > arquivo_unido.txt
```

### `head` — Ver as primeiras linhas de um arquivo

```bash
# Mostrar as 10 primeiras linhas (padrão)
head metadata.tsv

# Mostrar as 4 primeiras linhas (útil para ver um read FASTQ completo)
head -n 4 amostra.fastq

# Ver o primeiro read de um arquivo compactado
zcat amostra.fq.gz | head -n 4
```

### `tail` — Ver as últimas linhas de um arquivo

```bash
# Mostrar as 10 últimas linhas (padrão)
tail metadata.tsv

# Mostrar as últimas 20 linhas
tail -n 20 log_pipeline.txt

# Monitorar um arquivo em tempo real (útil para acompanhar logs)
tail -f log_pipeline.txt
```

### `less` — Navegar pelo conteúdo de arquivos grandes

```bash
# Abrir um arquivo para navegação interativa
less resultados.txt

# Abrir um arquivo compactado diretamente
zcat amostra.fq.gz | less
```

Dentro do `less`, os atalhos são:

| Atalho      | Ação                     |
|-------------|--------------------------|
| `Espaço`    | Próxima página           |
| `b`         | Página anterior          |
| `/padrão`   | Buscar texto             |
| `n`         | Próxima ocorrência       |
| `q`         | Sair                     |

---

## 4. Compressão e Arquivos Compactados

Arquivos de sequenciamento são grandes. Em bioinformática, quase sempre trabalhamos com dados comprimidos (`.gz`).

### `gzip` / `zip` — Comprimir arquivos

```bash
# Comprimir um arquivo (gera arquivo.fastq.gz e remove o original)
gzip amostra.fastq

# Comprimir mantendo o arquivo original (opção -k, keep)
gzip -k amostra.fastq

# Comprimir múltiplos arquivos
gzip *.fastq
```

### `gzip -d` / `unzip` — Descomprimir arquivos

```bash
# Descomprimir um arquivo .gz
gzip -d amostra.fastq.gz

# Alternativa para .gz
gunzip amostra.fastq.gz

# Descomprimir um arquivo .zip
unzip arquivo.zip

# Descomprimir em uma pasta específica
unzip arquivo.zip -d pasta_destino/
```

### `zcat` — Ler arquivos comprimidos sem descomprimir

```bash
# Ver o conteúdo de um .gz direto no terminal
zcat amostra.fq.gz

# Combinar com outros comandos via pipe (|)
zcat amostra.fq.gz | head -n 8
zcat amostra.fq.gz | wc -l
```

> 💡 **Dica:** Em bioinformática, quase nunca descomprimimos os arquivos FASTQ. Ferramentas como `FastQC`, `Trimmomatic` e `BWA` leem `.gz` diretamente, economizando disco.

---

## 5. Busca e Filtragem de Dados

### `grep` — Buscar padrões em arquivos

```bash
# Buscar uma palavra em um arquivo
grep "amazonas" metadata.tsv

# Contar quantas linhas contêm o padrão (-c)
grep -c "^@" amostra.fastq

# Busca sem diferenciar maiúsculas/minúsculas (-i)
grep -i "Amazonas" metadata.tsv

# Mostrar o número da linha (-n)
grep -n "ERROR" log_pipeline.txt

# Inverter: mostrar linhas que NÃO contêm o padrão (-v)
grep -v "#" metadata.tsv

# Buscar em arquivos comprimidos
zcat amostra.fq.gz | grep "^@" | head -n 5
```

### `find` — Procurar arquivos no sistema

```bash
# Encontrar um arquivo pelo nome no diretório atual e subpastas
find . -name "metadata.tsv"

# Encontrar todos os arquivos .fastq.gz
find . -name "*.fastq.gz"

# Encontrar apenas diretórios com um nome específico
find . -type d -name "raw"

# Encontrar arquivos modificados nas últimas 24h
find . -mtime -1
```

### `wc` — Contar linhas, palavras e caracteres

```bash
# Contar linhas de um arquivo
wc -l metadata.tsv

# Contar o número de reads em um FASTQ (dividir por 4)
zcat amostra.fq.gz | wc -l

# Contar palavras
wc -w relatorio.txt

# Contar caracteres
wc -c amostra.fastq
```

---

## 6. Manipulação de Texto

### `sort` — Ordenar o conteúdo de um arquivo

```bash
# Ordenar alfabeticamente
sort lista.txt

# Ordenar em ordem inversa
sort -r lista.txt

# Ordenar numericamente (opção -n)
sort -n contagens.txt

# Ordenar por uma coluna específica (coluna 2)
sort -k2 metadata.tsv
```

### `uniq` — Remover ou contar duplicatas

```bash
# Remover linhas duplicadas consecutivas (o arquivo precisa estar ordenado)
sort lista.txt | uniq

# Contar quantas vezes cada linha aparece (-c)
sort populacoes.txt | uniq -c

# Mostrar apenas as linhas duplicadas (-d)
sort lista.txt | uniq -d

# Mostrar apenas as linhas únicas (-u)
sort lista.txt | uniq -u
```

### `cut` — Extrair colunas de um arquivo

```bash
# Extrair a primeira coluna (delimitador tab, padrão)
cut -f1 metadata.tsv

# Extrair a segunda coluna com delimitador vírgula
cut -d',' -f2 arquivo.csv

# Extrair colunas 1 e 3 com delimitador ponto-e-vírgula
cut -d';' -f1,3 dados.txt

# Extrair os primeiros 10 caracteres de cada linha
cut -c1-10 arquivo.txt
```

### `sed` — Substituir texto em um arquivo

```bash
# Substituir a primeira ocorrência em cada linha
sed 's/oldtext/newtext/' arquivo.txt

# Substituir todas as ocorrências em cada linha (opção g, global)
sed 's/oldtext/newtext/g' arquivo.txt

# Salvar as alterações no mesmo arquivo (opção -i)
sed -i 's/oldtext/newtext/g' arquivo.txt

# Remover linhas que contêm um padrão
sed '/^#/d' arquivo.txt

# Substituir tabulações por vírgulas (útil para converter TSV em CSV)
sed 's/\t/,/g' metadata.tsv > metadata.csv
```

### `awk` — Processar e transformar texto estruturado

```bash
# Imprimir a segunda coluna de um arquivo
awk '{print $2}' metadata.tsv

# Imprimir colunas 1 e 2 com delimitador tab
awk '{print $1"\t"$2}' metadata.tsv

# Filtrar linhas onde a segunda coluna é "amazonas"
awk '$2 == "amazonas"' metadata.tsv

# Calcular a soma de uma coluna numérica
awk '{sum += $1} END {print sum}' contagens.txt

# Contar o número de reads em um FASTQ (1 read = 4 linhas)
awk 'NR % 4 == 1' amostra.fastq | wc -l
```

### `cmp` e `diff` — Comparar arquivos

```bash
# Verificar se dois arquivos são idênticos (cmp)
cmp arquivo1.txt arquivo2.txt

# Mostrar as diferenças linha a linha (diff)
diff arquivo1.txt arquivo2.txt

# Formato unificado, mais legível (opção -u)
diff -u arquivo1.txt arquivo2.txt

# Comparar dois diretórios recursivamente
diff -r pasta1/ pasta2/
```

---

## 7. Permissões e Administração

### `sudo` — Executar como superusuário

```bash
# Executar um comando com permissão de administrador
sudo apt update

# Abrir um arquivo protegido para edição
sudo nano /etc/hosts

# -----------------------------------------
# Use sudo com responsabilidade.
# Comandos executados como root afetam
# o sistema inteiro e não pedem confirmação.
# -----------------------------------------
```

### `chmod` — Alterar permissões de arquivos

```bash
# Tornar um script executável
chmod +x meu_script.sh

# Remover permissão de escrita (proteger um arquivo de edições acidentais)
chmod -w dados_brutos.fastq.gz

# Definir permissões completas para o dono, leitura para outros
chmod 644 metadata.tsv

# Tornar um diretório acessível a todos
chmod 755 resultados/
```

> 💡 **Entendendo as permissões:** O comando `ls -l` mostra a permissão de cada arquivo no formato `-rwxr-xr--`, onde:
> - `r` = leitura (read)
> - `w` = escrita (write)  
> - `x` = execução (execute)

### `apt` — Instalar e gerenciar programas

```bash
# Atualizar a lista de pacotes disponíveis
sudo apt update

# Instalar um programa
sudo apt install fastqc

# Instalar múltiplos programas de uma vez
sudo apt install bwa samtools bcftools

# Remover um programa
sudo apt remove fastqc

# Atualizar todos os programas instalados
sudo apt upgrade
```

### `man` — Acessar o manual de um comando

```bash
# Ver o manual completo do comando ls
man ls

# Ver o manual do grep
man grep

# -----------------------------------------
# Dentro do man:
# Espaço = próxima página
# b      = página anterior
# /termo = buscar
# q      = sair
# -----------------------------------------
```

### `history` — Ver o histórico de comandos

```bash
# Listar todos os comandos anteriores
history

# Ver os últimos 20 comandos
history 20

# Buscar um comando específico no histórico
history | grep "fastqc"

# Reexecutar o comando número 42 do histórico
!42

# Reexecutar o último comando que começa com "bwa"
!bwa
```

---

## 8. Monitoramento do Sistema

### `free` — Ver memória disponível

```bash
# Mostrar uso de memória RAM e swap
free

# Mostrar em formato legível (MB/GB)
free -h

# Atualizar a cada 2 segundos
watch -n 2 free -h
```

### `htop` — Monitorar processos em tempo real

```bash
# Abrir o monitor de processos (interface interativa)
htop

# -----------------------------------------
# Dentro do htop:
# F6     = ordenar por coluna
# F9     = matar processo
# F10/q  = sair
# -----------------------------------------
```

### `kill` — Encerrar processos

```bash
# Encerrar um processo pelo seu ID (PID)
kill 1234

# Forçar encerramento imediato (sinal -9)
kill -9 1234

# Encontrar o PID de um programa pelo nome
pgrep fastqc

# Matar todos os processos com um nome
pkill fastqc
```

---

## 9. Desligar e Reiniciar

### `reboot` e `shutdown`

```bash
# Reiniciar imediatamente
sudo reboot

# Desligar imediatamente
sudo shutdown -h now

# Desligar em 10 minutos
sudo shutdown -h +10

# Agendar reinicialização para horário específico
sudo shutdown -r 22:00

# Cancelar um shutdown agendado
sudo shutdown -c
```

---

## 📌 Referência Rápida

| Comando         | Função principal                                    |
|-----------------|-----------------------------------------------------|
| `ls`            | Listar conteúdo de diretório                        |
| `cd`            | Navegar entre diretórios                            |
| `clear`         | Limpar a tela do terminal                           |
| `mkdir`         | Criar diretórios                                    |
| `touch`         | Criar arquivo vazio                                 |
| `cp`            | Copiar arquivos ou diretórios                       |
| `mv`            | Mover ou renomear arquivos                          |
| `rm`            | Remover arquivos ou diretórios                      |
| `nano`          | Editar arquivos no terminal                         |
| `cat`           | Exibir conteúdo completo de um arquivo              |
| `head`          | Ver as primeiras N linhas                           |
| `tail`          | Ver as últimas N linhas                             |
| `less`          | Navegar por arquivos grandes                        |
| `gzip` / `zip`  | Comprimir arquivos                                  |
| `gzip -d` / `unzip` | Descomprimir arquivos                          |
| `zcat`          | Ler arquivos `.gz` sem descomprimir                 |
| `grep`          | Buscar padrões em arquivos                          |
| `find`          | Pesquisar arquivos no sistema                       |
| `wc`            | Contar linhas, palavras e caracteres                |
| `sort`          | Ordenar conteúdo                                    |
| `uniq`          | Remover ou contar duplicatas                        |
| `cut`           | Extrair colunas por delimitador                     |
| `sed`           | Substituir texto                                    |
| `awk`           | Processar texto estruturado por colunas             |
| `cmp`           | Verificar se dois arquivos são idênticos            |
| `diff`          | Mostrar diferenças entre arquivos                   |
| `sudo`          | Executar com permissão de superusuário              |
| `chmod`         | Alterar permissões de leitura/escrita/execução      |
| `apt`           | Instalar e gerenciar programas                      |
| `man`           | Acessar o manual de um comando                      |
| `history`       | Ver histórico de comandos                           |
| `free`          | Ver uso de memória RAM                              |
| `htop`          | Monitorar processos em tempo real                   |
| `kill`          | Encerrar um processo pelo PID                       |
| `reboot`        | Reiniciar a máquina                                 |
| `shutdown`      | Desligar a máquina                                  |

