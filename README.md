# 🧬 NGS Bioinformatics Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/username/ngs-bioinformatics-pipeline.svg)](https://github.com/username/ngs-bioinformatics-pipeline/stargazers)
[![GitHub issues](https://img.shields.io/github/issues/username/ngs-bioinformatics-pipeline.svg)](https://github.com/username/ngs-bioinformatics-pipeline/issues)

> **Полноценный конвейер для анализа данных NGS**: аннотация генома, вызов вариантов (SNV/indel/SV), филогенетический анализ и интерпретация.

---

## 📋 Содержание

- [🚀 Быстрый старт](#-быстрый-старт)
- [📚 Теоретическая база](#-теоретическая-база)
- [🔧 Пайплайны](#-пайплайны)
- [📊 Структура проекта](#-структура-проекта)
- [📦 Установка](#-установка)
- [🧪 Тестирование](#-тестирование)
- [📝 Лицензия](#-лицензия)

---

## 🚀 Быстрый старт

### 1. Клонирование репозитория

```bash
git clone https://github.com/username/ngs-bioinformatics-pipeline.git
cd ngs-bioinformatics-pipeline
```

### 2. Установка окружения

```bash
# Через conda (рекомендуется)
conda env create -f environment.yml
conda activate ngs-pipeline

# Или вручную: установить зависимости из requirements.txt
pip install -r requirements.txt
```

### 3. Запуск примера аннотации

```bash
bash pipelines/annotation/prokka_pipeline.sh \
  --input data/samples/example.fasta \
  --outdir results/annotation_test \
  --kingdom Bacteria
```

---

## 📚 Теоретическая база

Все теоретические материалы с формулами и объяснениями находятся в папке [`docs/`](docs/).

| Файл | Тема | Статус | Сложность | Теги |
|------|------|--------|-----------|------|
| [`📄 01_genome_annotation.md`](docs/01_genome_annotation.md) | Аннотация генома | ✅ Готово | 🟢 Базовый | `ORF` `Prokka` `Prodigal` `E-value` `KEGG` `GO` |
| [`📄 02_variant_calling.md`](docs/02_variant_calling.md) | Вызов вариантов | ✅ Готово | 🟡 Средний | `VCF` `QUAL` `PL` `GQ` `FreeBayes` `GATK` |
| [`📄 03_structural_variants.md`](docs/03_structural_variants.md) | Структурные варианты | ✅ Готово | 🔴 Продвинутый | `SVTYPE` `split-reads` `discordant-pairs` `CNV` |
| [`📄 04_variant_interpretation.md`](docs/04_variant_interpretation.md) | Интерпретация вариантов | ✅ Готово | 🔴 Продвинутый | `ACMG` `gnomAD` `ClinVar` `consequence` |
| [`📄 05_phylogenetics.md`](docs/05_phylogenetics.md) | Филогенетика | ✅ Готово | 🟡 Средний | `MSA` `NJ` `bootstrap` `ML` `Bayesian` |
| [`📐 formulas_reference.md`](docs/formulas_reference.md) | Сводка формул | ✅ Готово | 🟢 Базовый | `формулы` `справочник` |

> 💡 **Совет**: Начните с [`01_genome_annotation.md`](docs/01_genome_annotation.md) для знакомства с базовыми концепциями.

---

## 🔧 Пайплайны

### Доступные пайплайны

| Пайплайн | Описание | Входные данные | Выходные данные |
|----------|----------|----------------|-----------------|
| **🔬 Annotation** | Структурная и функциональная аннотация генома | `.fasta` (assembly) | `.gff`, `.faa`, `.gbk` |
| **🧬 Variant Calling** | Вызов SNV и indel вариантов | `.bam`, `.fasta` (ref) | `.vcf.gz` |
| **🔍 Structural Variants** | Детекция крупных структурных вариантов | `.bam`, `.fasta` (ref) | `.vcf.gz` (SV) |
| **🌳 Phylogenetics** | Построение филогенетических деревьев | `.fasta` (alignment) | `.tre`, `.nwk` |

### Запуск пайплайнов

#### Аннотация генома (Prokka)

```bash
bash pipelines/annotation/prokka_pipeline.sh \
  --input data/samples/genome.fasta \
  --outdir results/annotation \
  --kingdom Bacteria \
  --cpus 8
```

#### Вызов вариантов (FreeBayes)

```bash
bash pipelines/variant_calling/freebayes_pipeline.sh \
  --reference data/reference/ref.fasta \
  --bam sample.bam \
  --output results/variants/sample.vcf \
  --min-alternate-fraction 0.2
```

#### Филогенетический анализ (IQ-TREE)

```bash
bash pipelines/phylogenetics/iqtree_analysis.sh \
  --alignment data/samples/alignment.fasta \
  --outdir results/phylogeny \
  --bootstrap 1000 \
  --threads 4
```

> 📖 Подробная документация по каждому пайплайну находится в соответствующих подпапках [`pipelines/`](pipelines/).

---

## 📊 Структура проекта

```
ngs-bioinformatics-pipeline/
│
├── 📄 README.md                    # Главная страница проекта
├── 📄 LICENSE                      # Лицензия MIT
├── 📄 CITATION.cff                 # Как цитировать проект
├── 📄 .gitignore                   # Исключаемые файлы
├── 📄 requirements.txt             # Python-зависимости
├── 📄 environment.yml              # Conda-окружение
│
├── 📂 docs/                        # 📚 Теория и документация
│   ├── 📄 01_genome_annotation.md
│   ├── 📄 02_variant_calling.md
│   ├── 📄 03_structural_variants.md
│   ├── 📄 04_variant_interpretation.md
│   ├── 📄 05_phylogenetics.md
│   └── 📄 formulas_reference.md
│
├── 📂 pipelines/                   # 🔧 Готовые пайплайны
│   ├── 📂 annotation/
│   │   ├── 📄 prokka_pipeline.sh
│   │   ├── 📄 functional_annotation.sh
│   │   └── 📄 README.md
│   ├── 📂 variant_calling/
│   │   ├── 📄 freebayes_pipeline.sh
│   │   ├── 📄 bcftools_filtering.sh
│   │   ├── 📄 gatk_haplotypecaller.sh
│   │   └── 📄 README.md
│   └── 📂 phylogenetics/
│       ├── 📄 mafft_alignment.sh
│       ├── 📄 iqtree_analysis.sh
│       └── 📄 README.md
│
├── 📂 scripts/                     # 🛠️ Вспомогательные утилиты
│   ├── 📄 parse_vcf_stats.py
│   ├── 📄 calculate_n50.py
│   └── 📄 utils.sh
│
├── 📂 data/                        # 📦 Примеры данных
│   ├── 📂 reference/
│   │   └── 📄 example_ref.fasta
│   └── 📂 samples/
│       └── 📄 example.fasta
│
├── 📂 results/                     # 📊 Результаты (.gitignore!)
│   └── 📄 .gitkeep
│
└── 📂 tests/                       # 🧪 Автотесты
    ├── 📄 test_annotation.sh
    ├── 📄 test_variant_calling.sh
    └── 📂 test_data/
```

---

## 📦 Установка

### Системные требования

- **ОС**: Linux (Ubuntu 18.04+, CentOS 7+) или macOS
- **RAM**: минимум 8 GB (рекомендуется 16+ GB для WGS)
- **CPU**: 4+ ядер (рекомендуется 8+ для параллельных вычислений)
- **Диск**: 100+ GB свободного пространства

### Зависимости

#### Основные инструменты

| Инструмент | Версия | Назначение |
|------------|--------|------------|
| **Prokka** | ≥1.14 | Аннотация прокариот |
| **Prodigal** | ≥2.6 | Предсказание генов |
| **FreeBayes** | ≥1.3 | Вызов вариантов |
| **bcftools** | ≥1.9 | Работа с VCF/BCF |
| **GATK** | ≥4.0 | Вызов вариантов (eukaryotes) |
| **MAFFT** | ≥7.4 | Множественное выравнивание |
| **IQ-TREE** | ≥2.0 | Филогенетический анализ |
| **SAMtools** | ≥1.10 | Работа с BAM |
| **BWA** | ≥0.7.17 | Выравнивание чтений |

#### Python-библиотеки

```bash
pip install -r requirements.txt
```

**requirements.txt:**
```
pandas>=1.3.0
numpy>=1.20.0
biopython>=1.79
matplotlib>=3.4.0
seaborn>=0.11.0
scipy>=1.7.0
```

### Установка через Conda (рекомендуется)

```bash
conda env create -f environment.yml
conda activate ngs-pipeline
```

### Ручная установка

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y prokka samtools bcftools bedtools

# macOS (через Homebrew)
brew install prokka samtools bcftools bedtools
```

---

## 🧪 Тестирование

### Запуск тестов

```bash
# Все тесты
bash tests/run_all_tests.sh

# Отдельные тесты
bash tests/test_annotation.sh
bash tests/test_variant_calling.sh
bash tests/test_phylogenetics.sh
```

### Проверка результатов

```bash
# Проверка статистики аннотации
python scripts/calculate_n50.py results/annotation/example.gff

# Анализ VCF
bcftools stats results/variants/sample.vcf.gz
```

---

## 📝 Примеры использования

### Полный пайплайн анализа бактериального генома

```bash
#!/bin/bash

# 1. Аннотация
prokka --outdir annotation --prefix genome --kingdom Bacteria genome.fasta

# 2. Выравнивание чтений на референс
bwa index reference.fasta
bwa mem -t 8 reference.fasta reads_R1.fastq.gz reads_R2.fastq.gz | \
  samtools sort -@4 -o aligned.bam
samtools index aligned.bam

# 3. Вызов вариантов
freebayes -f reference.fasta -b aligned.bam --min-alternate-fraction 0.2 > variants.vcf

# 4. Фильтрация
bcftools filter -i 'QUAL>30 && DP>10' variants.vcf -Oz -o filtered.vcf.gz

# 5. Аннотация вариантов
bcftools csq -f reference.fasta -g genes.gff3 filtered.vcf.gz -Ov -o annotated.vcf

# 6. Построение дерева (если несколько образцов)
mafft --auto sequences.fasta > alignment.fasta
iqtree3 -s alignment.fasta -m MFP -B 1000 -T 4
```

---

## 🤝 Вклад в проект

Приветствуются pull requests, issues и предложения по улучшению!

### Как внести вклад

1. Fork репозиторий
2. Создайте ветку (`git checkout -b feature/AmazingFeature`)
3. Commit изменения (`git commit -m 'Add some AmazingFeature'`)
4. Push в ветку (`git push origin feature/AmazingFeature`)
5. Откройте Pull Request

---

## 📚 Дополнительные ресурсы

### Базы данных

- **[NCBI](https://www.ncbi.nlm.nih.gov/)** — основная база последовательностей
- **[Ensembl](https://www.ensembl.org/)** — аннотированные геномы эукариот
- **[KEGG](https://www.kegg.jp/)** — метаболические пути
- **[UniProt](https://www.uniprot.org/)** — белковые последовательности
- **[Pfam](https://pfam.xfam.org/)** — белковые домены
- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** — клинические варианты
- **[gnomAD](https://gnomad.broadinstitute.org/)** — популяционные частоты

### Туториалы и курсы

- **[Coursera: Bioinformatics](https://www.coursera.org/specializations/bioinformatics)** — специализация UC San Diego
- **[EMBL-EBI Training](https://www.ebi.ac.uk/training/)** — бесплатные онлайн-курсы
- **[Galaxy project](https://galaxyproject.org/)** — интерактивные туториалы

### Книги

- **"Bioinformatics Algorithms"** by Compeau & Pevzner
- **"Biological Sequence Analysis"** by Durbin et al.
- **"Statistical Methods in Bioinformatics"** by Ewens & Grant

---

## 📞 Контакты

- **Вопросы и предложения**: [Issues](https://github.com/username/ngs-bioinformatics-pipeline/issues)
- **Email**: your.email@example.com

---

## 📄 Лицензия

Этот проект распространяется под лицензией MIT — см. файл [LICENSE](LICENSE) для деталей.

---

## 🙏 Благодарности

- Разработчикам инструментов биоинформатики с открытым исходным кодом
- Сообществу GitHub за отличную платформу
- Курсу Blastim за лекционные материалы

---

<div align="center">

**⭐ Если этот проект был вам полезен, поставьте звезду! ⭐**

Made with ❤️ by [Your Name]

</div>
