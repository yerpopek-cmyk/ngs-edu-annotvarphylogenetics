# 2. Вызов вариантов (Variant Calling) — малые варианты (SNV/indel)

## 2.1. Философия и ключевые понятия

**Variant calling ≠ Интерпретация.** Calling отвечает на вопрос: «Есть ли отличие от референса?», интерпретация — «Каковы последствия и насколько это важно?».

Результат вызова — не список мутаций, а набор записей с метриками уверенности:

| Метрика | Описание |
|---------|----------|
| **QUAL** | Phred-качество варианта |
| **FILTER** | Статус фильтрации (PASS или причина) |
| **DP** | Глубина покрытия |
| **AD** | Количество чтений, поддерживающих референсный и альтернативный аллели |

### Основная формула Phred-шкалы

```
QUAL = -10 × log₁₀(P(вариант ошибочный))
```

**QUAL = 30** → вероятность ошибки 0.001 (1/1000)

---

## 2.2. Вероятностная эпоха и формат VCF

### PL (Genotype Likelihoods)

```
PL = -10 × log₁₀(L)
```

где `L` — вероятность наблюдать данные при данном генотипе. Нормализованы так, что лучший генотип имеет PL = 0.

### GQ (Genotype Quality)

Разница между лучшим и вторым лучшим PL (но не более 99).

### Структура VCF (обязательные поля)

| Поле | Описание |
|------|----------|
| **CHROM** | Хромосома/контиг |
| **POS** | Позиция (1-based) |
| **ID** | rsID или «.» |
| **REF** | Референсный аллель |
| **ALT** | Альтернативный аллель (через запятую) |
| **QUAL** | Phred-качество варианта |
| **FILTER** | PASS или причина |
| **INFO** | Дополнительные атрибуты (AN, AC, AF, DP…) |
| **FORMAT** | Перечень полей для образцов |
| **Образцы** | Значения полей FORMAT |

### Ключевые поля FORMAT

| Поле | Описание |
|------|----------|
| **GT** | Генотип (`0/0`, `0/1`, `1/1`, `1/2`; `/` — нефазированный, `|` — фазированный) |
| **DP** | Глубина в образце |
| **AD** | Массив поддержки: `[REF_count, ALT1_count, …]` |
| **GQ** | Качество генотипа (Phred, 0–99) |
| **PL** | Нормализованные логарифмические вероятности |

### Пример строки VCF

```
chr1  12345  .  A  G  50.2  PASS  DP=100  GT:AD:GQ:PL  0/1:45,55:99:0,50,500
```

---

## 2.3. Инструменты вызова вариантов

### 2.3.1. FreeBayes

Haplotype-based caller, строит локальные гаплотипы из чтений.

**Основные параметры:**

```bash
freebayes \
  -f ref.fasta \
  -b aln.bam \
  --min-alternate-count 2 \        # минимум ALT-чтений
  --min-alternate-fraction 0.2 \   # минимальная доля ALT
  --min-base-quality 20 \          # минимальное качество основания
  --pooled-continuous \            # для смешанных популяций
  -v out.vcf
```

### 2.3.2. bcftools

Набор утилит для работы с VCF/BCF.

**Основные команды:**

| Команда | Назначение |
|---------|------------|
| `bcftools view` | Конвертация, фильтрация, выборка образцов |
| `bcftools filter` | Фильтрация по выражениям |
| `bcftools query` | Извлечение полей в текст |
| `bcftools stats` | Статистика по VCF |
| `bcftools norm` | Нормализация |
| `bcftools annotate` | Добавление аннотаций |
| `bcftools merge` | Объединение VCF |
| `bcftools consensus` | Создание консенсусного генома |

**Примеры использования:**

```bash
# Просмотр первых записей
bcftools view calls.vcf.gz | head -20

# Фильтрация по QUAL и DP
bcftools filter -i 'QUAL>30 && INFO/DP>20' calls.vcf.gz -o high_qual.vcf

# Извлечение позиций и генотипов
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\n' calls.vcf.gz

# Статистика
bcftools stats calls.vcf.gz > stats.txt
plot-vcfstats stats.txt -p plots/
```

### 2.3.3. GATK HaplotypeCaller

Стандарт для диплоидных геномов. Строит локальные гаплотипы через граф де Брёйна.

**Минимальный пайплайн (один образец):**

```bash
# Добавление read groups
gatk AddOrReplaceReadGroups -I aln.bam -O rg.bam -ID sample1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample1

# Сортировка и индексация
samtools sort rg.bam -o sorted.bam
samtools index sorted.bam

# Маркировка дубликатов (опционально)
gatk MarkDuplicates -I sorted.bam -O dedup.bam -M metrics.txt

# Вызов вариантов
gatk HaplotypeCaller -R ref.fasta -I dedup.bam -O raw.vcf -ERC GVCF
```

**Жёсткие фильтры (Hard filters) GATK:**

| Фильтр | Порог | Описание |
|--------|-------|----------|
| `QD < 2.0` | QualByDepth | Качество относительно глубины |
| `FS > 60.0` | FisherStrand | Strand bias |
| `MQ < 40.0` | MappingQuality | Качество картирования |
| `MQRankSum < -12.5` | MQRankSum | Разница MQ между REF и ALT |
| `ReadPosRankSum < -8.0` | ReadPosRankSum | Позиция варианта в чтениях |

---

## 2.4. Нормализация вариантов

**Нормализация** — приведение записи варианта к стандартному виду (лево-выравнивание инделей, разделение мультиаллельных сайтов).

**Команды bcftools:**

```bash
# Лево-выравнивание инделей
bcftools norm -f ref.fasta calls.vcf -O z -o norm.vcf.gz

# Разделение мультиаллельных сайтов
bcftools norm -m -any calls.vcf -O v -o split.vcf
```

**Альтернативный инструмент — vt:**

```bash
vt decompose -s -e multiallelic.vcf > decomposed.vcf
vt normalize -r reference.fasta decomposed.vcf > normalized.vcf
```

---

## 2.5. Фильтрация и контроль качества

### Основные метрики для фильтрации

| Метрика | Порог | Описание |
|---------|-------|----------|
| **QUAL > 30** | Phred-качество варианта | Вероятность ошибки < 0.001 |
| **DP > 10** | Минимальная глубина | Достаточное покрытие |
| **MQ > 30** | Качество картирования | Хорошее позиционирование чтений |
| **Allele balance (AB)** | 0.4–0.6 для гетерозигот | `AB = ALT_count / DP` |

### Strand bias (SB)

```
SB = min(ALT_forward, ALT_reverse) / max(ALT_forward, ALT_reverse)
```

Чаще используют тест Фишера. GATK выводит **FS (Phred-scaled p-value):**

```
FS = -10 × log₁₀(p)
```

**FS > 60** → p < 1e⁻⁶ (сильный strand bias)

### Пять красных флажков 🚩

1. 🔴 Низкий MQ / повторы
2. 🔴 Сильный strand bias (FS > 60)
3. 🔴 ALT держится на 1–2 чтениях при низкой глубине
4. 🔴 Вариант на краю чтений / в гомополимере
5. 🔴 Несоответствие AD/DP (странный allele balance)

### Быстрый QC

```bash
bcftools stats calls.vcf.gz | grep "number of SNPs:"
bcftools stats -r TiTv calls.vcf.gz   # Ti/Tv ratio (для человека ~2.0–2.2)
```

---

## 2.6. Функциональная аннотация вариантов

### bcftools csq (быстрый, учитывает фазу)

```bash
bcftools csq -f reference.fasta -g genes.gff3.gz variants.vcf -Ov -o annotated.vcf
```

### snpEff

```bash
java -Xmx4g -jar snpEff.jar GRCh38.99 variants.vcf > annotated.vcf
```

---

## 2.7. Полный пример пайплайна (bash-скрипт)

```bash
#!/bin/bash
REF="hg38.fa"
READS1="sample_R1.fastq.gz"
READS2="sample_R2.fastq.gz"
SAMPLE="S1"

# 1. Выравнивание
bwa mem -t 8 $REF $READS1 $READS2 | samtools sort -@4 -o $SAMPLE.bam
samtools index $SAMPLE.bam

# 2. Вызов вариантов FreeBayes
freebayes -f $REF -b $SAMPLE.bam --min-alternate-fraction 0.2 --min-base-quality 20 > $SAMPLE.raw.vcf

# 3. Сжатие и индексация
bgzip $SAMPLE.raw.vcf
tabix -p vcf $SAMPLE.raw.vcf.gz

# 4. Фильтрация
bcftools filter -i 'QUAL>30 && INFO/DP>10 && INFO/MQ>30' $SAMPLE.raw.vcf.gz -O z -o $SAMPLE.filt.vcf.gz
tabix -p vcf $SAMPLE.filt.vcf.gz

# 5. Статистика
bcftools stats $SAMPLE.filt.vcf.gz > $SAMPLE.stats.txt
plot-vcfstats $SAMPLE.stats.txt -p $SAMPLE.stats_plots/

# 6. Просмотр первых вариантов
bcftools view -H $SAMPLE.filt.vcf.gz | head -20
```

---

## 🔗 Связанные документы

- [📐 Сводка формул](formulas_reference.md) — формулы QUAL, PL, FS
- [🧬 Аннотация генома](01_genome_annotation.md) — предыдущий этап
- [🔍 Структурные варианты](03_structural_variants.md) — следующий уровень
- [📊 Интерпретация](04_variant_interpretation.md) — клиническая значимость
- [🔧 Пайплайны](../pipelines/variant_calling/) — готовые скрипты
