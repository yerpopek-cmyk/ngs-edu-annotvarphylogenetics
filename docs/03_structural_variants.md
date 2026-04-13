# 🧬 Структурные варианты (Structural Variants)

## Обзор

Структурные варианты (SV) — крупные изменения в геноме (>50 п.н.): делеции, дупликации, инверсии, транслокации.

## Ключевые понятия

### SVTYPE
Тип структурного варианта в VCF формате:

| Код | Тип | Описание |
|-----|-----|----------|
| `DEL` | Делеция | Потеря участка ДНК |
| `DUP` | Дупликация | Удвоение участка ДНК |
| `INV` | Инверсия | Разворот участка ДНК |
| `INS` | Инсерция | Вставка нового материала |
| `BND` | Breakend | Транслокация/перестройка |
| `CNV` | Copy Number Variant | Изменение числа копий |

### Split-reads
Метод обнаружения SV по ридам, которые частично выравниваются в одном месте, а частично — в другом.

**Принцип:**
```
Референс:  AAAAAAAAAA----------TTTTTTTTTT
Рид:       AAAAAAAAAA          TTTTTTTTTT
           ||||||||||          ||||||||||
Разрыв указывает на точку разрыва SV
```

### Discordant pairs
Пары ридов с неожиданным расстоянием или ориентацией.

**Нормальная пара:**
```
R1 --->        <--- R2
     |--------|
     insert_size ~ 300-500 bp
```

**Дискордантная пара (делеция):**
```
R1 --->                    <--- R2
     |--------------------|
     insert_size >> 500 bp
```

**Дискордантная пара (инверсия):**
```
R1 --->        ---> R2    (неправильная ориентация)
```

## Формулы

### Оценка размера делеции
```
SV_size = observed_insert_size - expected_insert_size
```

### Поддержка split-read
```
Support = N(split_reads) / N(total_reads_at_locus)
```

## Методы обнаружения

| Метод | Преимущества | Недостатки |
|-------|-------------|------------|
| Split-reads | Точные breakpoints | Требует высокого покрытия |
| Discordant pairs | Обнаруживает крупные SV | Низкое разрешение |
| Read depth | Обнаружает CNV | Не определяет тип точно |
| De novo assembly | Полная реконструкция | Вычислительно затратно |

## Команды

```bash
# Manta для обнаружения SV
configManta.py \
  --bam reads.bam \
  --referenceFasta reference.fasta \
  --runDir manta_run
manta_run/runWorkflow.py -m local -j 8

# Delly для SV
delly call \
  -g reference.fasta \
  -o variants.bcf \
  reads.bam
```

## Ссылки

- [Manta](https://github.com/Illumina/manta)
- [Delly](https://github.com/dellytools/delly)
- [LUMPY](https://github.com/arq5x/lumpy-sv)
