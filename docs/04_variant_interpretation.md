# 📋 Интерпретация вариантов (Variant Interpretation)

## Обзор

Интерпретация вариантов — процесс определения клинической или функциональной значимости обнаруженных генетических вариантов.

## Ключевые понятия

### ACMG Guidelines
Рекомендации American College of Medical Genetics and Genomics для классификации вариантов.

**Категории патогенности:**

| Категория | Код | Описание |
|-----------|-----|----------|
| Pathogenic | P | Патогенный (>99% вероятность) |
| Likely Pathogenic | LP | Вероятно патогенный (90-99%) |
| Uncertain Significance | VUS | Неопределённая значимость |
| Likely Benign | LB | Вероятно безвредный (1-10%) |
| Benign | B | Безвредный (<1%) |

**Критерии ACMG:**
- **PVS1** — Null variant в гене, где LoF является механизмом заболевания
- **PS1-4** — Сильные доказательства патогенности
- **PM1-6** — Умеренные доказательства
- **PP1-5** — Поддерживающие доказательства
- **BA1-BS4** — Доказательства безвредности

### gnomAD (Genome Aggregation Database)
База данных частот аллелей в популяциях.

**Фильтрация по частоте:**
```
Редкий вариант: AF < 0.01 (1%)
Очень редкий: AF < 0.001 (0.1%)
```

### ClinVar
Архив клинически значимых вариантов с аннотациями экспертов.

### Consequence (Тип последствия)
Влияние варианта на ген/транскрипт.

**Иерархия последствий (по убыванию серьёзности):**

| Термин | Описание |
|--------|----------|
| `frameshift_variant` | Сдвиг рамки считывания |
| `stop_gained` | Преждевременный стоп-кодон |
| `splice_donor/acceptor_variant` | Нарушение сплайсинга |
| `missense_variant` | Замена аминокислоты |
| `synonymous_variant` | Синонимичная замена |
| `intron_variant` | В интроне |
| `intergenic_variant` | Между генами |

## Формулы

### Оценка патогенности (Bayesian)
```
Posterior_odds = Prior_odds × LR(ACMG_evidence)
```

где `LR` — likelihood ratio свидетельств ACMG.

### Частота аллелей
```
AF = N(alt_alleles) / N(total_alleles)
```

## Команды

```bash
# Аннотация с SnpEff
java -jar snpEff.jar \
  -v GRCh38.99 \
  input.vcf > annotated.vcf

# Аннотация с VEP (Variant Effect Predictor)
vep \
  --input input.vcf \
  --output output_annotated.vcf \
  --cache \
  --assembly GRCh38 \
  --plugin CADD
```

## Фильтрация по базам данных

```bash
# Отфильтровать варианты с высокой частотой в gnomAD
bcftools view -i 'INFO/gnomAD_AF<0.01' input.vcf > rare_variants.vcf

# Извлечь клинически значимые варианты из ClinVar
bcftools view -i 'INFO/CLNSIG="Pathogenic"' clinvar.vcf
```

## Ссылки

- [ACMG Guidelines](https://www.acmg.net/)
- [gnomAD](https://gnomad.broadinstitute.org/)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [SnpEff](http://snpeff.sourceforge.net/)
- [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)
