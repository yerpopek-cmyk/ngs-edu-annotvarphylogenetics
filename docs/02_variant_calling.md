# 🔬 Вызов вариантов (Variant Calling)

## Обзор

Вызов вариантов — процесс идентификации различий между исследуемым геномом и референсом.

## Ключевые понятия

### VCF (Variant Call Format)
Стандартный формат для хранения информации о вариантах.

**Структура VCF:**
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
chr1    100     .       A       G       30      PASS    .       GT:PL   0/1:0,30,300
```

### QUAL (Quality Score)
Логарифмическая оценка достоверности варианта.

**Формула:**
```
QUAL = -10 × log₁₀(P(error))
```
где `P(error)` — вероятность ошибки вызова варианта.

**Пример:** QUAL=30 означает вероятность ошибки 1/1000 (99.9% точность).

### PL (Phred-scaled Likelihoods)
Нормализованные вероятности генотипов.

**Формула:**
```
PL(g) = -10 × log₁₀(L(g))
```
где `L(g)` — правдоподобие генотипа `g`.

### GQ (Genotype Quality)
Уверенность в присвоенном генотипе.

**Формула:**
```
GQ = min(PL(other_genotypes)) - PL(best_genotype)
```

### FreeBayes
Байесовский инструмент для вызова вариантов.

### GATK (Genome Analysis Toolkit)
Набор инструментов от Broad Institute для анализа вариантов.

## Команды

```bash
# Вызов вариантов с FreeBayes
freebayes -f reference.fasta reads.bam > variants.vcf

# Вызов вариантов с GATK
gatk HaplotypeCaller \
  -R reference.fasta \
  -I reads.bam \
  -O variants.vcf
```

## Фильтрация вариантов

```bash
# Фильтрация по качеству
bcftools filter -i 'QUAL>30 && GQ>20' input.vcf > filtered.vcf
```

## Ссылки

- [VCF Specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [FreeBayes](https://github.com/freebayes/freebayes)
- [GATK](https://gatk.broadinstitute.org/)
