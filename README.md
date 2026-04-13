# 🧬 NGS Bioinformatics Pipeline

> Учебный репозиторий по анализу данных высокопроизводительного секвенирования: аннотация геномов, вызов вариантов и филогенетический анализ.

## 📋 Оглавление
- [🎯 Цель проекта](#-цель-проекта)
- [🚀 Быстрый старт](#-быстрый-старт)
- [📚 Теоретическая база](#-теоретическая-база)
- [🔧 Пайплайны](#-пайплайны)
- [📐 Формулы и метрики](#-формулы-и-метрики)
- [🧪 Тестирование](#-тестирование)
- [📦 Зависимости](#-зависимости)
- [🤝 Вклад в проект](#-вклад-в-проект)
- [📜 Лицензия](#-лицензия)

## 🎯 Цель проекта
Систематизация знаний и инструментов для анализа NGS-данных с акцентом на:
1. **Структурную и функциональную аннотацию** прокариотических геномов
2. **Вызов и интерпретацию вариантов** (SNV, indel, SV)
3. **Филогенетический анализ** на основе множественных выравниваний

## 🚀 Быстрый старт

### 1. Клонирование
```bash
git clone https://github.com/yourusername/ngs-bioinformatics-pipeline.git
cd ngs-bioinformatics-pipeline
```

---

## 📚 Теоретическая база

Все теоретические материалы с формулами и объяснениями находятся в папке [`docs/`](docs/).

[📥 Скачать все материалы](docs/) | [📐 Формулы](docs/formulas_reference.md) | [🏷️ По тегам](#-теоретическая-база)

| # | Файл | Тема | Статус | Сложность | Теги |
|---|------|------|--------|-----------|------|
| 1️⃣ | [`01_genome_annotation.md`](docs/01_genome_annotation.md) | Аннотация генома | ✅ Готово | 🟢 Базовый | `#ORF` `#Prokka` `#Prodigal` `#KEGG` `#GO` |
| 2️⃣ | [`02_variant_calling.md`](docs/02_variant_calling.md) | Вызов вариантов | ✅ Готово | 🟡 Средний | `#VCF` `#QUAL` `#GQ` `#FreeBayes` `#GATK` |
| 3️⃣ | [`03_structural_variants.md`](docs/03_structural_variants.md) | Структурные варианты | ✅ Готово | 🔴 Продвинутый | `#SVTYPE` `#split-reads` `#discordant-pairs` |
| 4️⃣ | [`04_variant_interpretation.md`](docs/04_variant_interpretation.md) | Интерпретация вариантов | ✅ Готово | 🔴 Продвинутый | `#ACMG` `#gnomAD` `#ClinVar` `#consequence` |
| 5️⃣ | [`05_phylogenetics.md`](docs/05_phylogenetics.md) | Филогенетика | ✅ Готово | 🔴 Продвинутый | `#MSA` `#NJ` `#bootstrap` `#ML` `#Bayesian` |
| 📐 | [`formulas_reference.md`](docs/formulas_reference.md) | Сводка формул | ✅ Готово | 🟢 Базовый | `#формулы` `#метрики` `#статистика` |

### 🏷️ Фильтр по тегам

| Тег | Описание | Файлы |
|-----|----------|-------|
| `#Базовый` | Вводные концепции | [01](docs/01_genome_annotation.md), [formulas](docs/formulas_reference.md) |
| `#Средний` | Практические пайплайны | [02](docs/02_variant_calling.md) |
| `#Продвинутый` | Сложный анализ | [03](docs/03_structural_variants.md), [04](docs/04_variant_interpretation.md), [05](docs/05_phylogenetics.md) |
| `#Формулы` | Математическая база | [01](docs/01_genome_annotation.md), [02](docs/02_variant_calling.md), [03](docs/03_structural_variants.md), [04](docs/04_variant_interpretation.md), [05](docs/05_phylogenetics.md), [formulas](docs/formulas_reference.md) |

> 💡 **Совет**: Используйте `Ctrl+F` для поиска по тегам или перейдите к [полному справочнику формул](docs/formulas_reference.md).
