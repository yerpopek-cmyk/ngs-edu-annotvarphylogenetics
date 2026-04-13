# 🧬 Аннотация генома

## Обзор

Аннотация генома — процесс идентификации функциональных элементов в последовательности ДНК.

## Ключевые понятия

### ORF (Open Reading Frame)
Последовательность ДНК между старт- и стоп-кодоном, потенциально кодирующая белок.

### Prokka
Инструмент для быстрой аннотации прокариотических геномов.

### Prodigal
Алгоритм предсказания генов в прокариотических геномах.

### E-value
Статистическая значимость выравнивания. Меньшее значение указывает на более значимое совпадение.

**Формула:**
```
E = m × n × 2^(-S)
```
где:
- `m` — длина запроса
- `n` — размер базы данных
- `S` — score выравнивания

### KEGG (Kyoto Encyclopedia of Genes and Genomes)
База данных путей и функциональной информации о генах.

### GO (Gene Ontology)
Система стандартизированной терминологии для описания функций генов.

## Команды

```bash
# Аннотация с Prokka
prokka --outdir annotation --prefix genome genome.fasta

# Предсказание генов с Prodigal
prodigal -i genome.fasta -a proteins.faa -d genes.fna -o output.gff
```

## Ссылки

- [Prokka GitHub](https://github.com/tseemann/prokka)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [KEGG](https://www.kegg.jp/)
- [Gene Ontology](http://geneontology.org/)
