#!/bin/bash
# ============================================
# 🧬 Genome Annotation Pipeline
# Тема: Аннотация генома (ORF, Prokka, Prodigal, E-value, KEGG, GO)
# Документация: docs/01_genome_annotation.md
# ============================================

set -e

echo "🧬 Запуск пайплайна аннотации генома..."

# Проверка входных данных
INPUT_GENOME="${1:-genome.fasta}"
OUTPUT_DIR="${2:-annotation_output}"

if [ ! -f "$INPUT_GENOME" ]; then
    echo "❌ Ошибка: Файл генома не найден: $INPUT_GENOME"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "📁 Входной файл: $INPUT_GENOME"
echo "📂 Выходная директория: $OUTPUT_DIR"

# Шаг 1: Контроль качества сборки
echo ""
echo "🔍 Шаг 1: Контроль качества сборки..."
if command -v quast &> /dev/null; then
    quast.py "$INPUT_GENOME" -o "$OUTPUT_DIR/quast_report" || echo "⚠️ QUAST не установлен, пропускаем..."
else
    echo "⚠️ QUAST не установлен, пропускаем контроль качества..."
fi

# Шаг 2: Предсказание ORF с помощью Prodigal
echo ""
echo "🔬 Шаг 2: Предсказание ORF (Prodigal)..."
if command -v prodigal &> /dev/null; then
    prodigal -i "$INPUT_GENOME" -a "$OUTPUT_DIR/proteins.faa" -d "$OUTPUT_DIR/genes.fna" -o "$OUTPUT_DIR/prodigal.gff" -f gff || echo "⚠️ Prodigal не установлен, пропускаем..."
else
    echo "⚠️ Prodigal не установлен, пропускаем предсказание ORF..."
fi

# Шаг 3: Аннотация с помощью Prokka
echo ""
echo "🏷️ Шаг 3: Полная аннотация (Prokka)..."
if command -v prokka &> /dev/null; then
    prokka --outdir "$OUTPUT_DIR/prokka" --prefix genome "$INPUT_GENOME" || echo "⚠️ Prokka не установлен, пропускаем..."
else
    echo "⚠️ Prokka не установлен, пропускаем полную аннотацию..."
fi

# Шаг 4: Поиск функциональных доменов (HMMER)
echo ""
echo "🔎 Шаг 4: Поиск функциональных доменов..."
if [ -f "$OUTPUT_DIR/proteins.faa" ] && command -v hmmscan &> /dev/null; then
    hmmscan --domtblout "$OUTPUT_DIR/pfam_domains.tbl" /path/to/Pfam-A.hmm "$OUTPUT_DIR/proteins.faa" || echo "⚠️ HMMER не настроен, пропускаем..."
else
    echo "⚠️ HMMER не доступен или белки не найдены, пропускаем..."
fi

# Шаг 5: Аннотация путей KEGG и GO терминов
echo ""
echo "🛤️ Шаг 5: Аннотация KEGG и GO..."
if command -v eggnoegger-mapper &> /dev/null; then
    emapper.py -i "$OUTPUT_DIR/proteins.faa" --output "$OUTPUT_DIR/eggnog" --go_evidence electronic || echo "⚠️ eggNOG-mapper не установлен, пропускаем..."
else
    echo "⚠️ eggNOG-mapper не установлен, пропускаем KEGG/GO аннотацию..."
fi

# Шаг 6: Генерация отчёта
echo ""
echo "📊 Шаг 6: Генерация итогового отчёта..."
cat > "$OUTPUT_DIR/annotation_summary.txt" << EOF
============================================
🧬 GENOME ANNOTATION SUMMARY
============================================
Input genome: $INPUT_GENOME
Output directory: $OUTPUT_DIR
Date: $(date)

Files generated:
$(ls -lh "$OUTPUT_DIR" | tail -n +2)

Statistics:
- Total genes: $(grep -c "^>" "$OUTPUT_DIR/proteins.faa" 2>/dev/null || echo "N/A")
- Total CDS: $(grep -c "CDS" "$OUTPUT_DIR/prodigal.gff" 2>/dev/null || echo "N/A")
EOF

echo ""
echo "✅ Аннотация завершена!"
echo "📄 Отчёт: $OUTPUT_DIR/annotation_summary.txt"
echo ""
echo "📚 Дополнительная информация: docs/01_genome_annotation.md"
