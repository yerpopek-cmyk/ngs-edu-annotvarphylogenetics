#!/bin/bash
# ============================================
# 🧩 Structural Variants Detection Pipeline
# Тема: Структурные варианты (SVTYPE, split-reads, discordant pairs)
# Документация: docs/03_structural_variants.md
# ============================================

set -e

echo "🧩 Запуск пайплайна поиска структурных вариантов..."

# Проверка входных данных
REF_GENOME="${1:-reference.fasta}"
BAM_FILE="${2:-sample.bam}"
OUTPUT_DIR="${3:-sv_output}"

if [ ! -f "$REF_GENOME" ]; then
    echo "❌ Ошибка: Референсный геном не найден: $REF_GENOME"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo "❌ Ошибка: BAM файл не найден: $BAM_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "📁 Референс: $REF_GENOME"
echo "📁 BAM файл: $BAM_FILE"
echo "📂 Выходная директория: $OUTPUT_DIR"

# Шаг 1: Подготовка BAM файла
echo ""
echo "🔧 Шаг 1: Подготовка BAM файла..."
if command -v samtools &> /dev/null; then
    samtools sort -n -o "$OUTPUT_DIR/name_sorted.bam" "$BAM_FILE" || echo "⚠️ Сортировка по имени не удалась..."
    samtools index "$OUTPUT_DIR/name_sorted.bam" || echo "⚠️ Индексация не удалась..."
else
    echo "⚠️ samtools не установлен, пропускаем..."
    cp "$BAM_FILE" "$OUTPUT_DIR/name_sorted.bam"
fi

# Шаг 2: Поиск SV через Delly
echo ""
echo "🔍 Шаг 2: Поиск SV (Delly)..."
if command -v delly &> /dev/null; then
    delly call -g "$REF_GENOME" -o "$OUTPUT_DIR/delly.bcf" "$BAM_FILE" || echo "⚠️ Delly не удался..."
    bcftools convert "$OUTPUT_DIR/delly.bcf" -Ov -o "$OUTPUT_DIR/delly.vcf" 2>/dev/null || true
else
    echo "⚠️ Delly не установлен, пропускаем..."
fi

# Шаг 3: Поиск SV через Manta
echo ""
echo "🌀 Шаг 3: Поиск SV (Manta)..."
if command -v configManta.py &> /dev/null; then
    configManta.py --bamFile "$BAM_FILE" --referenceFasta "$REF_GENOME" --runDir "$OUTPUT_DIR/manta_run" || echo "⚠️ Manta настройка не удалась..."
    if [ -d "$OUTPUT_DIR/manta_run" ]; then
        cd "$OUTPUT_DIR/manta_run" && ./runWorkflow.py -m local || echo "⚠️ Manta запуск не удался..."
        cd - > /dev/null
    fi
else
    echo "⚠️ Manta не установлена, пропускаем..."
fi

# Шаг 4: Поиск SV через Lumpy
echo ""
echo "🧪 Шаг 4: Поиск SV (Lumpy)..."
if command -v lumpyexpress &> /dev/null; then
    # Извлечение discordant pairs и split reads
    if command -v samblaster &> /dev/null; then
        samblaster -i "$BAM_FILE" \
            -s "$OUTPUT_DIR/split.sam" \
            -d "$OUTPUT_DIR/discordant.sam" || echo "⚠️ samblaster не удался..."
    fi
    
    lumpyexpress \
        -B "$BAM_FILE" \
        -S "$OUTPUT_DIR/split.sam" \
        -D "$OUTPUT_DIR/discordant.sam" \
        -o "$OUTPUT_DIR/lumpy.vcf" || echo "⚠️ Lumpy не удался..."
else
    echo "⚠️ Lumpy не установлен, пропускаем..."
fi

# Шаг 5: Фильтрация и объединение результатов
echo ""
echo "🔗 Шаг 5: Объединение и фильтрация SV..."
if command -v SURVIVOR &> /dev/null; then
    # Создание списка VCF файлов
    find "$OUTPUT_DIR" -name "*.vcf" -type f > "$OUTPUT_DIR/vcflist.txt"
    SURVIVOR merge "$OUTPUT_DIR/vcflist.txt" 1000 2 1 0 0 50 "$OUTPUT_DIR/consensus.vcf" || echo "⚠️ SURVIVOR merge не удался..."
    
    # Фильтрация по качеству
    SURVIVOR filter "$OUTPUT_DIR/consensus.vcf" 10 1000 0 0 0 0 "$OUTPUT_DIR/filtered.vcf" || echo "⚠️ SURVIVOR filter не удался..."
else
    echo "⚠️ SURVIVOR не установлен, пропускаем объединение..."
fi

# Шаг 6: Аннотация SV
echo ""
echo "🏷️ Шаг 6: Аннотация структурных вариантов..."
if command -v svanno &> /dev/null; then
    svanno -i "$OUTPUT_DIR/filtered.vcf" -r "$REF_GENOME" -o "$OUTPUT_DIR/annotated.vcf" || echo "⚠️ Аннотация SV не удалась..."
elif command -v vep &> /dev/null; then
    vep -i "$OUTPUT_DIR/filtered.vcf" -o "$OUTPUT_DIR/annotated.vcf" --vcf --offline || echo "⚠️ VEP аннотация не удалась..."
else
    echo "⚠️ Инструменты аннотации SV не найдены, пропускаем..."
fi

# Шаг 7: Генерация отчёта
echo ""
echo "📊 Шаг 7: Генерация итогового отчёта..."
cat > "$OUTPUT_DIR/sv_summary.txt" << EOF
============================================
🧩 STRUCTURAL VARIANTS SUMMARY
============================================
Reference: $REF_GENOME
BAM file: $BAM_FILE
Output directory: $OUTPUT_DIR
Date: $(date)

Detected SV types:
EOF

# Подсчёт по типам SV
for vcf in "$OUTPUT_DIR"/*.vcf; do
    if [ -f "$vcf" ]; then
        echo "" >> "$OUTPUT_DIR/sv_summary.txt"
        echo "📄 Файл: $(basename $vcf)" >> "$OUTPUT_DIR/sv_summary.txt"
        
        del_count=$(grep -v "^#" "$vcf" 2>/dev/null | grep -c "DEL" || echo "0")
        dup_count=$(grep -v "^#" "$vcf" 2>/dev/null | grep -c "DUP" || echo "0")
        inv_count=$(grep -v "^#" "$vcf" 2>/dev/null | grep -c "INV" || echo "0")
        tra_count=$(grep -v "^#" "$vcf" 2>/dev/null | grep -c "TRA" || echo "0")
        ins_count=$(grep -v "^#" "$vcf" 2>/dev/null | grep -c "INS" || echo "0")
        
        echo "  - DEL (deletions): $del_count" >> "$OUTPUT_DIR/sv_summary.txt"
        echo "  - DUP (duplications): $dup_count" >> "$OUTPUT_DIR/sv_summary.txt"
        echo "  - INV (inversions): $inv_count" >> "$OUTPUT_DIR/sv_summary.txt"
        echo "  - TRA (translocations): $tra_count" >> "$OUTPUT_DIR/sv_summary.txt"
        echo "  - INS (insertions): $ins_count" >> "$OUTPUT_DIR/sv_summary.txt"
    fi
done

echo ""
echo "✅ Поиск структурных вариантов завершён!"
echo "📄 Отчёт: $OUTPUT_DIR/sv_summary.txt"
echo ""
echo "📚 Дополнительная информация: docs/03_structural_variants.md"
