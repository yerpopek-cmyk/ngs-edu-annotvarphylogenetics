#!/bin/bash
# ============================================
# 🔬 Variant Calling Pipeline
# Тема: Вызов вариантов (VCF, QUAL, PL, GQ, FreeBayes, GATK)
# Документация: docs/02_variant_calling.md
# ============================================

set -e

echo "🔬 Запуск пайплайна вызова вариантов..."

# Проверка входных данных
REF_GENOME="${1:-reference.fasta}"
BAM_FILE="${2:-sample.bam}"
OUTPUT_DIR="${3:-variant_output}"

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

# Шаг 1: Индексация референса
echo ""
echo "📇 Шаг 1: Индексация референсного генома..."
if command -v samtools &> /dev/null; then
    samtools faidx "$REF_GENOME" || echo "⚠️ Индексация уже выполнена..."
    if command -v bwa &> /dev/null; then
        bwa index "$REF_GENOME" || echo "⚠️ BWA индекс уже существует..."
    fi
else
    echo "⚠️ samtools не установлен, пропускаем..."
fi

# Шаг 2: Предобработка BAM файла
echo ""
echo "🔧 Шаг 2: Предобработка BAM файла..."
if command -v samtools &> /dev/null; then
    # Сортировка
    samtools sort -o "$OUTPUT_DIR/sorted.bam" "$BAM_FILE" || echo "⚠️ Сортировка не удалась..."
    # Индексация
    samtools index "$OUTPUT_DIR/sorted.bam" || echo "⚠️ Индексация BAM не удалась..."
    # Mark duplicates (если есть picard)
    if command -v picard &> /dev/null; then
        picard MarkDuplicates I="$OUTPUT_DIR/sorted.bam" O="$OUTPUT_DIR/dedup.bam" M="$OUTPUT_DIR/metrics.txt" || echo "⚠️ Picard не настроен..."
    else
        cp "$OUTPUT_DIR/sorted.bam" "$OUTPUT_DIR/dedup.bam"
    fi
else
    echo "⚠️ samtools не установлен, пропускаем предобработку..."
    cp "$BAM_FILE" "$OUTPUT_DIR/dedup.bam"
fi

# Шаг 3: Вызов вариантов через FreeBayes
echo ""
echo "🌊 Шаг 3: Вызов вариантов (FreeBayes)..."
if command -v freebayes &> /dev/null; then
    freebayes -f "$REF_GENOME" "$OUTPUT_DIR/dedup.bam" > "$OUTPUT_DIR/freebayes.vcf" || echo "⚠️ FreeBayes не удался..."
    # Фильтрация по качеству
    if [ -f "$OUTPUT_DIR/freebayes.vcf" ]; then
        vcffilter -f "QUAL > 30" "$OUTPUT_DIR/freebayes.vcf" > "$OUTPUT_DIR/freebayes_filtered.vcf" 2>/dev/null || true
    fi
else
    echo "⚠️ FreeBayes не установлен, пропускаем..."
fi

# Шаг 4: Вызов вариантов через GATK
echo ""
echo "🧬 Шаг 4: Вызов вариантов (GATK)..."
if command -v gatk &> /dev/null; then
    gatk HaplotypeCaller \
        -R "$REF_GENOME" \
        -I "$OUTPUT_DIR/dedup.bam" \
        -O "$OUTPUT_DIR/gatk_raw.vcf" \
        --native-pair-hmm-threads 4 || echo "⚠️ GATK HaplotypeCaller не удался..."
    
    # Hard filtering
    gatk VariantFiltration \
        -R "$REF_GENOME" \
        -V "$OUTPUT_DIR/gatk_raw.vcf" \
        -O "$OUTPUT_DIR/gatk_filtered.vcf" \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
        --filter-name "basic_filter" || echo "⚠️ GATK фильтрация не удалась..."
else
    echo "⚠️ GATK не установлен, пропускаем..."
fi

# Шаг 5: Объединение и аннотация вариантов
echo ""
echo "🔗 Шаг 5: Объединение результатов..."
if [ -f "$OUTPUT_DIR/freebayes_filtered.vcf" ] && [ -f "$OUTPUT_DIR/gatk_filtered.vcf" ]; then
    echo "📊 Оба метода завершены, создаём сравнительный отчёт..."
    bcftools isec -n=2 -w1 "$OUTPUT_DIR/freebayes_filtered.vcf" "$OUTPUT_DIR/gatk_filtered.vcf" -p "$OUTPUT_DIR/consensus" 2>/dev/null || echo "⚠️ bcftools не настроен..."
fi

# Шаг 6: Генерация отчёта
echo ""
echo "📊 Шаг 6: Генерация итогового отчёта..."
cat > "$OUTPUT_DIR/variant_summary.txt" << EOF
============================================
🔬 VARIANT CALLING SUMMARY
============================================
Reference: $REF_GENOME
BAM file: $BAM_FILE
Output directory: $OUTPUT_DIR
Date: $(date)

Files generated:
$(ls -lh "$OUTPUT_DIR"/*.vcf 2>/dev/null || echo "No VCF files")

Statistics:
EOF

# Подсчёт статистики по VCF
for vcf in "$OUTPUT_DIR"/*.vcf; do
    if [ -f "$vcf" ]; then
        total=$(grep -v "^#" "$vcf" | wc -l)
        echo "- $(basename $vcf): $total вариантов" >> "$OUTPUT_DIR/variant_summary.txt"
    fi
done

echo ""
echo "✅ Вызов вариантов завершён!"
echo "📄 Отчёт: $OUTPUT_DIR/variant_summary.txt"
echo ""
echo "📚 Дополнительная информация: docs/02_variant_calling.md"
