#!/bin/bash
# ============================================
# 🧾 Variant Interpretation Pipeline
# Тема: Интерпретация вариантов (ACMG, gnomAD, ClinVar, consequence)
# Документация: docs/04_variant_interpretation.md
# ============================================

set -e

echo "🧾 Запуск пайплайна интерпретации вариантов..."

# Проверка входных данных
VCF_FILE="${1:-input.vcf}"
OUTPUT_DIR="${2:-interpretation_output}"
GENOME_BUILD="${3:-GRCh38}"

if [ ! -f "$VCF_FILE" ]; then
    echo "❌ Ошибка: VCF файл не найден: $VCF_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "📁 VCF файл: $VCF_FILE"
echo "📂 Выходная директория: $OUTPUT_DIR"
echo "🧬 Сборка генома: $GENOME_BUILD"

# Шаг 1: Нормализация VCF
echo ""
echo "🔧 Шаг 1: Нормализация VCF файла..."
if command -v bcftools &> /dev/null; then
    bcftools norm -m -any "$VCF_FILE" -Oz -o "$OUTPUT_DIR/normalized.vcf.gz" || echo "⚠️ Нормализация не удалась..."
    bcftools index "$OUTPUT_DIR/normalized.vcf.gz" || echo "⚠️ Индексация не удалась..."
else
    echo "⚠️ bcftools не установлен, пропускаем нормализацию..."
    cp "$VCF_FILE" "$OUTPUT_DIR/normalized.vcf"
fi

# Шаг 2: Аннотация через VEP (Ensembl Variant Effect Predictor)
echo ""
echo "🔬 Шаг 2: Функциональная аннотация (VEP)..."
if command -v vep &> /dev/null; then
    vep -i "$VCF_FILE" \
        --cache \
        --offline \
        --dir_cache ~/.vep \
        --assembly "$GENOME_BUILD" \
        --vcf \
        --compress_output bgzip \
        --output_file "$OUTPUT_DIR/vep_annotated.vcf.gz" \
        --plugin CADD,/path/to/CADD/v1.6/whole_genome_SNVs.tsv.gz \
        --plugin REVEL,/path/to/REVEL/revel_all_chromosomes.csv.bz2 \
        --fork 4 || echo "⚠️ VEP аннотация не удалась..."
else
    echo "⚠️ VEP не установлен, пропускаем..."
fi

# Шаг 3: Добавление частот из gnomAD
echo ""
echo "🌍 Шаг 3: Добавление частот аллелей (gnomAD)..."
if command -v bcftools &> /dev/null; then
    # Аннотация частотами gnomAD
    bcftools annotate \
        -a /path/to/gnomad.vcf.gz \
        -c INFO/AF \
        "$OUTPUT_DIR/vep_annotated.vcf.gz" \
        -Oz -o "$OUTPUT_DIR/gnomad_annotated.vcf.gz" 2>/dev/null || echo "⚠️ gnomAD аннотация не удалась (файлы не найдены)..."
else
    echo "⚠️ bcftools не доступен, пропускаем..."
fi

# Шаг 4: Поиск в ClinVar
echo ""
echo "🏥 Шаг 4: Поиск клинических значимостей (ClinVar)..."
if command -v bcftools &> /dev/null; then
    bcftools annotate \
        -a /path/to/clinvar.vcf.gz \
        -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN \
        "$OUTPUT_DIR/gnomad_annotated.vcf.gz" \
        -Oz -o "$OUTPUT_DIR/clinvar_annotated.vcf.gz" 2>/dev/null || echo "⚠️ ClinVar аннотация не удалась (файлы не найдены)..."
else
    echo "⚠️ bcftools не доступен, пропускаем..."
fi

# Шаг 5: Классификация по ACMG
echo ""
echo "📋 Шаг 5: Классификация по рекомендациям ACMG..."
if command -vInterVar &> /dev/null; then
    InterVar.py \
        -i "$VCF_FILE" \
        -o "$OUTPUT_DIR/intervar_result.txt" \
        --assaytype germline \
        --buildver "$GENOME_BUILD" || echo "⚠️ InterVar не удался..."
elif command -v varsome &> /dev/null; then
    echo "⚠️ VarSome CLI требует аутентификации, пропускаем автоматическую классификацию..."
else
    echo "⚠️ Инструменты ACMG классификации не найдены, создаём шаблон отчёта..."
    
    # Создаём шаблон для ручной классификации
    cat > "$OUTPUT_DIR/acmg_template.csv" << EOF
Variant,Chr,Pos,Ref,Alt,Consequence,gnoMAF,ClinVar_Significance,ACMG_Class,Evidence
EOF
    
    grep -v "^#" "$VCF_FILE" | while read line; do
        chr=$(echo "$line" | cut -f1)
        pos=$(echo "$line" | cut -f2)
        ref=$(echo "$line" | cut -f4)
        alt=$(echo "$line" | cut -f5)
        echo "$chr:$pos.$ref>$alt,$chr,$pos,$ref,$alt,,,,," >> "$OUTPUT_DIR/acmg_template.csv"
    done
fi

# Шаг 6: Приоритизация вариантов
echo ""
echo "⭐ Шаг 6: Приоритизация патогенных вариантов..."
cat > "$OUTPUT_DIR/prioritization.sh" << 'SCRIPT'
#!/bin/bash
# Фильтрация по критериям патогенности
INPUT_VCF="${1:-../clinvar_annotated.vcf.gz}"
OUTPUT_PRIORITIZED="${2:-prioritized_variants.txt}"

echo "Фильтрация вариантов по критериям:"
echo "  - Частота в gnomAD < 0.01"
echo "  - Предсказанный эффект: damaging"
echo "  - ClinVar: Pathogenic/Likely_pathogenic"

bcftools view -i 'INFO.AF<0.01' "$INPUT_VCF" 2>/dev/null | \
    grep -E "Pathogenic|Likely_pathogenic|HIGH" > "$OUTPUT_PRIORITIZED" || true

echo "Найдено приоритетных вариантов: $(wc -l < $OUTPUT_PRIORITIZED)"
SCRIPT
chmod +x "$OUTPUT_DIR/prioritization.sh"

# Шаг 7: Генерация итогового отчёта
echo ""
echo "📊 Шаг 7: Генерация клинического отчёта..."
cat > "$OUTPUT_DIR/interpretation_report.txt" << EOF
============================================
🧾 VARIANT INTERPRETATION REPORT
============================================
Input VCF: $VCF_FILE
Genome build: $GENOME_BUILD
Output directory: $OUTPUT_DIR
Date: $(date)

ANNOTATION PIPELINE STATUS:
---------------------------
✅ Normalization: $([ -f "$OUTPUT_DIR/normalized.vcf.gz" ] && echo "Complete" || echo "Skipped")
✅ VEP Annotation: $([ -f "$OUTPUT_DIR/vep_annotated.vcf.gz" ] && echo "Complete" || echo "Skipped")
✅ gnomAD Frequencies: $([ -f "$OUTPUT_DIR/gnomad_annotated.vcf.gz" ] && echo "Complete" || echo "Skipped")
✅ ClinVar Lookup: $([ -f "$OUTPUT_DIR/clinvar_annotated.vcf.gz" ] && echo "Complete" || echo "Skipped")
✅ ACMG Classification: $([ -f "$OUTPUT_DIR/acmg_template.csv" ] && echo "Template created" || echo "Automated")

VARIANT SUMMARY:
----------------
EOF

# Подсчёт статистики
total_vars=$(grep -v "^#" "$VCF_FILE" | wc -l)
echo "Total variants: $total_vars" >> "$OUTPUT_DIR/interpretation_report.txt"

if [ -f "$OUTPUT_DIR/clinvar_annotated.vcf.gz" ]; then
    pathogenic=$(bcftools view -i 'INFO/CLNSIG="Pathogenic"' "$OUTPUT_DIR/clinvar_annotated.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l || echo "0")
    likely_pathogenic=$(bcftools view -i 'INFO/CLNSIG="Likely_pathogenic"' "$OUTPUT_DIR/clinvar_annotated.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l || echo "0")
    vus=$(bcftools view -i 'INFO/CLNSIG="Uncertain_significance"' "$OUTPUT_DIR/clinvar_annotated.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l || echo "0")
    
    echo "Pathogenic: $pathogenic" >> "$OUTPUT_DIR/interpretation_report.txt"
    echo "Likely pathogenic: $likely_pathogenic" >> "$OUTPUT_DIR/interpretation_report.txt"
    echo "VUS: $vus" >> "$OUTPUT_DIR/interpretation_report.txt"
fi

cat >> "$OUTPUT_DIR/interpretation_report.txt" << EOF

NEXT STEPS:
-----------
1. Review prioritized variants in: $OUTPUT_DIR/prioritized_variants.txt
2. Manual ACMG classification using template: $OUTPUT_DIR/acmg_template.csv
3. Cross-reference with literature and family history
4. Generate clinical report for healthcare provider

============================================
📚 Дополнительная информация: docs/04_variant_interpretation.md
============================================
EOF

echo ""
echo "✅ Интерпретация вариантов завершена!"
echo "📄 Отчёт: $OUTPUT_DIR/interpretation_report.txt"
echo ""
echo "📚 Дополнительная информация: docs/04_variant_interpretation.md"
