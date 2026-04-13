#!/bin/bash
# ============================================
# 🚀 Bioinformatics Pipelines Runner
# Запуск всех пайплайнов или выборочно
# ============================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCS_DIR="$(dirname "$SCRIPT_DIR")/docs"

echo "🚀 Bioinformatics Pipelines Suite"
echo "================================="
echo ""
echo "Доступные пайплайны:"
echo "  1. 🧬 Genome Annotation (01_genome_annotation.sh)"
echo "  2. 🔬 Variant Calling (02_variant_calling.sh)"
echo "  3. 🧩 Structural Variants (03_structural_variants.sh)"
echo "  4. 🧾 Variant Interpretation (04_variant_interpretation.sh)"
echo "  5. 🌳 Phylogenetics (05_phylogenetics.sh)"
echo "  6. 📐 Formulas Reference (formulas_reference.sh)"
echo ""

# Проверка аргументов
if [ $# -eq 0 ]; then
    echo "Использование:"
    echo "  ./run_all.sh                    # Запустить все пайплайны"
    echo "  ./run_all.sh 1                  # Только аннотация генома"
    echo "  ./run_all.sh 2                  # Только variant calling"
    echo "  ./run_all.sh 1 3 5              # Несколько пайплайнов"
    echo "  ./run_all.sh all                # Все пайплайны"
    echo ""
    echo "Примеры запуска отдельных скриптов:"
    echo "  ./01_genome_annotation.sh genome.fasta output_dir"
    echo "  ./02_variant_calling.sh reference.fasta sample.bam output_dir"
    echo "  ./03_structural_variants.sh reference.fasta sample.bam output_dir"
    echo "  ./04_variant_interpretation.sh input.vcf output_dir GRCh38"
    echo "  ./05_phylogenetics.sh sequences.fasta output_dir all"
    echo "  ./formulas_reference.sh output_dir"
    exit 0
fi

# Функция запуска пайплайна
run_pipeline() {
    local pipeline="$1"
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "🔄 Запуск: $pipeline"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    case "$pipeline" in
        1|genome|annotation)
            echo "⚠️  Требуется входной файл генома"
            echo "   Пример: ./01_genome_annotation.sh genome.fasta annotation_output"
            ;;
        2|variant|calling)
            echo "⚠️  Требуется референс и BAM файл"
            echo "   Пример: ./02_variant_calling.sh reference.fasta sample.bam variant_output"
            ;;
        3|sv|structural)
            echo "⚠️  Требуется референс и BAM файл"
            echo "   Пример: ./03_structural_variants.sh reference.fasta sample.bam sv_output"
            ;;
        4|interpret|interpretation)
            echo "⚠️  Требуется VCF файл"
            echo "   Пример: ./04_variant_interpretation.sh input.vcf interpretation_output GRCh38"
            ;;
        5|phylo|phylogenetics)
            echo "⚠️  Требуется файл последовательностей"
            echo "   Пример: ./05_phylogenetics.sh sequences.fasta phylo_output all"
            ;;
        6|formulas|reference)
            echo "📐 Генерация справочника формул..."
            "$SCRIPT_DIR/formulas_reference.sh" "$SCRIPT_DIR/../formulas_output"
            ;;
        all)
            echo "🔄 Запуск всех доступных демо-пайплайнов..."
            "$SCRIPT_DIR/formulas_reference.sh" "$SCRIPT_DIR/../formulas_output"
            ;;
        *)
            echo "❌ Неизвестный пайплайн: $pipeline"
            return 1
            ;;
    esac
}

# Запуск выбранных пайплайнов
for arg in "$@"; do
    run_pipeline "$arg"
done

echo ""
echo "================================="
echo "✅ Все запрошенные пайплайны завершены!"
echo ""
echo "📚 Документация:"
echo "   - docs/01_genome_annotation.md"
echo "   - docs/02_variant_calling.md"
echo "   - docs/03_structural_variants.md"
echo "   - docs/04_variant_interpretation.md"
echo "   - docs/05_phylogenetics.md"
echo "   - docs/formulas_reference.md"
echo ""
echo "🔗 GitHub README: ../../README.md"
echo "================================="
