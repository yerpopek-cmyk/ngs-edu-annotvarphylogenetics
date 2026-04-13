#!/bin/bash
# =============================================================================
# FreeBayes Pipeline - Вызов вариантов (SNV/indel)
# =============================================================================

set -euo pipefail

# Параметры по умолчанию
REFERENCE=""
BAM=""
OUTPUT="variants.vcf"
MIN_ALT_COUNT=2
MIN_ALT_FRAC=0.2
MIN_BASE_QUAL=20
CPUS=4

usage() {
    cat << HELP
Использование: $0 [опции]

Обязательные опции:
  --reference|-r <файл>  Референсный геном (fasta)
  --bam|-b <файл>       Выровненные чтения (BAM)

Опциональные опции:
  --output|-o <файл>    Выходной VCF (по умолчанию: variants.vcf)
  --min-alt-count <N>   Минимум ALT-чтений (по умолчанию: 2)
  --min-alt-frac <F>    Минимальная доля ALT (по умолчанию: 0.2)
  --min-base-qual <Q>   Минимальное качество основания (по умолчанию: 20)
  --cpus <N>            Количество ядер (по умолчанию: 4)
  --help|-h             Показать эту справку

Пример:
  $0 -r ref.fasta -b sample.bam -o sample.vcf --min-alt-frac 0.2
HELP
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --reference|-r) REFERENCE="$2"; shift 2 ;;
        --bam|-b) BAM="$2"; shift 2 ;;
        --output|-o) OUTPUT="$2"; shift 2 ;;
        --min-alt-count) MIN_ALT_COUNT="$2"; shift 2 ;;
        --min-alt-frac) MIN_ALT_FRAC="$2"; shift 2 ;;
        --min-base-qual) MIN_BASE_QUAL="$2"; shift 2 ;;
        --cpus|-c) CPUS="$2"; shift 2 ;;
        --help|-h) usage ;;
        *) echo "Неизвестная опция: $1"; usage ;;
    esac
done

if [[ -z "$REFERENCE" ]] || [[ ! -f "$REFERENCE" ]]; then
    echo "❌ Ошибка: Референсный файл не указан или не существует"
    usage
fi

if [[ -z "$BAM" ]] || [[ ! -f "$BAM" ]]; then
    echo "❌ Ошибка: BAM файл не указан или не существует"
    usage
fi

check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "❌ Ошибка: $1 не найден"
        exit 1
    fi
}

echo "🔍 Проверка зависимостей..."
check_dependency freebayes
check_dependency samtools
check_dependency bcftools

echo "🚀 Запуск FreeBayes..."
freebayes \
    -f "$REFERENCE" \
    -b "$BAM" \
    --min-alternate-count "$MIN_ALT_COUNT" \
    --min-alternate-fraction "$MIN_ALT_FRAC" \
    --min-base-quality "$MIN_BASE_QUAL" \
    --threads "$CPUS" \
    -v "$OUTPUT"

echo "📦 Сжатие и индексация VCF..."
bgzip -f "$OUTPUT"
tabix -p vcf "${OUTPUT}.gz"

echo "📊 Статистика вариантов..."
bcftools stats "${OUTPUT}.gz" > "${OUTPUT}.stats.txt"

echo "✅ Вызов вариантов завершён!"
echo "📄 Результат: ${OUTPUT}.gz"
