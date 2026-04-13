#!/bin/bash
# =============================================================================
# IQ-TREE Pipeline - Филогенетический анализ
# =============================================================================

set -euo pipefail

ALIGNMENT=""
OUTDIR="./phylogeny"
BOOTSTRAP=1000
SEQTYPE="DNA"
THREADS=4
MODEL="MFP"

usage() {
    cat << HELP
Использование: $0 [опции]

Обязательные опции:
  --alignment|-a <файл>  Выравнивание (fasta)

Опциональные опции:
  --outdir|-o <дир>      Выходная директория (по умолчанию: ./phylogeny)
  --bootstrap|-b <N>     Количество bootstrap репликаций (по умолчанию: 1000)
  --seqtype|-s <тип>     Тип данных: DNA, AA, CODON (по умолчанию: DNA)
  --threads|-t <N>       Количество потоков (по умолчанию: 4)
  --model|-m <модель>    Модель эволюции или MFP для автовыбора (по умолчанию: MFP)
  --help|-h              Показать эту справку

Пример:
  $0 -a alignment.fasta -b 1000 -t 8
HELP
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --alignment|-a) ALIGNMENT="$2"; shift 2 ;;
        --outdir|-o) OUTDIR="$2"; shift 2 ;;
        --bootstrap|-b) BOOTSTRAP="$2"; shift 2 ;;
        --seqtype|-s) SEQTYPE="$2"; shift 2 ;;
        --threads|-t) THREADS="$2"; shift 2 ;;
        --model|-m) MODEL="$2"; shift 2 ;;
        --help|-h) usage ;;
        *) echo "Неизвестная опция: $1"; usage ;;
    esac
done

if [[ -z "$ALIGNMENT" ]] || [[ ! -f "$ALIGNMENT" ]]; then
    echo "❌ Ошибка: Файл выравнивания не указан или не существует"
    usage
fi

check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "❌ Ошибка: $1 не найден"
        exit 1
    fi
}

echo "🔍 Проверка зависимостей..."
check_dependency iqtree3

mkdir -p "$OUTDIR"

echo "🚀 Запуск IQ-TREE..."
echo "   Выравнивание: $ALIGNMENT"
echo "   Bootstrap: $BOOTSTRAP"
echo "   Тип: $SEQTYPE"
echo "   Потоки: $THREADS"
echo "   Модель: $MODEL"

iqtree3 \
    -s "$ALIGNMENT" \
    -m "$MODEL" \
    -B "$BOOTSTRAP" \
    --seqtype "$SEQTYPE" \
    -T "$THREADS" \
    -pre "$OUTDIR/tree"

echo ""
echo "✅ Филогенетический анализ завершён!"
echo "📄 Дерево: $OUTDIR/tree.treefile"
echo "📊 Статистика: $OUTDIR/tree.iqtree"
