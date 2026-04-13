#!/bin/bash
# =============================================================================
# Prokka Pipeline - Структурная аннотация генома
# =============================================================================

set -euo pipefail

# Параметры по умолчанию
OUTDIR="./annotation"
PREFIX="genome"
KINGDOM="Bacteria"
CPUS=4
EVALUE="1e-9"
FORCE=false

# Функция помощи
usage() {
    cat << HELP
Использование: $0 [опции] <входной_файл.fasta>

Опции:
  --outdir <дир>    Папка для результатов (по умолчанию: ./annotation)
  --prefix <преф>   Префикс для выходных файлов (по умолчанию: genome)
  --kingdom <тип>   Тип организма: Archaea, Bacteria, Mitochondria, Viruses (по умолчанию: Bacteria)
  --cpus <N>        Количество ядер (по умолчанию: 4)
  --evalue <порог>  Порог E-value для BLAST (по умолчанию: 1e-9)
  --force           Перезаписать существующую директорию
  --help            Показать эту справку

Пример:
  $0 --input genome.fasta --outdir results --kingdom Bacteria --cpus 8
HELP
    exit 1
}

# Парсинг аргументов
INPUT=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --input|-i) INPUT="$2"; shift 2 ;;
        --outdir|-o) OUTDIR="$2"; shift 2 ;;
        --prefix|-p) PREFIX="$2"; shift 2 ;;
        --kingdom|-k) KINGDOM="$2"; shift 2 ;;
        --cpus|-c) CPUS="$2"; shift 2 ;;
        --evalue|-e) EVALUE="$2"; shift 2 ;;
        --force|-f) FORCE=true; shift ;;
        --help|-h) usage ;;
        *) INPUT="$1"; shift ;;
    esac
done

# Проверка входного файла
if [[ -z "$INPUT" ]] || [[ ! -f "$INPUT" ]]; then
    echo "❌ Ошибка: Входной файл не указан или не существует: $INPUT"
    usage
fi

# Проверка зависимостей
check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "❌ Ошибка: $1 не найден. Установите его и попробуйте снова."
        exit 1
    fi
}

echo "🔍 Проверка зависимостей..."
check_dependency prokka
check_dependency prodigal
check_dependency blastp

# Создание директории
if [[ -d "$OUTDIR" ]]; then
    if [[ "$FORCE" == true ]]; then
        echo "⚠️  Директория существует, удаляем: $OUTDIR"
        rm -rf "$OUTDIR"
    else
        echo "❌ Ошибка: Директория уже существует: $OUTDIR (используйте --force)"
        exit 1
    fi
fi

mkdir -p "$OUTDIR"

# Запуск Prokka
echo "🚀 Запуск Prokka..."
echo "   Входной файл: $INPUT"
echo "   Выходная директория: $OUTDIR"
echo "   Префикс: $PREFIX"
echo "   Kingdom: $KINGDOM"
echo "   CPU: $CPUS"
echo "   E-value: $EVALUE"

prokka \
    --outdir "$OUTDIR" \
    --prefix "$PREFIX" \
    --kingdom "$KINGDOM" \
    --cpus "$CPUS" \
    --evalue "$EVALUE" \
    "$INPUT"

# Анализ результатов
echo ""
echo "📊 Анализ результатов..."
FAA_FILE="$OUTDIR/${PREFIX}.faa"
GFF_FILE="$OUTDIR/${PREFIX}.gff"

if [[ -f "$FAA_FILE" ]]; then
    TOTAL_CDS=$(grep -c ">" "$FAA_FILE")
    HYPOTHETIC=$(grep -c "hypothetical protein" "$FAA_FILE" || echo 0)
    RIBOSOMAL=$(grep -c "ribosomal" "$FAA_FILE" || echo 0)
    
    echo "   ✅ Всего CDS: $TOTAL_CDS"
    echo "   ⚪ Гипотетические белки: $HYPOTHETIC"
    echo "   🔬 Рибосомальные белки: $RIBOSOMAL"
fi

# Вывод списка файлов
echo ""
echo "📁 Выходные файлы:"
ls -lh "$OUTDIR/"

echo ""
echo "✅ Аннотация завершена успешно!"
echo "📄 Основной результат: $GFF_FILE"
echo "🧬 Белковые последовательности: $FAA_FILE"
