#!/bin/bash
# ============================================
# 🌳 Phylogenetics Analysis Pipeline
# Тема: Филогенетика (MSA, NJ, bootstrap, ML, Bayesian)
# Документация: docs/05_phylogenetics.md
# ============================================

set -e

echo "🌳 Запуск пайплайна филогенетического анализа..."

# Проверка входных данных
INPUT_SEQ="${1:-sequences.fasta}"
OUTPUT_DIR="${2:-phylo_output}"
ANALYSIS_TYPE="${3:-all}"

if [ ! -f "$INPUT_SEQ" ]; then
    echo "❌ Ошибка: Файл последовательностей не найден: $INPUT_SEQ"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "📁 Входные последовательности: $INPUT_SEQ"
echo "📂 Выходная директория: $OUTPUT_DIR"
echo "🔬 Тип анализа: $ANALYSIS_TYPE"

# Шаг 1: Контроль качества и выравнивание
echo ""
echo "🔍 Шаг 1: Множественное выравнивание (MSA)..."
if command -v mafft &> /dev/null; then
    mafft --auto "$INPUT_SEQ" > "$OUTPUT_DIR/alignment.fasta" || echo "⚠️ MAFFT не удался..."
elif command -v muscle &> /dev/null; then
    muscle -in "$INPUT_SEQ" -out "$OUTPUT_DIR/alignment.fasta" || echo "⚠️ MUSCLE не удался..."
elif command -v clustalo &> /dev/null; then
    clustalo -i "$INPUT_SEQ" -o "$OUTPUT_DIR/alignment.fasta" --force || echo "⚠️ Clustal Omega не удался..."
else
    echo "⚠️ Инструменты MSA не найдены, пропускаем выравнивание..."
fi

# Шаг 2: Обрезка выравнивания (trimming)
echo ""
echo "✂️ Шаг 2: Обрезка выравнивания..."
if [ -f "$OUTPUT_DIR/alignment.fasta" ] && command -v trimal &> /dev/null; then
    trimal -in "$OUTPUT_DIR/alignment.fasta" -out "$OUTPUT_DIR/alignment_trimmed.fasta" -automated1 || echo "⚠️ trimAl не удался..."
elif [ -f "$OUTPUT_DIR/alignment.fasta" ] && command -v Gblocks &> /dev/null; then
    Gblocks "$OUTPUT_DIR/alignment.fasta" -t=d -b4=5 -b5=h || echo "⚠️ Gblocks не удался..."
    if [ -f "$OUTPUT_DIR/alignment.fasta-gb" ]; then
        mv "$OUTPUT_DIR/alignment.fasta-gb" "$OUTPUT_DIR/alignment_trimmed.fasta"
    fi
else
    echo "⚠️ Инструменты trimming не найдены, используем исходное выравнивание..."
    if [ -f "$OUTPUT_DIR/alignment.fasta" ]; then
        cp "$OUTPUT_DIR/alignment.fasta" "$OUTPUT_DIR/alignment_trimmed.fasta"
    fi
fi

# Шаг 3: Построение дерева методом Neighbor-Joining (NJ)
echo ""
echo "🌿 Шаг 3: Дерево Neighbor-Joining (NJ)..."
if [ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ] && command -v fasttree &> /dev/null; then
    FastTree -nt "$OUTPUT_DIR/alignment_trimmed.fasta" > "$OUTPUT_DIR/tree_nj.newick" || echo "⚠️ FastTree NJ не удался..."
elif [ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ] && command -v neighbor &> /dev/null; then
    # PHYLIP neighbor требует конвертации
    echo "⚠️ PHYLIP формат требует дополнительной конвертации..."
else
    echo "⚠️ Инструменты NJ не найдены, пропускаем..."
fi

# Шаг 4: Maximum Likelihood (ML) дерево с bootstrap
echo ""
echo "📊 Шаг 4: Maximum Likelihood (ML) с bootstrap..."
if [ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ] && command -v iqtree &> /dev/null; then
    iqtree -s "$OUTPUT_DIR/alignment_trimmed.fasta" \
           -m TEST \
           -B 1000 \
           -T AUTO \
           -pre "$OUTPUT_DIR/ml_tree" || echo "⚠️ IQ-TREE не удался..."
elif [ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ] && command -v raxml &> /dev/null; then
    raxmlHPC -s "$OUTPUT_DIR/alignment_trimmed.fasta" \
             -n ml_tree \
             -m GTRGAMMA \
             -p 12345 \
             -b 12345 \
             -N 1000 || echo "⚠️ RAxML не удался..."
elif [ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ] && command -v fasttree &> /dev/null; then
    FastTree -nt -gtr -gamma "$OUTPUT_DIR/alignment_trimmed.fasta" > "$OUTPUT_DIR/tree_ml.newick" || echo "⚠️ FastTree ML не удался..."
else
    echo "⚠️ Инструменты ML не найдены, пропускаем..."
fi

# Шаг 5: Bayesian inference (MrBayes/BEAST)
echo ""
echo "🎲 Шаг 5: Bayesian inference..."
if [ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ] && command -v mb &> /dev/null; then
    # Создаём файл команд для MrBayes
    cat > "$OUTPUT_DIR/mrbayes_commands.nex" << EOF
#NEXUS
begin data;
    dimensions ntax=0 nchar=0;
    format datatype=dna;
    matrix
    ;
end;

begin mrbayes;
    lset nst=6 rates=invgamma;
    mcmc ngen=1000000 samplefreq=1000 printfreq=1000 nchains=4;
    sumt burnin=250;
end;
EOF
    
    echo "⚠️ MrBayes требует ручной настройки NEXUS файла..."
    echo "   Шаблон создан: $OUTPUT_DIR/mrbayes_commands.nex"
elif command -v beast &> /dev/null; then
    echo "⚠️ BEAST требует XML конфигурации через BEAUti..."
else
    echo "⚠️ Bayesian инструменты не найдены, пропускаем..."
fi

# Шаг 6: Визуализация и аннотация дерева
echo ""
echo "🎨 Шаг 6: Подготовка к визуализации..."
if command -v figtree &> /dev/null; then
    echo "✅ FigTree доступен для визуализации"
elif command -v itol &> /dev/null; then
    echo "✅ iTOL CLI доступен"
else
    echo "ℹ️  Для визуализации рекомендуется:"
    echo "   - FigTree (desktop)"
    echo "   - iTOL (online: https://itol.embl.de)"
    echo "   - ggtree (R пакет)"
fi

# Создание файла для загрузки в iTOL
if [ -f "$OUTPUT_DIR/ml_tree.treefile" ]; then
    cp "$OUTPUT_DIR/ml_tree.treefile" "$OUTPUT_DIR/tree_for_itol.newick"
elif [ -f "$OUTPUT_DIR/tree_ml.newick" ]; then
    cp "$OUTPUT_DIR/tree_ml.newick" "$OUTPUT_DIR/tree_for_itol.newick"
elif [ -f "$OUTPUT_DIR/tree_nj.newick" ]; then
    cp "$OUTPUT_DIR/tree_nj.newick" "$OUTPUT_DIR/tree_for_itol.newick"
fi

# Шаг 7: Генерация отчёта
echo ""
echo "📊 Шаг 7: Генерация итогового отчёта..."
cat > "$OUTPUT_DIR/phylo_report.txt" << EOF
============================================
🌳 PHYLOGENETICS ANALYSIS REPORT
============================================
Input sequences: $INPUT_SEQ
Output directory: $OUTPUT_DIR
Analysis type: $ANALYSIS_TYPE
Date: $(date)

PIPELINE STATUS:
----------------
✅ MSA Alignment: $([ -f "$OUTPUT_DIR/alignment.fasta" ] && echo "Complete" || echo "Skipped")
✅ Trimming: $([ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ] && echo "Complete" || echo "Skipped")
✅ NJ Tree: $([ -f "$OUTPUT_DIR/tree_nj.newick" ] && echo "Complete" || echo "Skipped")
✅ ML Tree: $([ -f "$OUTPUT_DIR/ml_tree.treefile" ] || [ -f "$OUTPUT_DIR/tree_ml.newick" ] && echo "Complete" || echo "Skipped")
✅ Bayesian: $([ -f "$OUTPUT_DIR/mrbayes_commands.nex" ] && echo "Template created" || echo "Skipped")

ALIGNMENT STATISTICS:
---------------------
EOF

if [ -f "$OUTPUT_DIR/alignment.fasta" ]; then
    seq_count=$(grep -c "^>" "$OUTPUT_DIR/alignment.fasta")
    align_length=$(head -2 "$OUTPUT_DIR/alignment.fasta" | tail -1 | wc -c)
    echo "Number of sequences: $seq_count" >> "$OUTPUT_DIR/phylo_report.txt"
    echo "Alignment length: ~$((align_length - 1)) bp" >> "$OUTPUT_DIR/phylo_report.txt"
fi

if [ -f "$OUTPUT_DIR/alignment_trimmed.fasta" ]; then
    trimmed_length=$(head -2 "$OUTPUT_DIR/alignment_trimmed.fasta" | tail -1 | wc -c)
    echo "Trimmed length: ~$((trimmed_length - 1)) bp" >> "$OUTPUT_DIR/phylo_report.txt"
fi

cat >> "$OUTPUT_DIR/phylo_report.txt" << EOF

OUTPUT FILES:
-------------
$(ls -lh "$OUTPUT_DIR"/*.fasta "$OUTPUT_DIR"/*.newick "$OUTPUT_DIR"/*.treefile 2>/dev/null || echo "No tree files generated")

RECOMMENDED NEXT STEPS:
-----------------------
1. Visualize tree in FigTree or upload to iTOL
2. Root the tree using appropriate outgroup
3. Annotate clades with metadata
4. Calculate divergence times (if molecular clock needed)
5. Export publication-quality figures (SVG/PDF)

VISUALIZATION TOOLS:
--------------------
- FigTree: http://tree.bio.ed.ac.uk/software/figtree/
- iTOL: https://itol.embl.de
- ggtree (R): https://guangchuangyu.github.io/software/ggtree/
- TreeGraph 2: https://treegraph2.bioinfweb.info

============================================
📚 Дополнительная информация: docs/05_phylogenetics.md
============================================
EOF

echo ""
echo "✅ Филогенетический анализ завершён!"
echo "📄 Отчёт: $OUTPUT_DIR/phylo_report.txt"
echo "🌳 Дерево для визуализации: $OUTPUT_DIR/tree_for_itol.newick"
echo ""
echo "📚 Дополнительная информация: docs/05_phylogenetics.md"
