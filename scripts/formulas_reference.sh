#!/bin/bash
# ============================================
# 📐 Formulas Reference & Calculator
# Тема: Сводка формул с пояснениями переменных
# Документация: docs/formulas_reference.md
# ============================================

set -e

echo "📐 Формулы и калькулятор биоинформатики..."

OUTPUT_DIR="${1:-formulas_output}"
mkdir -p "$OUTPUT_DIR"

cat > "$OUTPUT_DIR/formulas_reference.txt" << 'EOF'
============================================
📐 BIOINFORMATICS FORMULAS REFERENCE
============================================

1. E-VALUE (Expect Value)
-------------------------
Формула: E = K * m * n * e^(-λ*S)

Где:
  E — ожидаемое число совпадений со score ≥ S
  K, λ — параметры статистики экстремальных значений (для матрицы замен)
  m — длина запроса (query length)
  n — размер базы данных
  S — raw alignment score

Интерпретация:
  - E < 1e-5: высокозначимое совпадение
  - E < 0.01: значимое совпадение
  - E > 10: вероятно случайное совпадение

2. BIT SCORE
------------
Формула: S' = (λ*S - ln(K)) / ln(2)

Где:
  S' — bit score (нормализованный score)
  S — raw score
  λ, K — параметры матрицы замен

Преимущество: не зависит от размера базы данных

3. SEQUENCE IDENTITY
--------------------
Формула: Identity (%) = (M / L) * 100

Где:
  M — число идентичных позиций
  L — длина выравнивания (без gaps) или общая длина

4. GC CONTENT
-------------
Формула: GC% = ((G + C) / (A + T + G + C)) * 100

Где:
  A, T, G, C — количество соответствующих нуклеотидов

5. MELTING TEMPERATURE (Tm)
---------------------------
Базовая формула (Wallace rule):
Tm = 2*(A+T) + 4*(G+C)

Улучшенная формула (Nearest Neighbor):
Tm = ΔH / (ΔS + R*ln(C/4)) - 273.15

Где:
  ΔH — энтальпия (kcal/mol)
  ΔS — энтропия (cal/mol*K)
  R — газовая постоянная (1.987 cal/mol*K)
  C — концентрация олигонуклеотида

6. COVERAGE (DEPTH)
-------------------
Формула: Coverage = (N * L) / G

Где:
  N — число ридов
  L — средняя длина рида
  G — размер генома

7. PHRED QUALITY SCORE
----------------------
Формула: Q = -10 * log10(P)

Где:
  Q — Phred score
  P — вероятность ошибки базового вызова

Обратно: P = 10^(-Q/10)

Примеры:
  Q20 = 1% ошибка (точность 99%)
  Q30 = 0.1% ошибка (точность 99.9%)
  Q40 = 0.01% ошибка (точность 99.99%)

8. ALLELE FREQUENCY
-------------------
Формула: AF = (count_alt) / (count_ref + count_alt)

Где:
  count_alt — число аллелей варианта
  count_ref — число референсных аллелей

9. HARDY-WEINBERG EQUILIBRIUM
-----------------------------
Формула: p² + 2pq + q² = 1

Где:
  p — частота доминантного аллеля
  q — частота рецессивного аллеля
  p² — частота гомозигот AA
  2pq — частота гетерозигот Aa
  q² — частота гомозигот aa

10. FST (FIXATION INDEX)
------------------------
Формула: Fst = (Ht - Hs) / Ht

Где:
  Ht — ожидаемая гетерозиготность в общей популяции
  Hs — средняя ожидаемая гетерозиготность в субпопуляциях

Интерпретация:
  0-0.05: низкая дифференциация
  0.05-0.15: умеренная
  0.15-0.25: высокая
  >0.25: очень высокая

11. JUKES-CANTOR DISTANCE
-------------------------
Формула: d = -(3/4) * ln(1 - (4/3)*p)

Где:
  d — эволюционное расстояние (замен на сайт)
  p — доля несовпадающих сайтов

12. KIMURA 2-PARAMETER DISTANCE
-------------------------------
Формула: d = -(1/2)*ln(1-2P-Q) - (1/4)*ln(1-2Q)

Где:
  P — доля транзиций
  Q — доля трансверсий

13. TAJIMA'S D
--------------
Формула: D = (π - θw) / √(Var(π - θw))

Где:
  π — среднее число попарных различий
  θw — оценка Вараттера (на основе числа сегрегирующих сайтов)

Интерпретация:
  D ≈ 0: нейтральная эволюция
  D < 0: избыток редких аллелей (селективное сканирование, expansion)
  D > 0: избыток частых аллелей (balancing selection, bottleneck)

14. LIKELIHOOD RATIO TEST (LRT)
-------------------------------
Формула: LRT = 2 * (ln(L1) - ln(L0))

Где:
  L1 — likelihood альтернативной модели
  L0 — likelihood нулевой модели

Распределение: χ² с k степенями свободы (k = разница в параметрах)

15. BAYESIAN POSTERIOR PROBABILITY
----------------------------------
Формула: P(H|D) = P(D|H) * P(H) / P(D)

Где:
  P(H|D) — апостериорная вероятность гипотезы
  P(D|H) — likelihood данных при данной гипотезе
  P(H) — априорная вероятность
  P(D) — маргинальная вероятность данных

============================================
🧮 QUICK CALCULATOR
============================================
EOF

# Функция для расчёта GC-content
calculate_gc() {
    local seq="$1"
    local g=$(echo "$seq" | grep -o "[Gg]" | wc -l)
    local c=$(echo "$seq" | grep -o "[Cc]" | wc -l)
    local a=$(echo "$seq" | grep -o "[Aa]" | wc -l)
    local t=$(echo "$seq" | grep -o "[Tt]" | wc -l)
    local total=$((a + t + g + c))
    
    if [ $total -gt 0 ]; then
        local gc_percent=$(awk "BEGIN {printf \"%.2f\", (($g + $c) / $total) * 100}")
        echo "GC-content: $gc_percent%"
    else
        echo "Некорректная последовательность"
    fi
}

# Функция для расчёта Tm (Wallace rule)
calculate_tm() {
    local seq="$1"
    local at=$(echo "$seq" | grep -o "[AaTt]" | wc -l)
    local gc=$(echo "$seq" | grep -o "[GgCc]" | wc -l)
    
    local tm=$((2 * at + 4 * gc))
    echo "Tm (Wallace): ${tm}°C"
}

# Функция для расчёта Phred quality
calculate_phred() {
    local q="$1"
    local p=$(awk "BEGIN {printf \"%.10f\", 10^(-$q/10)}")
    local accuracy=$(awk "BEGIN {printf \"%.4f\", (1 - 10^(-$q/10)) * 100}")
    echo "Q$q → Error probability: $p, Accuracy: ${accuracy}%"
}

# Примеры использования
cat >> "$OUTPUT_DIR/formulas_reference.txt" << EOF

ПРИМЕРЫ РАСЧЁТОВ:
-----------------

1. GC-Content для последовательности "ATGCGCAT":
$(calculate_gc "ATGCGCAT")

2. Tm для праймера "ATGCGCATGC":
$(calculate_tm "ATGCGCATGC")

3. Phred quality scores:
$(calculate_phred 20)
$(calculate_phred 30)
$(calculate_phred 40)

============================================
📊 USAGE EXAMPLES
============================================

# Расчёт GC-content
echo "ATGCGCATGCAT" | ./scripts/formulas_reference.sh gc

# Расчёт Tm праймера
echo "ATGCGCATGC" | ./scripts/formulas_reference.sh tm

# Конвертация Phred score
./scripts/formulas_reference.sh phred 30

============================================
📚 Дополнительная информация: docs/formulas_reference.md
============================================
EOF

echo ""
echo "✅ Справочник формул создан!"
echo "📄 Файл: $OUTPUT_DIR/formulas_reference.txt"
echo ""
echo "📚 Дополнительная информация: docs/formulas_reference.md"
