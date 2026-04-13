# 🌳 Филогенетика (Phylogenetics)

## Обзор

Филогенетический анализ — изучение эволюционных взаимоотношений между организмами на основе генетических данных.

## Ключевые понятия

### MSA (Multiple Sequence Alignment)
Множественное выравнивание последовательностей для выявления гомологичных позиций.

**Популярные инструменты:**
- MAFFT — быстрое и точное выравнивание
- MUSCLE — баланс скорости и точности
- Clustal Omega — для больших наборов данных

### NJ (Neighbor-Joining)
Дистанционный метод построения деревьев.

**Алгоритм:**
1. Вычислить матрицу попарных расстояний
2. Найти ближайших соседей
3. Объединить в узел
4. Повторять до полного дерева

**Формула расстояния NJ:**
```
D(i,j) = d(i,j) - (r(i) + r(j)) / (N-2)
```
где:
- `d(i,j)` — исходное расстояние между таксонами i и j
- `r(i) = Σ d(i,k)` — сумма расстояний от i до всех других таксонов
- `N` — число таксонов

### Bootstrap
Метод оценки поддержки ветвей дерева.

**Принцип:**
1. Создать N псевдовыборок (обычно 100-1000)
2. Построить дерево для каждой выборки
3. Подсчитать, как часто каждая клада встречается

**Интерпретация:**
| Bootstrap | Надёжность |
|-----------|------------|
| ≥95% | Очень высокая |
| 70-94% | Высокая |
| 50-69% | Умеренная |
| <50% | Низкая |

### ML (Maximum Likelihood)
Статистический метод поиска дерева с максимальной вероятностью данных.

**Формула правдоподобия:**
```
L(Tree) = P(Data | Tree, Model)
```

**Популярные инструменты:**
- RAxML — быстрый ML анализ
- IQ-TREE — современный ML с подбором моделей
- PhyML — классический ML инструмент

### Bayesian Inference
Байесовский подход с использованием MCMC (Markov Chain Monte Carlo).

**Формула Байеса:**
```
P(Tree | Data) = [P(Data | Tree) × P(Tree)] / P(Data)
```

где:
- `P(Tree | Data)` — апостериорная вероятность дерева
- `P(Data | Tree)` — правдоподобие
- `P(Tree)` — априорная вероятность
- `P(Data)` — маргинальная вероятность данных

**Инструменты:**
- MrBayes — классический байесовский анализ
- BEAST — филогенетика с датировкой

## Модели эволюции

### Нуклеотидные модели

| Модель | Описание | Параметры |
|--------|----------|-----------|
| JC69 | Jukes-Cantor | 0 |
| K80 | Kimura 2-parameter | 1 (transition/transversion) |
| HKY85 | Hasegawa-Kishino-Yano | 1 + частоты оснований |
| GTR | General Time Reversible | 5 + частоты оснований |

### Формула поправки Джукса-Кантора
```
d = -(3/4) × ln(1 - (4/3) × p)
```
где `p` — доля различающихся сайтов.

## Команды

```bash
# Множественное выравнивание с MAFFT
mafft --auto input.fasta > alignment.fasta

# Построение дерева с IQ-TREE (ML + bootstrap)
iqtree \
  -s alignment.fasta \
  -m MFP \
  -bb 1000 \
  -nt AUTO

# Построение NJ дерева с FastME
fastme -i alignment.fasta -o tree.nwk

# Байесовский анализ с MrBayes
mb << EOF
execute alignment.fasta
lset nst=6 rates=invgamma
mcmc ngen=1000000 samplefreq=1000
sumt
EOF
```

## Визуализация

```bash
# Конвертация форматов с Newick Utilities
nw_reroot tree.nwk outgroup > rerooted.nwk

# Просмотр дерева
figtree tree.nwk
```

## Ссылки

- [IQ-TREE](http://www.iqtree.org/)
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [MrBayes](http://nbisweden.github.io/MrBayes/)
- [BEAST](http://beast.community/)
- [iTOL](https://itol.embl.de/) — онлайн визуализация деревьев
