#!/usr/bin/env python3
"""
N50 Calculator - Вычисление метрик качества сборки генома
"""

import sys
from pathlib import Path

def parse_fasta(filepath):
    """Парсинг FASTA файла и возврат длин последовательностей."""
    lengths = []
    current_length = 0
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_length > 0:
                    lengths.append(current_length)
                current_length = 0
            else:
                current_length += len(line)
        if current_length > 0:
            lengths.append(current_length)
    
    return lengths

def calculate_n50(lengths):
    """Расчёт N50 и других метрик."""
    if not lengths:
        return None
    
    lengths.sort(reverse=True)
    total = sum(lengths)
    cumsum = 0
    n50 = 0
    n90 = 0
    l50 = 0
    l90 = 0
    
    for i, length in enumerate(lengths):
        cumsum += length
        if n50 == 0 and cumsum >= total * 0.5:
            n50 = length
            l50 = i + 1
        if n90 == 0 and cumsum >= total * 0.9:
            n90 = length
            l90 = i + 1
    
    return {
        'total_length': total,
        'num_contigs': len(lengths),
        'max_length': max(lengths),
        'min_length': min(lengths),
        'mean_length': total / len(lengths),
        'n50': n50,
        'l50': l50,
        'n90': n90,
        'l90': l90
    }

def main():
    if len(sys.argv) < 2:
        print("Использование: calculate_n50.py <файл.fasta>")
        sys.exit(1)
    
    filepath = Path(sys.argv[1])
    if not filepath.exists():
        print(f"❌ Файл не найден: {filepath}")
        sys.exit(1)
    
    print(f"📊 Анализ файла: {filepath}")
    lengths = parse_fasta(filepath)
    metrics = calculate_n50(lengths)
    
    if metrics:
        print("\n" + "="*50)
        print("📈 МЕТРИКИ СБОРКИ")
        print("="*50)
        print(f"   Общая длина:       {metrics['total_length']:,} bp")
        print(f"   Количество контигов: {metrics['num_contigs']}")
        print(f"   Макс. длина:       {metrics['max_length']:,} bp")
        print(f"   Мин. длина:        {metrics['min_length']:,} bp")
        print(f"   Средняя длина:     {metrics['mean_length']:,.0f} bp")
        print("-"*50)
        print(f"   🎯 N50:            {metrics['n50']:,} bp")
        print(f"   🎯 L50:            {metrics['l50']} контигов")
        print(f"   🎯 N90:            {metrics['n90']:,} bp")
        print(f"   🎯 L90:            {metrics['l90']} контигов")
        print("="*50)
    else:
        print("❌ Не удалось рассчитать метрики")
        sys.exit(1)

if __name__ == "__main__":
    main()
