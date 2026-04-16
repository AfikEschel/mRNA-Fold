#!/usr/bin/env python3
"""Test quartet preprocessing functions."""

import sys
sys.path.insert(0, 'src')

from mrnafold.quartets import generate_quartets, find_conflicting_quartets, find_stackable_quartets

def main():
    seq = 'AUCUGCAUGGCCAAGAGGGUUA'
    quartets = generate_quartets(seq)
    print(f'Q: {len(quartets)} valid quartets')

    conflicts = find_conflicting_quartets(quartets)
    print(f'Total conflicting pairs: {len(conflicts)}')

    stackable_dict = find_stackable_quartets(quartets)
    total_stackable = sum(len(stackable) for stackable in stackable_dict.values())
    print(f'Total stackable relationships: {total_stackable}')

    if quartets:
        q0 = quartets[0]
        q0_conflicts = [pair for pair in conflicts if q0 in pair]
        q0_stackable = stackable_dict.get(q0, [])
        print(f'QC({q0}): {len(q0_conflicts)} conflicting quartets')
        print(f'QS({q0}): {len(q0_stackable)} stackable quartets')

if __name__ == '__main__':
    main()