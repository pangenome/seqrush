#!/usr/bin/env python3
import sys

with open('d.paf') as f, open('d_no_self.paf', 'w') as out:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 6 and parts[0] != parts[5]:
            out.write(line)