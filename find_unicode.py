#!/usr/bin/env python3
import os
for root, dirs, files in os.walk('.'):
    if '.git' in root:
        continue
    for f in files:
        if f.endswith('.nf') or f.endswith('.config') or f.endswith('.py'):
            path = os.path.join(root, f)
            with open(path, 'rb') as fh:
                for i, line in enumerate(fh, 1):
                    for j, b in enumerate(line):
                        if b > 127:
                            context = line[max(0,j-10):j+10]
                            print(f"{path}:{i}:{j} byte=0x{b:02x} context={context}")
