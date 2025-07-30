#!/usr/bin/env python3
import os
import re

def needs_fix(content, line_num):
    """Check if this validate_paf line needs the missing fields added."""
    lines = content.split('\n')
    # Look at the next few lines to see if paf: and seqwish_style: are already there
    for i in range(line_num, min(line_num + 5, len(lines))):
        if 'paf:' in lines[i]:
            return False
    return True

def fix_file(filepath):
    """Fix a single file by adding missing fields after validate_paf lines."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    fixed_lines = []
    i = 0
    
    while i < len(lines):
        line = lines[i]
        fixed_lines.append(line)
        
        # Check if this line has validate_paf: and ends with a comma
        if 'validate_paf:' in line and line.strip().endswith(','):
            # Check if we need to add the fields
            if needs_fix(content, i):
                # Add the missing fields
                indent = len(line) - len(line.lstrip())
                fixed_lines.append(' ' * indent + 'paf: None,')
                fixed_lines.append(' ' * indent + 'seqwish_style: false,')
        
        i += 1
    
    return '\n'.join(fixed_lines)

# Process all files
files_to_fix = [
    'tests/test_programmatic_variations.rs',
    'tests/test_edge_traversal.rs', 
    'tests/test_inversion_detection.rs',
    'tests/mathematical_bidirected_tests.rs',
    'tests/test_complex_structural_variations.rs',
    'tests/integration_tests.rs',
    'tests/test_real_inversion.rs',
    'tests/bidirected_inversion_test.rs',
    'tests/bidirected_tests.rs',
    'src/lib.rs',
    'src/reverse_complement_tests.rs'
]

for filepath in files_to_fix:
    if os.path.exists(filepath):
        print(f"Fixing {filepath}")
        fixed_content = fix_file(filepath)
        with open(filepath, 'w') as f:
            f.write(fixed_content)

print("Done!")