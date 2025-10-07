#!/usr/bin/env python3
import os
import re

def fix_args_struct(content):
    """Add missing paf and seqwish_style fields to Args structs."""
    # Pattern to find Args struct initialization
    pattern = r'(let args = Args \{[^}]+validate_paf: (?:true|false),)\s*\};'
    
    def replacer(match):
        return match.group(1) + '\n        paf: None,\n        seqwish_style: false,\n    };'
    
    # Apply the replacement
    fixed = re.sub(pattern, replacer, content, flags=re.DOTALL)
    return fixed

# Fix test files
test_dirs = ['tests', 'src']
for test_dir in test_dirs:
    if not os.path.exists(test_dir):
        continue
        
    for filename in os.listdir(test_dir):
        if filename.endswith('.rs'):
            filepath = os.path.join(test_dir, filename)
            with open(filepath, 'r') as f:
                content = f.read()
            
            if 'Args {' in content and 'validate_paf:' in content:
                # Check if already has paf field
                if 'paf:' not in content:
                    print(f"Fixing {filepath}")
                    fixed_content = fix_args_struct(content)
                    with open(filepath, 'w') as f:
                        f.write(fixed_content)

print("Done!")