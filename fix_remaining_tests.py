#!/usr/bin/env python3
import re
import glob

def fix_args_structs(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Pattern to match Args structs that end with seqwish_style: false, followed by whitespace and run_
    pattern = r'(let args = Args \{[\s\S]*?seqwish_style: false,)(\s*)(run_(?:inversion_aware_)?seqrush)'
    
    def replacer(match):
        # Add closing brace before the run_ call
        return match.group(1) + '\n    };\n    \n    ' + match.group(3)
    
    fixed = re.sub(pattern, replacer, content)
    
    # Also fix the pattern where it's just missing a closing brace before println
    pattern2 = r'(let args = Args \{[\s\S]*?seqwish_style: false,)(\s*)(println!)'
    fixed = re.sub(pattern2, lambda m: m.group(1) + '\n    };\n    \n    ' + m.group(3), fixed)
    
    if fixed != content:
        with open(filepath, 'w') as f:
            f.write(fixed)
        return True
    return False

# Fix specific files with issues
problem_files = [
    'tests/test_inversion_detection.rs',
    'tests/test_real_inversion.rs',
    'tests/integration_tests.rs'
]

for filepath in problem_files:
    if fix_args_structs(filepath):
        print(f"Fixed {filepath}")

print("Done!")