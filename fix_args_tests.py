#!/usr/bin/env python3
import os
import re

# List of files to fix
files_to_fix = [
    "tests/bidirected_inversion_test.rs",
    "tests/bidirected_tests.rs",
    "tests/integration_tests.rs",
    "tests/mathematical_bidirected_tests.rs",
    "tests/test_complex_structural_variations.rs",
    "tests/test_edge_traversal.rs",
    "tests/test_inversion_detection.rs",
    "tests/test_programmatic_variations.rs",
    "tests/test_rc_node_grouping.rs",
    "tests/test_real_inversion.rs",
    "tests/test_topological_sort.rs",
    "src/reverse_complement_tests.rs",
    "src/lib.rs"
]

# New fields to add
new_fields = """        groom: false,
        iterative_groom: None,
        odgi_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,"""

for filepath in files_to_fix:
    if not os.path.exists(filepath):
        print(f"Skipping {filepath} - does not exist")
        continue

    with open(filepath, 'r') as f:
        content = f.read()

    # Pattern to match Args struct initialization ending with }
    # Look for patterns like "no_sort: false,\n    };" or "no_sort: true,\n    };"
    pattern = r'(no_sort: (?:true|false)),(\s*\n\s*)\};'

    if re.search(pattern, content):
        # Replace with the fields added
        replacement = r'\1,\2' + new_fields + r'\2};'
        new_content = re.sub(pattern, replacement, content)

        with open(filepath, 'w') as f:
            f.write(new_content)
        print(f"Fixed {filepath}")
    else:
        print(f"Pattern not found in {filepath}")

print("Done")