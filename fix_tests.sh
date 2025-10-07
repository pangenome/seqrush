#!/bin/bash

# Find all test files with Args initialization
for file in tests/*.rs; do
    if grep -q "Args {" "$file"; then
        echo "Fixing $file"
        # Add the missing fields before the closing brace
        sed -i '/validate_paf: .*,$/a\        paf: None,\n        seqwish_style: false,' "$file"
    fi
done

# Also fix src/lib.rs if needed
if grep -q "Args {" "src/lib.rs"; then
    echo "Fixing src/lib.rs"
    sed -i '/validate_paf: .*,$/a\        paf: None,\n        seqwish_style: false,' "src/lib.rs"
fi

# Also fix src/reverse_complement_tests.rs if needed  
if grep -q "Args {" "src/reverse_complement_tests.rs"; then
    echo "Fixing src/reverse_complement_tests.rs"
    sed -i '/validate_paf: .*,$/a\        paf: None,\n        seqwish_style: false,' "src/reverse_complement_tests.rs"
fi

echo "Done!"