"""Tool to automatically add the Fortran tests to the fpm.toml file.

To use, simply run this script from any directory. It will add the tests to the
fpm.toml file. The version control must be done manually.

python tools/dev/test_to_toml.py
"""

import toml
import re
from pathlib import Path

root_path = Path(__file__).parent.parent.parent
fpm_path = root_path / "fpm.toml"
tests_path = root_path / "test"


# Read fpm.toml
with open(fpm_path, "r") as f:
    fpm_text = f.read()

fpm = toml.loads(fpm_text)


# Discover all tests in the filesystem
found_tests = {}
for f90 in sorted(tests_path.rglob("test_*.f90")):
    content = f90.read_text()
    if not (
        re.search(r"^\s*program\b", content, re.IGNORECASE | re.MULTILINE)
        and re.search(r"^\s*end\s+program\b", content, re.IGNORECASE | re.MULTILINE)
    ):
        continue
    name = f90.stem.removeprefix("test_")
    found_tests[name] = {
        "name": name,
        "source-dir": f90.parent.relative_to(root_path).as_posix(),
        "main": f90.name,
    }

# Calculate differences
existing_tests = {t["name"]: t for t in fpm.get("test", [])}
found_names = set(found_tests.keys())
existing_names = set(existing_tests.keys())

to_add = found_names - existing_names
to_remove = existing_names - found_names

# Tests with same name but different source-dir or main (moved/renamed file)
to_update = {
    name
    for name in found_names & existing_names
    if found_tests[name]["source-dir"] != existing_tests[name]["source-dir"]
    or found_tests[name]["main"] != existing_tests[name]["main"]
}

if not to_add and not to_remove and not to_update:
    print("fpm.toml is already up to date.")
    exit(0)


# Remove [[test]] blocks that no longer exist
for name in to_remove:
    fpm_text = re.sub(
        rf'\[\[test\]\][^\[]*?name\s*=\s*"{re.escape(name)}".*?(?=\[\[|\Z)',
        "",
        fpm_text,
        flags=re.DOTALL,
    )


# Serialize the new blocks
def test_block(t: dict) -> str:
    return (
        f"[[test]]\n"
        f'name       = "{t["name"]}"\n'
        f'source-dir = "{t["source-dir"]}"\n'
        f'main       = "{t["main"]}"\n'
    )


# Merge existing and new tests, sorted alphabetically
all_tests = {t["name"]: t for t in fpm.get("test", [])}
all_tests.update({name: found_tests[name] for name in to_add | to_update})
for name in to_remove:
    all_tests.pop(name, None)

new_test_section = "\n".join(
    test_block(t) for t in sorted(all_tests.values(), key=lambda t: t["name"])
)


# Replace or insert the [[test]] section in the text
if re.search(r"\[\[test\]\]", fpm_text):
    # Remove all existing [[test]] blocks and append the new section
    fpm_text = re.sub(
        r"(\[\[test\]\].*?)(?=\[\[(?!test)|\Z)",
        "",
        fpm_text,
        flags=re.DOTALL,
    )
    fpm_text = fpm_text.rstrip() + "\n\n" + new_test_section
else:
    # No [[test]] blocks found: append at the end
    fpm_text = fpm_text.rstrip() + "\n\n" + new_test_section + "\n"


# Write
with open(fpm_path, "w") as f:
    f.write(fpm_text)

if to_add:
    print(f"Added {len(to_add)} test(s):")
    for name in sorted(to_add):
        print(
            f"  + {name}  ({found_tests[name]['source-dir']}/{found_tests[name]['main']})"
        )

if to_remove:
    print(f"Removed {len(to_remove)} test(s):")
    for name in sorted(to_remove):
        print(f"  - {name}")

if to_update:
    print(f"Updated {len(to_update)} test(s):")
    for name in sorted(to_update):
        print(
            f"  ~ {name}  ({found_tests[name]['source-dir']}/{found_tests[name]['main']})"
        )
