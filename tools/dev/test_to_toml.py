"""Tool to automatically add the Fortran tests to the fpm.toml file.

To use, simply run this scipt from any directory. It will add the tests to the
fpm.toml file. The version control must be done manually.
"""

import toml
from pathlib import Path


root_path = Path(__file__).parent.parent.parent
fpm_path = root_path / "fpm.toml"
tests_path = root_path / "tests"


# -- Read toml ----------------------------------------------------------------
with open(fpm_path, "r") as f:
    fpm = toml.load(f)


# -- Get all tests in filesystem ----------------------------------------------
found_tests = {}
for f90 in sorted(tests_path.rglob("test_*.f90")):
    name = f90.stem.removeprefix("test_")
    found_tests[name] = {
        "name": name,
        "source-dir": f90.parent.relative_to(root_path).as_posix(),
        "main": f90.name,
    }


# -- Get differences in toml and filesystem -----------------------------------
existing_names = {t["name"] for t in fpm.get("test", [])}
found_names = set(found_tests.keys())

to_add = found_names - existing_names
to_remove = existing_names - found_names


# -- Apply changes ------------------------------------------------------------
fpm.setdefault("test", [])

fpm["test"] = [t for t in fpm["test"] if t["name"] not in to_remove]
fpm["test"].extend(found_tests[name] for name in to_add)
fpm["test"].sort(key=lambda t: t["name"])


# -- Write changes ------------------------------------------------------------
with open(fpm_path, "w") as f:
    toml.dump(fpm, f)

if to_add:
    print(f"Added {len(to_add)} test(s):")
    for name in sorted(to_add):
        print(f"  + {name}  ({found_tests[name]['source-dir']}/{found_tests[name]['main']})")

if to_remove:
    print(f"Removed {len(to_remove)} test(s):")
    for name in sorted(to_remove):
        print(f"  - {name}")

if not to_add and not to_remove:
    print("fpm.toml is already up to date.")