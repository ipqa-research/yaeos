#!/usr/bin/env python3
"""Align FORD doc-comment blocks (`!!`, `!>`) with the code they document.

Only the leading spaces of comment-only lines change. Rules:
* `!>` (predoc) block: indentation of the next code line.
* `!!` block: indentation of the next code line, unless that line closes a
  scope (end/contains/else/case); then the previous line's indentation, plus
  3 when the previous line opens a scope (module, type, procedure, ...).
"""
import pathlib
import re
import sys

END = re.compile(r"\s*(end\b|contains\b|else\b|elseif\b|case\b)", re.I)
OPENER = re.compile(
    r"\s*("
    r"(module|submodule|program|block\s*data)\b(?!\s+procedure)"
    r"|(abstract\s+)?interface\b"
    r"|type\s*(,|::|\s+[a-z_])"
    r"|((pure|impure|elemental|recursive|module|real|integer|logical|complex|character|double\s+precision|type\s*\([^)]*\)|class\s*\([^)]*\))(\([^)]*\))?\s+)*(subroutine|function)\b"
    r")",
    re.I,
)


def ind(line):
    return len(line) - len(line.lstrip(" "))


def marker(line):
    s = line.lstrip(" ")
    if s.startswith("!!"):
        return "!!"
    if s.startswith("!>"):
        return "!>"
    return None


def fix_text(text):
    lines = text.split("\n")
    changed = 0
    i = 0
    while i < len(lines):
        m = marker(lines[i])
        if not m:
            i += 1
            continue
        j = i
        while j < len(lines) and marker(lines[j]) == m:
            j += 1
        k = j
        while k < len(lines) and not lines[k].strip():
            k += 1
        nxt = lines[k] if k < len(lines) else ""
        p = i - 1
        while p >= 0 and not lines[p].strip():
            p -= 1
        prev = lines[p] if p >= 0 else ""
        if m == "!>" or (nxt.strip() and not END.match(nxt) and not marker(nxt)):
            exp = ind(nxt) if nxt.strip() else None
        elif prev.strip():
            exp = ind(prev) + (3 if OPENER.match(prev) else 0)
        else:
            exp = None
        if exp is not None:
            for n in range(i, j):
                new = " " * exp + lines[n].lstrip(" ")
                if new != lines[n]:
                    lines[n] = new
                    changed += 1
        i = j
    return "\n".join(lines), changed


def main():
    write = "--write" in sys.argv
    roots = [a for a in sys.argv[1:] if not a.startswith("--")] or ["src"]
    total = files = 0
    for root in roots:
        for f in sorted(pathlib.Path(root).rglob("*.f90")):
            text = f.read_text()
            new, n = fix_text(text)
            if n:
                files += 1
                total += n
                if write:
                    f.write_text(new)
                else:
                    print(f"{f}: {n} line(s)")
    print(f"{'re-indented' if write else 'would re-indent'} {total} line(s) in {files} file(s)")


if __name__ == "__main__":
    main()
