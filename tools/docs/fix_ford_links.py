#!/usr/bin/env python3
"""Repair known broken links in a FORD-generated site (FORD 7.0.13).

Run after ``ford ford.md`` and before publishing::

    python3 tools/docs/fix_ford_links.py doc/ford_site

Two FORD bugs are corrected in place:

1. Type-bound procedure bindings (``procedure :: name => impl``) are rendered
   in module and type pages with the absolute filesystem path of the output
   directory, e.g. ``href='/home/runner/.../doc/ford_site/proc/impl.html'``.
   FORD's ``macros.html`` joins ``tb.bindings`` without its ``relurl``
   filter. These hrefs are rewritten relative to the page.
2. Generic-interface pages link each module procedure to an anchor
   (``interface/abs.html#moduleprocedure-abshyperdual``) that FORD never
   writes. When the procedure has its own page (``proc/abshyperdual.html``),
   the link is pointed there. FORD's ``~N`` disambiguation suffix is dropped
   when no ``proc/<name>~N.html`` page exists.

The script is idempotent, uses only the standard library, and touches nothing
but ``href`` attribute values. Exit status is 0; a summary is printed.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

HREF = re.compile(r"""(?P<attr>\b(?:xlink:)?href=)(?P<q>["'])(?P<url>.*?)(?P=q)""")
MODPROC = re.compile(r"^(?P<path>[^#]*interface/[^#]+\.html)#moduleprocedure-(?P<name>[^#]+)$")


def fix_page(page: Path, root: Path, ids_cache: dict[Path, set[str]]) -> tuple[str, int, int]:
    text = page.read_text(encoding="utf-8")
    root_str = str(root)
    n_abs = 0
    n_modproc = 0

    def page_ids(path: Path) -> set[str]:
        if path not in ids_cache:
            try:
                ids_cache[path] = set(re.findall(r"""\bid=["']([^"']+)["']""",
                                                 path.read_text(encoding="utf-8")))
            except OSError:
                ids_cache[path] = set()
        return ids_cache[path]

    def repl(match: re.Match) -> str:
        nonlocal n_abs, n_modproc
        url = match.group("url")
        new = url
        # 1. Absolute output-directory paths -> relative.
        if url.startswith(root_str + os.sep):
            target, _, frag = url.partition("#")
            new = os.path.relpath(target, page.parent).replace(os.sep, "/")
            if frag:
                new += "#" + frag
            n_abs += 1
        # 2. Interface module-procedure anchors that do not exist.
        m = MODPROC.match(new)
        if m:
            iface = (page.parent / m.group("path")).resolve()
            anchor = "moduleprocedure-" + m.group("name")
            proc_page = root / "proc" / (m.group("name") + ".html")
            if not proc_page.is_file():
                # FORD suffixes "~N" to names it has already used; the
                # procedure page itself usually carries the plain name.
                base = re.sub(r"~\d+$", "", m.group("name"))
                proc_page = root / "proc" / (base + ".html")
            if anchor not in page_ids(iface) and proc_page.is_file():
                new = os.path.relpath(proc_page, page.parent).replace(os.sep, "/")
                n_modproc += 1
        if new == url:
            return match.group(0)
        return f"{match.group('attr')}{match.group('q')}{new}{match.group('q')}"

    return HREF.sub(repl, text), n_abs, n_modproc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("site", type=Path, help="FORD output directory, e.g. doc/ford_site")
    args = ap.parse_args()
    root = args.site.resolve()
    if not (root / "index.html").is_file():
        ap.error(f"{root} does not look like a FORD site (no index.html)")

    ids_cache: dict[Path, set[str]] = {}
    total_abs = total_modproc = changed = 0
    for page in sorted(root.rglob("*.html")):
        if (root / "page" / "python-api") in page.parents:
            continue  # Sphinx output, not FORD
        new_text, n_abs, n_modproc = fix_page(page, root, ids_cache)
        if n_abs or n_modproc:
            page.write_text(new_text, encoding="utf-8")
            changed += 1
            total_abs += n_abs
            total_modproc += n_modproc
    print(f"fix_ford_links: {total_abs} absolute binding link(s) made relative, "
          f"{total_modproc} interface anchor(s) redirected, {changed} page(s) changed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
