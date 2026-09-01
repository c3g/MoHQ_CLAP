#!/usr/bin/env python3
"""check_r_vectors.py -- catch the comma errors a bracket count cannot.

    python3 bin/check_r_vectors.py bin/mohq_common.R

WHY THIS EXISTS
---------------
A missing comma between two string literals in an R vector,

    "^FBXW10"
    "^PRAMEF",

is a syntax error that balanced-bracket checks and "extract all the quoted
strings" checks both pass happily -- both were run, both said fine, and the
file then failed inside a Nextflow task with

    mohq_common.R:465:3: unexpected string constant

An adjacent-comma error costs a whole scheduler round-trip to discover. This
checks the one structure that keeps breaking: top-level `name <- c(...)`
character vectors, verifying their contents are comma-separated literals.

Not an R parser. It only checks what has actually gone wrong here.
"""
import re
import sys

def strip_comments(s: str) -> str:
    """Remove # comments, respecting BOTH quote styles.

    R strings may be single- or double-quoted, and a single-quoted string
    frequently contains double quotes -- regexes like

        sub('.*gene_name "([^"]+)".*', "\\\\1", attr)

    hold three of them. Tracking only `"` desynchronises on the first such
    line and every vector after it is misread, which produced a confident
    MISSING COMMA report against a perfectly valid palette. A checker that
    cries wolf is worse than no checker, so it tracks the opening quote
    character and only that character closes the string.
    """
    # Line by line, quote state reset per line. R string literals do not span
    # lines in this codebase, and scanning per line means one malformed line
    # cannot desynchronise the rest of the file -- which is exactly how the
    # previous version turned one unusual regex into a false report three
    # hundred lines later.
    kept = []
    for line in s.split("\n"):
        quote, esc, cut = None, False, len(line)
        for i, ch in enumerate(line):
            if esc:
                esc = False
            elif ch == "\\" and quote:
                esc = True
            elif quote:
                if ch == quote:
                    quote = None
            elif ch in ('"', "'"):
                quote = ch
            elif ch == "#":
                cut = i
                break
        kept.append(line[:cut])
    return "\n".join(kept)


def blocks(src: str):
    """Yield (name, body, start_line) for each top-level `name <- c(...)`."""
    for m in re.finditer(r'^([A-Za-z_.][\w.]*)\s*<-\s*c\(', src, re.M):
        name, i = m.group(1), m.end() - 1
        depth, j = 0, i
        while j < len(src):
            if src[j] == "(":
                depth += 1
            elif src[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        yield name, src[i + 1:j], src[:m.start()].count("\n") + 1


ELEMENT = re.compile(r'\s*(?:[\w.]+\s*=\s*)?"(?:[^"\\]|\\.)*"\s*')

bad = 0
for path in sys.argv[1:]:
    src = strip_comments(open(path, encoding="utf-8").read())
    for name, body, line in blocks(src):
        if '"' not in body:
            continue                     # not a character vector
        pos, n, first_line = 0, 0, line
        while pos < len(body):
            m = ELEMENT.match(body, pos)
            if not m:
                snippet = body[max(0, pos - 40):pos + 40].strip().replace("\n", " ")
                print(f"{path}:{line} [{name}] element {n + 1} is not a plain "
                      f"string literal, or a comma is missing before it:\n"
                      f"    ...{snippet}...")
                bad += 1
                break
            n += 1
            pos = m.end()
            if pos < len(body):
                if body[pos] != ",":
                    ctx = body[max(0, pos - 60):pos + 40].strip().replace("\n", " ")
                    print(f"{path}:{line} [{name}] MISSING COMMA after element "
                          f"{n}:\n    ...{ctx}...")
                    bad += 1
                    break
                pos += 1
                if ELEMENT.match(body, pos) is None and body[pos:].strip():
                    if body[pos:].lstrip().startswith(","):
                        print(f"{path}:{line} [{name}] DOUBLE COMMA after "
                              f"element {n}")
                        bad += 1
                        break
        else:
            print(f"  ok  {name}: {n} element(s)")

sys.exit(1 if bad else 0)
