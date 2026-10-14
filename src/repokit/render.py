"""Token filling and the banner every generated file carries."""
from __future__ import annotations

import re

TOKEN = re.compile(r"\{\{\s*([A-Za-z][A-Za-z0-9_.-]*)\s*\}\}")
COMMAND = "make generate"


def tokens(text: str) -> list[str]:
    return TOKEN.findall(text)


def fill(text: str, mapping: dict[str, str], where: str) -> str:
    missing: list[str] = []

    def replace(match: re.Match) -> str:
        key = match.group(1)
        if key not in mapping:
            missing.append(key)
            return match.group(0)
        return mapping[key]

    out = TOKEN.sub(replace, text)
    if missing:
        raise KeyError(f"{where}: undefined tokens {sorted(set(missing))}")
    return out


def banner(source: str, style: str) -> str:
    words = f"GENERATED from {source} by `{COMMAND}`; edit the source, not this file."
    if style == "markdown":
        return f"<!-- {words} -->\n"
    if style == "hash":
        return f"# {words}\n"
    if style == "tex":
        return f"% {words}\n"
    raise ValueError(style)
