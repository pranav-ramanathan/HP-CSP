"""
HP Sequence Parser (Compressed Notation)

Parses compressed HP sequences supporting ASCII digits and Unicode subscript digits.
Accepts whitespace/newlines/commas and ignores punctuation by default.
Optionally enforces strict mode.
"""

from typing import Tuple, List

SUBSCRIPT_MAP = {
    '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4',
    '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
}


class ParseError(ValueError):
    """Raised for invalid compressed HP sequences."""


def _is_ascii_digit(ch: str) -> bool:
    return '0' <= ch <= '9'


def _is_sub_digit(ch: str) -> bool:
    return ch in SUBSCRIPT_MAP


def expand_hp_sequence(raw: str, strict: bool = False) -> Tuple[str, List[str]]:
    """
    Expand a compressed HP sequence string.

    Grammar:
      sequence := item+
      item     := letter count?
      letter   := 'H' | 'P'
      count    := ascii_digits | unicode_subscript_digits

    Rules:
      - Mixed count styles (e.g., H₂3) are an error.
      - Zero repeats (H0 or H₀) are not allowed.
      - Unknown characters are ignored with a warning by default; in strict mode, error.

    Args:
        raw: Input string (may include whitespace/newlines/commas/punctuation)
        strict: If True, error on any char outside [H P 0-9 subscripts whitespace]

    Returns:
        (expanded_sequence, warnings)
    """
    s = raw
    i = 0
    n = len(s)
    out: List[str] = []
    warnings: List[str] = []

    def err(msg: str) -> None:
        raise ParseError(msg)

    while i < n:
        ch = s[i]

        # Skip whitespace and common separators by default
        if ch.isspace() or ch in {',', ';', '.', '-', '_', '(', ')', '[', ']', '{', '}'}:
            i += 1
            continue

        if ch == 'H' or ch == 'P':
            letter = ch
            i += 1
            # Parse optional count: ascii or subscript, but not mixed
            count_str = ''
            count_style = None  # 'ascii' or 'sub'
            while i < n:
                c2 = s[i]
                if _is_ascii_digit(c2):
                    if count_style is None:
                        count_style = 'ascii'
                    elif count_style != 'ascii':
                        err("Mixed digit styles in count (ASCII + subscript)")
                    count_str += c2
                    i += 1
                    continue
                if _is_sub_digit(c2):
                    if count_style is None:
                        count_style = 'sub'
                    elif count_style != 'sub':
                        err("Mixed digit styles in count (ASCII + subscript)")
                    count_str += SUBSCRIPT_MAP[c2]
                    i += 1
                    continue
                break

            if count_str == '':
                count = 1
            else:
                # Leading zeros are allowed but zero itself is not
                try:
                    count = int(count_str)
                except ValueError:
                    err(f"Invalid count after {letter}: '{count_str}'")
                if count == 0:
                    err("Zero repeats not allowed")
            out.append(letter * count)
            continue

        # Standalone digits (without preceding H/P)
        if _is_ascii_digit(ch) or _is_sub_digit(ch):
            if strict:
                err("Count without preceding letter (H/P)")
            warnings.append("Ignoring standalone digits without preceding H/P")
            i += 1
            continue

        # Unknown characters
        if strict and not ch.isspace():
            err(f"Illegal character in strict mode: '{ch}'")
        warnings.append(f"Ignoring unknown character: '{ch}'")
        i += 1

    expanded = ''.join(out)
    if len(expanded) == 0:
        err("Empty sequence.")
    # Validate alphabet
    if any(c not in ('H', 'P') for c in expanded):
        err("Expanded sequence contains invalid letters")
    return expanded, warnings


def hydrophobic_fraction(expanded: str) -> float:
    """Return fraction of H in expanded sequence (0..1)."""
    if not expanded:
        return 0.0
    return expanded.count('H') / len(expanded)


