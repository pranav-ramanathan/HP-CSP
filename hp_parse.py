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
    Expand a compressed HP sequence string with group repeats.

    Grammar:
      sequence := item+
      item     := atom count?
      atom     := letter | '(' sequence ')'
      letter   := 'H' | 'P'
      count    := ascii_digits | unicode_subscript_digits

    Rules:
      - Mixed digit styles within one count (e.g., H₂3) are an error.
      - Zero repeats (…0 or …₀) are not allowed.
      - Unknown characters are ignored with a warning by default; in strict mode, error.

    Args:
        raw: Input string (may include whitespace/newlines/commas/punctuation)
        strict: If True, error on any char outside [H P 0-9 subscripts whitespace]

    Returns:
        (expanded_sequence, warnings)
    """
    s = raw
    n = len(s)
    warnings: List[str] = []

    def err(msg: str) -> None:
        raise ParseError(msg)

    def read_count(i: int) -> Tuple[int, int]:
        count_str = ''
        count_style = None  # 'ascii' or 'sub'
        j = i
        while j < n:
            c = s[j]
            if _is_ascii_digit(c):
                if count_style is None:
                    count_style = 'ascii'
                elif count_style != 'ascii':
                    err("Mixed digit styles in count (ASCII + subscript)")
                count_str += c
                j += 1
                continue
            if _is_sub_digit(c):
                if count_style is None:
                    count_style = 'sub'
                elif count_style != 'sub':
                    err("Mixed digit styles in count (ASCII + subscript)")
                count_str += SUBSCRIPT_MAP[c]
                j += 1
                continue
            break
        if count_str == '':
            return 1, i
        try:
            val = int(count_str)
        except ValueError:
            err(f"Invalid count: '{count_str}'")
        if val == 0:
            err("Zero repeats not allowed")
        return val, j

    def parse_seq(i: int, end_char: str = None) -> Tuple[str, int]:
        out_parts: List[str] = []
        j = i
        while j < n:
            ch = s[j]
            # End of a parenthesized group
            if end_char is not None and ch == end_char:
                j += 1
                return ''.join(out_parts), j

            # Whitespace and common separators (parentheses handled explicitly)
            if ch.isspace() or ch in {',', ';', '.', '-', '_'}:
                j += 1
                continue

            # Single letters with optional count
            if ch == 'H' or ch == 'P':
                letter = ch
                j += 1
                cnt, j = read_count(j)
                out_parts.append(letter * cnt)
                continue

            # Parenthesized group
            if ch == '(':
                j += 1
                inner, j = parse_seq(j, end_char=')')
                cnt, j = read_count(j)
                out_parts.append(inner * cnt)
                continue

            # Unmatched ')'
            if ch == ')':
                if strict:
                    err("Unmatched ')'")
                warnings.append("Ignoring unmatched ')'")
                j += 1
                continue

            # Standalone digits (no preceding letter/group)
            if _is_ascii_digit(ch) or _is_sub_digit(ch):
                if strict:
                    err("Count without preceding letter/group")
                warnings.append("Ignoring standalone digits without preceding letter/group")
                j += 1
                continue

            # Unknown characters
            if strict and not ch.isspace():
                err(f"Illegal character in strict mode: '{ch}'")
            warnings.append(f"Ignoring unknown character: '{ch}'")
            j += 1

        # If we were inside a group and hit end-of-string, it's an error in strict mode
        if end_char is not None:
            if strict:
                err("Unclosed '(' group")
            warnings.append("Unclosed '(' group (ignored)")
        return ''.join(out_parts), j

    expanded, _ = parse_seq(0, end_char=None)
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


