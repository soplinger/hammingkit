#!/usr/bin/env python3
"""hammingkit: offline Hamming and extended SECDED error-correcting codes.

Standard library only. No network, no accounts.
"""
from __future__ import annotations

import argparse
import random
import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path

MAX_R = 6          # largest base Hamming code is (63,57)
VERIFY_MAX_K = 11  # verify brute-forces 2^k data words


class HammingError(ValueError):
    """Bad input. The message is meant for the person at the terminal."""


@dataclass(frozen=True)
class Code:
    """A Hamming code with r parity bits, optionally extended with an overall parity bit."""

    r: int
    extended: bool = False

    @property
    def n(self) -> int:
        """Length of the base Hamming code (positions 1..n)."""
        return (1 << self.r) - 1

    @property
    def length(self) -> int:
        """Transmitted block length."""
        return self.n + (1 if self.extended else 0)

    @property
    def k(self) -> int:
        """Data bits per block."""
        return self.n - self.r

    @property
    def min_distance(self) -> int:
        return 4 if self.extended else 3

    @property
    def name(self) -> str:
        return f"({self.length},{self.k})"

    def parity_positions(self) -> list[int]:
        return [1 << i for i in range(self.r)]

    def data_positions(self) -> list[int]:
        return [p for p in range(1, self.n + 1) if p & (p - 1)]


@dataclass(frozen=True)
class Decoded:
    status: str                 # "ok", "corrected" or "detected"
    data: list[int] | None      # None when the block cannot be trusted
    error_position: int | None  # 1-based transmitted position that was flipped back


def parse_code(text: str) -> Code:
    """Turn 'N,K' into a code. Base codes are (2^r-1, ...), extended are (2^r, ...)."""
    try:
        n_text, k_text = text.split(",")
        n, k = int(n_text), int(k_text)
    except ValueError:
        raise HammingError(f"code must look like N,K such as 7,4 or 8,4, got {text!r}")
    for r in range(2, MAX_R + 1):
        if n == (1 << r) - 1 and k == n - r:
            return Code(r, False)
        if n == (1 << r) and k == (1 << r) - 1 - r:
            return Code(r, True)
    raise HammingError(
        f"unsupported code {n},{k}; supported: 3,1 4,1 7,4 8,4 15,11 16,11 31,26 32,26 ..."
    )


def parse_bits(text: str) -> list[int]:
    """Accept 0/1 strings, ignoring underscores and spaces."""
    cleaned = text.replace("_", "").replace(" ", "")
    if not cleaned or set(cleaned) - {"0", "1"}:
        raise HammingError(f"expected only 0 and 1 bits, got {text!r}")
    return [int(ch) for ch in cleaned]


def bits_to_str(bits: list[int]) -> str:
    return "".join(str(b) for b in bits)


def encode_block(code: Code, data: list[int]) -> list[int]:
    if len(data) != code.k:
        raise HammingError(f"a {code.name} block takes {code.k} data bits, got {len(data)}")
    word = [0] * code.n
    for pos, bit in zip(code.data_positions(), data):
        word[pos - 1] = bit
    # Each parity bit makes the XOR over its covered positions even.
    for p in code.parity_positions():
        acc = 0
        for pos in range(1, code.n + 1):
            if pos & p and pos != p:
                acc ^= word[pos - 1]
        word[p - 1] = acc
    if code.extended:
        word.append(sum(word) % 2)
    return word


def syndrome(code: Code, word: list[int]) -> int:
    """XOR of the 1-based positions holding a 1. Zero means no single-bit error."""
    s = 0
    for pos in range(1, code.n + 1):
        if word[pos - 1]:
            s ^= pos
    return s


def decode_block(code: Code, word: list[int]) -> Decoded:
    if len(word) != code.length:
        raise HammingError(f"a {code.name} block is {code.length} bits, got {len(word)}")
    bits = list(word)
    s = syndrome(code, bits)
    if not code.extended:
        if s == 0:
            return Decoded("ok", _data_of(code, bits), None)
        bits[s - 1] ^= 1
        return Decoded("corrected", _data_of(code, bits), s)
    overall_odd = sum(bits) % 2 == 1
    if s == 0 and not overall_odd:
        return Decoded("ok", _data_of(code, bits), None)
    if overall_odd:
        # An odd number of flips: assume one. Syndrome 0 means the overall bit itself.
        pos = code.length if s == 0 else s
        bits[pos - 1] ^= 1
        return Decoded("corrected", _data_of(code, bits), pos)
    # Syndrome non-zero with even parity: two errors, which we cannot fix.
    return Decoded("detected", None, None)


def _data_of(code: Code, bits: list[int]) -> list[int]:
    return [bits[p - 1] for p in code.data_positions()]


def encode_bits(code: Code, bits: list[int], pad: bool = False) -> list[int]:
    if pad:
        bits = bits + [0] * ((-len(bits)) % code.k)
    if not bits or len(bits) % code.k:
        raise HammingError(
            f"data must be a multiple of {code.k} bits for {code.name}, got {len(bits)} "
            "(use --pad to zero-fill)"
        )
    out: list[int] = []
    for i in range(0, len(bits), code.k):
        out += encode_block(code, bits[i : i + code.k])
    return out


def split_blocks(code: Code, bits: list[int]) -> list[list[int]]:
    if not bits or len(bits) % code.length:
        raise HammingError(
            f"codeword stream must be a multiple of {code.length} bits for {code.name}, "
            f"got {len(bits)}"
        )
    return [bits[i : i + code.length] for i in range(0, len(bits), code.length)]


def flip_positions(bits: list[int], positions: list[int]) -> list[int]:
    out = list(bits)
    for pos in positions:
        if not 1 <= pos <= len(out):
            raise HammingError(f"position {pos} is outside 1..{len(out)}")
        out[pos - 1] ^= 1
    return out


def coverage(code: Code) -> list[tuple[int, list[int]]]:
    """Which transmitted positions each parity bit checks."""
    return [
        (p, [pos for pos in range(1, code.n + 1) if pos & p])
        for p in code.parity_positions()
    ]


def simulate(code: Code, p: float, trials: int, seed: int) -> dict:
    """Monte Carlo: each transmitted bit flips independently with probability p."""
    if not 0 <= p <= 1:
        raise HammingError("p must be between 0 and 1")
    if trials < 1:
        raise HammingError("trials must be at least 1")
    rng = random.Random(seed)
    tally = {"clean": 0, "single": 0, "multi_recovered": 0, "detected": 0, "miscorrected": 0}
    for _ in range(trials):
        data = [rng.randint(0, 1) for _ in range(code.k)]
        word = encode_block(code, data)
        flips = [i for i in range(code.length) if rng.random() < p]
        for i in flips:
            word[i] ^= 1
        result = decode_block(code, word)
        if not flips:
            tally["clean"] += 1
        elif result.status == "detected":
            tally["detected"] += 1
        elif result.data == data:
            tally["single" if len(flips) == 1 else "multi_recovered"] += 1
        else:
            tally["miscorrected"] += 1
    n = code.length
    expected = {
        "clean": (1 - p) ** n,
        "single": n * p * (1 - p) ** (n - 1),
    }
    expected["multi"] = 1 - expected["clean"] - expected["single"]
    return {"tally": tally, "expected": expected, "trials": trials}


def verify(code: Code) -> tuple[list[str], bool]:
    """Exhaustively check minimum distance and single-error correction."""
    if code.k > VERIFY_MAX_K:
        raise HammingError(f"verify brute-forces 2^k words; use a code with k <= {VERIFY_MAX_K}")
    lines = [
        f"code {code.name}: r={code.r} n={code.length} k={code.k} "
        f"extended={'yes' if code.extended else 'no'}"
    ]
    min_weight = None
    bad_corrections = 0
    for data in product([0, 1], repeat=code.k):
        word = encode_block(code, list(data))
        if any(data):
            weight = sum(word)
            min_weight = weight if min_weight is None else min(min_weight, weight)
        for pos in range(code.length):
            broken = list(word)
            broken[pos] ^= 1
            result = decode_block(code, broken)
            if result.data != list(data):
                bad_corrections += 1
    d_ok = min_weight == code.min_distance
    c_ok = bad_corrections == 0
    words = (2 ** code.k) * code.length
    lines.append(
        f"minimum distance {min_weight} (expected {code.min_distance}): "
        + ("PASS" if d_ok else "FAIL")
    )
    lines.append(
        f"single-bit correction over {words} error cases: "
        + ("PASS" if c_ok else f"FAIL ({bad_corrections} bad)")
    )
    return lines, d_ok and c_ok


def run_batch(path: str) -> tuple[int, list[str]]:
    """Run a batch file. Returns (checks run, problems as 'line N: message')."""
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise HammingError(f"cannot read {path}: {exc.strerror}")
    code: Code | None = None
    checks = 0
    problems: list[str] = []
    for lineno, raw in enumerate(lines, start=1):
        text = raw.strip()
        if not text or text.startswith("#"):
            continue
        tokens = text.split()
        try:
            if tokens[0] == "code" and len(tokens) == 2:
                code = parse_code(tokens[1])
                continue
            if code is None:
                raise HammingError("no 'code N,K' line before this check")
            checks += 1
            _check_line(code, tokens)
        except HammingError as exc:
            problems.append(f"line {lineno}: {exc}")
    if checks == 0 and not problems:
        problems.append("no checks found in file")
    return checks, problems


def _check_line(code: Code, tokens: list[str]) -> None:
    """Handle 'encode BITS [expect CODEWORD]' and 'decode WORD [expect STATUS [DATA]]'."""
    cmd = tokens[0]
    if cmd == "encode" and len(tokens) in (2, 4):
        encoded = encode_bits(code, parse_bits(tokens[1]))
        if len(tokens) == 4:
            _expect_word(tokens, "expect")
            want = parse_bits(tokens[3])
            if encoded != want:
                raise HammingError(
                    f"encode {tokens[1]} gave {bits_to_str(encoded)}, expected {tokens[3]}"
                )
        return
    if cmd == "decode" and len(tokens) in (2, 4, 5):
        word = parse_bits(tokens[1])
        if len(tokens) == 2:
            _check_block(code, word)
            return
        if tokens[2] != "expect":
            raise HammingError(f"expected 'expect' after the word, got {tokens[2]!r}")
        result = _check_block(code, word)
        if tokens[3] != result.status:
            raise HammingError(f"decode {tokens[1]} gave {result.status}, expected {tokens[3]}")
        if len(tokens) == 5:
            want = parse_bits(tokens[4])
            if result.data != want:
                got = "none" if result.data is None else bits_to_str(result.data)
                raise HammingError(f"decode {tokens[1]} gave data {got}, expected {tokens[4]}")
        return
    raise HammingError(f"unknown or malformed directive {cmd!r}")


def _expect_word(tokens: list[str], word: str) -> None:
    if tokens[2] != word:
        raise HammingError(f"expected {word!r} after the bits, got {tokens[2]!r}")


def _check_block(code: Code, word: list[int]) -> Decoded:
    if len(word) != code.length:
        raise HammingError(f"a {code.name} block is {code.length} bits, got {len(word)}")
    return decode_block(code, word)


# ---------------------------------------------------------------- command line

def _cmd_info(args) -> int:
    code = parse_code(args.code)
    kind = "Extended Hamming (SECDED)" if code.extended else "Hamming"
    rate = code.k / code.length
    overhead = (code.length - code.k) / code.k * 100
    parity = f"{code.r}" + (" + 1 overall" if code.extended else "")
    behaviour = (
        "corrects 1 error, detects 2" if code.extended
        else "corrects 1 error (2 errors may be miscorrected)"
    )
    print(f"{kind} {code.name}")
    print(f"  {'parity bits':<20}{parity}")
    print(f"  {'data bits per block':<20}{code.k}")
    print(f"  {'block length':<20}{code.length}")
    print(f"  {'code rate':<20}{rate:.4f}")
    print(f"  {'overhead':<20}{overhead:.1f}%")
    print(f"  {'minimum distance':<20}{code.min_distance}")
    print(f"  {'behaviour':<20}{behaviour}")
    print(f"  {'data positions':<20}{','.join(str(x) for x in code.data_positions())}")
    print("  parity coverage (1-based transmitted positions):")
    for p, covered in coverage(code):
        print(f"    p{p:<3} covers {','.join(str(x) for x in covered)}")
    if code.extended:
        print(f"    overall covers 1..{code.n} and itself")
    return 0


def _cmd_encode(args) -> int:
    code = parse_code(args.code)
    print(bits_to_str(encode_bits(code, parse_bits(args.bits), pad=args.pad)))
    return 0


def _cmd_decode(args) -> int:
    code = parse_code(args.code)
    blocks = split_blocks(code, parse_bits(args.bits))
    worst = 0
    for index, block in enumerate(blocks, start=1):
        result = decode_block(code, block)
        if result.status == "detected":
            worst = 1
            print(f"block {index}: detected   double error, no data recovered")
        elif result.status == "corrected":
            print(
                f"block {index}: corrected  data {bits_to_str(result.data)} "
                f"(flipped position {result.error_position})"
            )
        else:
            print(f"block {index}: ok         data {bits_to_str(result.data)}")
    return worst


def _cmd_corrupt(args) -> int:
    positions = []
    for piece in args.positions.split(","):
        try:
            positions.append(int(piece))
        except ValueError:
            raise HammingError(f"positions must be integers, got {piece!r}")
    print(bits_to_str(flip_positions(parse_bits(args.bits), positions)))
    return 0


def _cmd_simulate(args) -> int:
    code = parse_code(args.code)
    out = simulate(code, args.p, args.trials, args.seed)
    tally, expected, trials = out["tally"], out["expected"], out["trials"]
    multi_seen = tally["multi_recovered"] + tally["detected"] + tally["miscorrected"]
    print(f"code {code.name}, p={args.p}, {trials} trials, seed {args.seed}")
    print(f"{'outcome':<26}{'observed':>12}{'rate':>10}{'expected':>12}")
    rows = [
        ("clean (no flips)", tally["clean"], expected["clean"]),
        ("one flip, corrected", tally["single"], expected["single"]),
        ("two or more flips", multi_seen, expected["multi"]),
    ]
    for label, seen, exp in rows:
        print(f"{label:<26}{seen:>12}{seen / trials:>10.4f}{exp:>12.4f}")
    print("two or more flips split:")
    print(f"  recovered by chance {tally['multi_recovered']:>8}")
    print(f"  detected (no data)  {tally['detected']:>8}")
    print(f"  silently wrong      {tally['miscorrected']:>8}")
    return 0


def _cmd_verify(args) -> int:
    code = parse_code(args.code)
    lines, ok = verify(code)
    for line in lines:
        print(line)
    print("VERIFIED" if ok else "FAILED")
    return 0 if ok else 1


def _cmd_batch(args) -> int:
    checks, problems = run_batch(args.file)
    for problem in problems:
        print(problem)
    print(f"{checks} checks, {len(problems)} problems")
    return 1 if problems else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hammingkit",
        description="Offline Hamming and extended SECDED error-correcting codes.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_code(p, default="7,4"):
        p.add_argument("-c", "--code", default=default, help="N,K such as 7,4 or 8,4 (default 7,4)")

    p = sub.add_parser("info", help="describe a code and its parity coverage")
    p.add_argument("code", nargs="?", default="7,4")
    p.set_defaults(func=_cmd_info)

    p = sub.add_parser("encode", help="encode data bits into codewords")
    p.add_argument("bits")
    add_code(p)
    p.add_argument("--pad", action="store_true", help="zero-fill the last block")
    p.set_defaults(func=_cmd_encode)

    p = sub.add_parser("decode", help="decode a codeword stream, correcting single errors")
    p.add_argument("bits")
    add_code(p)
    p.set_defaults(func=_cmd_decode)

    p = sub.add_parser("corrupt", help="flip bits at 1-based positions, for testing")
    p.add_argument("bits")
    p.add_argument("positions", help="comma-separated positions, e.g. 3,6")
    p.set_defaults(func=_cmd_corrupt)

    p = sub.add_parser("simulate", help="Monte Carlo error rates against exact binomial odds")
    add_code(p)
    p.add_argument("--p", type=float, required=True, help="per-bit flip probability")
    p.add_argument("--trials", type=int, default=10000)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=_cmd_simulate)

    p = sub.add_parser("verify", help="exhaustively prove distance and correction")
    add_code(p)
    p.set_defaults(func=_cmd_verify)

    p = sub.add_parser("batch", help="run a checks file; exit 1 on any failure")
    p.add_argument("file")
    p.set_defaults(func=_cmd_batch)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except HammingError as exc:
        print(f"hammingkit: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
