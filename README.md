# hammingkit

**Offline Hamming and extended SECDED error-correcting codes, in one stdlib Python file.**

Encode bits, corrupt them on purpose, decode them back, measure how often a noisy channel
gets through, and prove the code does what it claims, all from the terminal. No servers,
no accounts, no network calls.

## Features

- **Hamming codes** from (3,1) up to (63,57), by `N,K` notation such as `7,4` or `15,11`.
- **Extended SECDED codes** such as `8,4` and `16,11`: correct one error and detect two.
- **Encoding** of any whole number of blocks, with optional `--pad` zero-fill for the tail.
- **Decoding** of codeword streams block by block, reporting `ok`, `corrected` (with the flipped
  position) or `detected` (double error, no data trusted).
- **Code anatomy** with `info`: rate, overhead, minimum distance, data positions and each
  parity bit's coverage.
- **Corruption tool** that flips chosen 1-based positions, so you can test by hand.
- **Monte Carlo simulation** with a seed, comparing observed outcome rates against exact
  binomial probabilities.
- **Exhaustive verification** that proves minimum distance and single-error correction for
  every data word of small codes.
- **Batch checks files** with `encode` and `decode` expectations, line-numbered problem
  reports, and exit codes for CI.
- **Strict input validation** with messages that say what was expected, never a traceback.
- **Pure standard library.** Tests use pytest only.

## Quick start

```bash
python3 hammingkit.py encode 1011
# 0110011

python3 hammingkit.py decode 0110011
# block 1: ok         data 1011

python3 hammingkit.py corrupt 0110011 2
# 0010011

python3 hammingkit.py decode 0010011
# block 1: corrected  data 1011 (flipped position 2)
```

Use the extended code when you also need to know about double errors:

```bash
python3 hammingkit.py info 8,4
python3 hammingkit.py encode 1011 -c 8,4          # 01100110
python3 hammingkit.py decode 00101110 -c 8,4      # block 1: detected   double error, no data recovered
```

A batch file turns the same checks into a CI gate:

```bash
python3 hammingkit.py batch examples/checks.txt
# 7 checks, 0 problems
```

A failing expectation is reported with its line number and makes the exit code 1:

```
line 4: decode 01100111 gave corrected, expected detected
```

## Commands

| Command | What it does |
| --- | --- |
| `info [N,K]` | Describe a code and show parity coverage |
| `encode BITS [-c N,K] [--pad]` | Encode data bits into codewords |
| `decode BITS [-c N,K]` | Decode and correct a codeword stream |
| `corrupt BITS POS[,POS...]` | Flip bits at 1-based positions |
| `simulate --p P [--trials T] [--seed S] [-c N,K]` | Monte Carlo error rates |
| `verify [-c N,K]` | Prove minimum distance and single-error correction (k ≤ 11) |
| `batch FILE` | Run `code`, `encode` and `decode` checks; exit 1 on any problem |

Exit codes: `0` success, `1` a decode detected an error, a verify failed or a batch had
problems, `2` bad input.

## Batch file format

```
# comments start with #
code 8,4
encode 1011 expect 01100110
decode 00100110 expect corrected 1011
decode 00101110 expect detected
```

Each `decode` may check the status (`ok`, `corrected`, `detected`) and, optionally, the data.
Each `encode` may check the codeword. Bits are written without spaces in batch files.

## Requirements

- Python 3.9 or newer
- pytest, only to run the tests

## Testing

```bash
python3 -m pip install pytest
python3 -m pytest -q
```

The tests cover code parsing, known codewords, zero syndromes for every codeword, single-error
correction at every position, double-error detection for every pair of positions in an
extended code, padding and block-size errors, reproducible simulation, exhaustive
verification, batch line-number reporting, the bundled example file and CLI exit codes.

## Project layout

```
hammingkit.py          the whole tool: codes, decoder, simulator, verifier, CLI
examples/checks.txt    a sample batch file for the (8,4) code
tests/test_hammingkit.py
README.md
```
