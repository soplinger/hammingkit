"""Tests for hammingkit's core logic and command line."""
import sys
from itertools import combinations, product
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import hammingkit as hk  # noqa: E402


def all_data(code):
    return [list(bits) for bits in product([0, 1], repeat=code.k)]


def test_parse_code_accepts_base_and_extended():
    assert hk.parse_code("7,4") == hk.Code(3, False)
    assert hk.parse_code("8,4") == hk.Code(3, True)
    assert hk.parse_code("15,11") == hk.Code(4, False)
    assert hk.parse_code("16,11").extended


def test_parse_code_rejects_unsupported_shapes():
    for text in ["7,3", "9,4", "seven,four", "7"]:
        with pytest.raises(hk.HammingError):
            hk.parse_code(text)


def test_known_codewords_for_7_4_and_8_4():
    assert hk.bits_to_str(hk.encode_block(hk.parse_code("7,4"), [1, 0, 1, 1])) == "0110011"
    assert hk.bits_to_str(hk.encode_block(hk.parse_code("8,4"), [1, 0, 1, 1])) == "01100110"


def test_every_encoded_word_has_zero_syndrome():
    for text in ["7,4", "8,4", "15,11"]:
        code = hk.parse_code(text)
        for data in all_data(code)[:256]:
            assert hk.syndrome(code, hk.encode_block(code, data)) == 0


def test_coverage_for_7_4_matches_textbook():
    code = hk.parse_code("7,4")
    assert hk.coverage(code) == [
        (1, [1, 3, 5, 7]),
        (2, [2, 3, 6, 7]),
        (4, [4, 5, 6, 7]),
    ]


def test_single_error_corrected_at_every_position():
    code = hk.parse_code("15,11")
    data = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0]
    word = hk.encode_block(code, data)
    for pos in range(1, code.length + 1):
        result = hk.decode_block(code, hk.flip_positions(word, [pos]))
        assert result.status == "corrected"
        assert result.error_position == pos
        assert result.data == data


def test_extended_code_detects_every_double_error():
    code = hk.parse_code("8,4")
    for data in all_data(code):
        word = hk.encode_block(code, data)
        for a, b in combinations(range(1, code.length + 1), 2):
            result = hk.decode_block(code, hk.flip_positions(word, [a, b]))
            assert result.status == "detected"
            assert result.data is None


def test_extended_overall_parity_bit_error_is_corrected_at_last_position():
    code = hk.parse_code("8,4")
    word = hk.encode_block(code, [1, 0, 1, 1])
    result = hk.decode_block(code, hk.flip_positions(word, [8]))
    assert result.status == "corrected"
    assert result.error_position == 8
    assert result.data == [1, 0, 1, 1]


def test_encode_bits_requires_whole_blocks_unless_padded():
    code = hk.parse_code("7,4")
    with pytest.raises(hk.HammingError):
        hk.encode_bits(code, [1, 0, 1])
    padded = hk.encode_bits(code, [1, 0, 1], pad=True)
    assert len(padded) == 7
    assert hk.split_blocks(code, padded) == [padded]


def test_split_blocks_rejects_partial_block():
    code = hk.parse_code("8,4")
    with pytest.raises(hk.HammingError):
        hk.split_blocks(code, [0] * 10)


def test_flip_positions_rejects_out_of_range():
    with pytest.raises(hk.HammingError):
        hk.flip_positions([0, 1, 1], [4])


def test_simulate_is_reproducible_and_tallies_every_trial():
    code = hk.parse_code("7,4")
    first = hk.simulate(code, 0.05, 2000, seed=7)
    second = hk.simulate(code, 0.05, 2000, seed=7)
    assert first == second
    assert sum(first["tally"].values()) == 2000
    assert abs(sum(first["expected"].values()) - 1) < 1e-12


def test_simulate_with_zero_probability_is_all_clean():
    out = hk.simulate(hk.parse_code("7,4"), 0.0, 100, seed=1)
    assert out["tally"]["clean"] == 100


def test_verify_passes_for_small_codes():
    for text in ["3,1", "4,1", "7,4", "8,4", "15,11"]:
        lines, ok = hk.verify(hk.parse_code(text))
        assert ok, lines


def test_verify_refuses_large_k():
    with pytest.raises(hk.HammingError):
        hk.verify(hk.parse_code("31,26"))


def test_batch_reports_line_numbers(tmp_path):
    checks = tmp_path / "checks.txt"
    checks.write_text(
        "# header\n"
        "code 8,4\n"
        "encode 1011 expect 01100110\n"
        "decode 01100111 expect detected\n",
        encoding="utf-8",
    )
    count, problems = hk.run_batch(str(checks))
    assert count == 2
    assert problems == ["line 4: decode 01100111 gave corrected, expected detected"]


def test_batch_rejects_check_before_code_line(tmp_path):
    checks = tmp_path / "checks.txt"
    checks.write_text("encode 1011\n", encoding="utf-8")
    _, problems = hk.run_batch(str(checks))
    assert problems == ["line 1: no 'code N,K' line before this check"]


def test_bundled_example_file_passes():
    count, problems = hk.run_batch(str(ROOT / "examples" / "checks.txt"))
    assert count >= 5
    assert problems == []


def test_cli_returns_two_for_bad_input(capsys):
    assert hk.main(["encode", "10x1"]) == 2
    assert "hammingkit: error" in capsys.readouterr().err


def test_cli_batch_exit_codes(tmp_path, capsys):
    good = ROOT / "examples" / "checks.txt"
    assert hk.main(["batch", str(good)]) == 0
    bad = tmp_path / "bad.txt"
    bad.write_text("code 7,4\ndecode 0110011 expect detected\n", encoding="utf-8")
    assert hk.main(["batch", str(bad)]) == 1
    assert "line 2" in capsys.readouterr().out
