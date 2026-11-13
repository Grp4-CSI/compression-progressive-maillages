"""Lecture OBJ et quantification."""
import numpy as np

from io_utils import Quantizer, read_obj


def test_quantization_error_is_at_most_half_a_step(sphere):
    v, _ = sphere
    q = Quantizer.fit(v, 12)
    err = np.abs(q.dequantize(q.quantize(v)) - v)
    assert err.max() <= q.step / 2 + 1e-12
    assert q.quantize(v).min() >= 0 and q.quantize(v).max() <= 4095


def test_quantizer_header_roundtrip(sphere):
    q = Quantizer.fit(sphere[0], 10)
    assert Quantizer.from_bytes(q.to_bytes()) == q


def test_read_obj_triangulates_and_drops_unused(tmp_path):
    path = tmp_path / "quad.obj"
    path.write_text("v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 9 9 9\nf 1/1/1 2/2/2 3/3/3 4/4/4\n")
    v, f = read_obj(str(path))
    assert len(v) == 4
    assert f.tolist() == [[0, 1, 2], [0, 2, 3]]
