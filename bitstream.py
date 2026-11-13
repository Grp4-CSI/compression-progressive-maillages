"""Flux binaire compact d'un ProgressiveMesh (fichier .pmc).

Disposition (petit-boutiste) :
    signature "PMC3" | quantificateur | num_vertices, n0, m0, nombre de niveaux
    bloc(zlib positions de M0, uint16) | bloc(zlib faces de M0, uint32)
    bloc(zlib faces mises de côté, uint32 ; vide si le maillage est une variété)
    pour chaque niveau : bloc(schéma + zlib étiquettes) | bloc(Huffman résidus)
où bloc = longueur sur 4 octets + contenu. Le nombre d'étiquettes d'un niveau
n'est pas écrit : chaque sommet réinséré ajoute 2 faces, le décodeur le recalcule.
"""
from __future__ import annotations

import heapq
import struct
import zlib
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from coloring import pack_labels, unpack_labels
from io_utils import Quantizer
from levels import Level, ProgressiveMesh

MAGIC = b"PMC3"
COUNTS = struct.Struct("<IIIH")
U32 = struct.Struct("<I")


# ---------------------------------------------------------- entiers variables
def _zigzag(n: int) -> int:
    """Entier signé -> entier positif : 0, -1, 1, -2, 2 ... -> 0, 1, 2, 3, 4 ..."""
    return 2 * n if n >= 0 else -2 * n - 1


def _unzigzag(z: int) -> int:
    return z // 2 if z % 2 == 0 else -(z + 1) // 2


def _put_varint(out: bytearray, n: int) -> None:
    """Écrit n sur 7 bits utiles par octet ; le bit de poids fort signale une suite."""
    while True:
        byte = n & 0x7F
        n >>= 7
        out.append(byte | (0x80 if n else 0))
        if not n:
            return


def _get_varint(data: bytes, pos: int) -> tuple[int, int]:
    n = shift = 0
    while True:
        byte = data[pos]
        pos += 1
        n |= (byte & 0x7F) << shift
        shift += 7
        if not byte & 0x80:
            return n, pos


# ------------------------------------------------------------------ Huffman
def _code_lengths(freq: Counter) -> dict[int, int]:
    """Longueurs des codes de Huffman : on fusionne à chaque tour les deux groupes
    les moins fréquents, dont chaque symbole descend d'un niveau dans l'arbre."""
    if len(freq) == 1:
        return {next(iter(freq)): 1}
    heap = [(f, i, [s]) for i, (s, f) in enumerate(sorted(freq.items()))]
    heapq.heapify(heap)
    lengths = dict.fromkeys(freq, 0)
    tie = len(heap)
    while len(heap) > 1:
        f1, _, s1 = heapq.heappop(heap)
        f2, _, s2 = heapq.heappop(heap)
        for s in s1 + s2:
            lengths[s] += 1
        heapq.heappush(heap, (f1 + f2, tie, s1 + s2))
        tie += 1
    return lengths


def _canonical_codes(lengths: dict[int, int]) -> dict[int, tuple[int, int]]:
    """Symbole -> (code, longueur), attribution canonique : seules les longueurs
    ont besoin d'être transmises."""
    codes, code, prev = {}, 0, 0
    for sym, length in sorted(lengths.items(), key=lambda t: (t[1], t[0])):
        code <<= length - prev
        codes[sym] = (code, length)
        code += 1
        prev = length
    return codes


def huffman_encode(values: np.ndarray) -> bytes:
    """Nombre de valeurs, table (symbole en zigzag, longueur), puis les codes collés."""
    values = [int(x) for x in np.asarray(values).reshape(-1)]
    out = bytearray()
    _put_varint(out, len(values))
    if not values:
        return bytes(out)
    lengths = _code_lengths(Counter(values))
    _put_varint(out, len(lengths))
    for sym, length in sorted(lengths.items(), key=lambda t: (t[1], t[0])):
        _put_varint(out, _zigzag(sym))
        out.append(length)
    codes = _canonical_codes(lengths)
    bits = "".join(format(codes[v][0], f"0{codes[v][1]}b") for v in values)
    bits += "0" * (-len(bits) % 8)
    out += int(bits, 2).to_bytes(len(bits) // 8, "big") if bits else b""
    return bytes(out)


def huffman_decode(data: bytes) -> np.ndarray:
    """Relit la table, puis les bits un par un jusqu'au nombre de valeurs annoncé."""
    count, pos = _get_varint(data, 0)
    if count == 0:
        return np.zeros(0, dtype=np.int64)
    n_sym, pos = _get_varint(data, pos)
    lengths = {}
    for _ in range(n_sym):
        z, pos = _get_varint(data, pos)
        lengths[_unzigzag(z)] = data[pos]
        pos += 1
    decode_table = {(length, code): sym for sym, (code, length) in _canonical_codes(lengths).items()}
    payload = data[pos:]
    bits = format(int.from_bytes(payload, "big"), f"0{8 * len(payload)}b") if payload else ""
    values, code, length, i = [], 0, 0, 0
    while len(values) < count:
        code = (code << 1) | (bits[i] == "1")
        length += 1
        i += 1
        sym = decode_table.get((length, code))
        if sym is not None:
            values.append(sym)
            code = length = 0
    return np.array(values, dtype=np.int64)


# ------------------------------------------------------------------- flux
@dataclass
class StreamStats:
    """Taille en octets de chaque partie du flux, niveau par niveau."""
    header_bytes: int = 0
    base_bytes: int = 0
    connectivity_bytes: list[int] = field(default_factory=list)
    geometry_bytes: list[int] = field(default_factory=list)

    @property
    def total_bytes(self) -> int:
        return self.header_bytes + self.base_bytes + sum(self.connectivity_bytes) + sum(self.geometry_bytes)

    def cumulative_bits(self) -> list[int]:
        """Bits reçus après le maillage de base, puis après chaque niveau."""
        acc = 8 * (self.header_bytes + self.base_bytes)
        out = [acc]
        for c, g in zip(self.connectivity_bytes, self.geometry_bytes):
            acc += 8 * (c + g)
            out.append(acc)
        return out


def _block(payload: bytes) -> bytes:
    return U32.pack(len(payload)) + payload


def encode(pm: ProgressiveMesh) -> tuple[bytes, StreamStats]:
    """Écrit le flux binaire et renvoie aussi la taille de chaque partie."""
    stats = StreamStats()
    header = MAGIC + pm.quantizer.to_bytes() + COUNTS.pack(
        pm.num_vertices, len(pm.base_positions), len(pm.base_faces), len(pm.levels))
    base = _block(zlib.compress(pm.base_positions.astype("<u2").tobytes(), 9)) \
        + _block(zlib.compress(pm.base_faces.astype("<u4").tobytes(), 9)) \
        + _block(zlib.compress(pm.fixed_faces.astype("<u4").tobytes(), 9))
    stats.header_bytes, stats.base_bytes = len(header), len(base)
    parts = [header, base]
    for level in pm.levels:
        conn = _block(bytes([level.scheme]) + zlib.compress(pack_labels(level.labels, level.scheme), 9))
        geom = _block(huffman_encode(level.residuals))
        stats.connectivity_bytes.append(len(conn))
        stats.geometry_bytes.append(len(geom))
        parts += [conn, geom]
    return b"".join(parts), stats


def decode(data: bytes) -> ProgressiveMesh:
    """Relit un flux binaire (sans vertex_order, qui n'est jamais transmis)."""
    if data[:4] != MAGIC:
        raise ValueError("ce n'est pas un flux PMC3")
    pos = 4
    quantizer = Quantizer.from_bytes(data[pos:pos + Quantizer.HEADER.size])
    pos += Quantizer.HEADER.size
    num_vertices, n0, m0, num_levels = COUNTS.unpack_from(data, pos)
    pos += COUNTS.size

    def read_block() -> bytes:
        nonlocal pos                    # curseur de lecture de la fonction englobante
        (length,) = U32.unpack_from(data, pos)
        pos += U32.size + length
        return data[pos - length:pos]

    base_positions = np.frombuffer(zlib.decompress(read_block()), "<u2").astype(np.int64).reshape(n0, 3)
    base_faces = np.frombuffer(zlib.decompress(read_block()), "<u4").astype(np.int64).reshape(m0, 3)
    fixed_faces = np.frombuffer(zlib.decompress(read_block()), "<u4").astype(np.int64).reshape(-1, 3)
    levels, num_faces = [], m0
    for _ in range(num_levels):
        conn = read_block()
        scheme = conn[0]
        labels = unpack_labels(zlib.decompress(conn[1:]), num_faces, scheme)
        residuals = huffman_decode(read_block()).reshape(-1, 3)
        levels.append(Level(scheme, labels, residuals))
        num_faces += 2 * len(residuals)
    return ProgressiveMesh(quantizer, base_positions, base_faces, levels, num_vertices, fixed_faces)
