"""Append-only committed ledger for the mat_theory pilot.

Every batch design is committed (Merkle root of its rows) BEFORE execution;
results, claims, and shard operations are appended with hash chaining so the
trail is tamper-evident. This is the minimal integrity layer from
02_zk_verifiable_computation.md: commitments + set-membership proofs, no ZK.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path


def _h(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def leaf_hash(obj) -> str:
    return _h(b"leaf:" + json.dumps(obj, sort_keys=True).encode())


def merkle_root(leaves: list[str]) -> str:
    """Merkle root over already-hashed leaves (duplicate last on odd levels)."""
    if not leaves:
        return _h(b"empty")
    level = list(leaves)
    while len(level) > 1:
        if len(level) % 2:
            level.append(level[-1])
        level = [_h(b"node:" + bytes.fromhex(a) + bytes.fromhex(b))
                 for a, b in zip(level[::2], level[1::2])]
    return level[0]


def merkle_proof(leaves: list[str], index: int) -> list[tuple[str, str]]:
    """Audit path for leaves[index]: list of (side, sibling_hash)."""
    path, level, i = [], list(leaves), index
    while len(level) > 1:
        if len(level) % 2:
            level.append(level[-1])
        sib = i + 1 if i % 2 == 0 else i - 1
        path.append(("R" if i % 2 == 0 else "L", level[sib]))
        level = [_h(b"node:" + bytes.fromhex(a) + bytes.fromhex(b))
                 for a, b in zip(level[::2], level[1::2])]
        i //= 2
    return path


def verify_membership(leaf: str, path: list[tuple[str, str]], root: str) -> bool:
    h = leaf
    for side, sib in path:
        pair = (bytes.fromhex(h) + bytes.fromhex(sib)) if side == "R" \
            else (bytes.fromhex(sib) + bytes.fromhex(h))
        h = _h(b"node:" + pair)
    return h == root


class Ledger:
    """Hash-chained append-only event log written to a JSONL file."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.entries: list[dict] = []
        if self.path.exists():
            with self.path.open() as f:
                self.entries = [json.loads(line) for line in f]

    @property
    def head(self) -> str:
        return self.entries[-1]["entry_hash"] if self.entries else _h(b"genesis")

    def append(self, kind: str, payload: dict) -> dict:
        entry = {
            "seq": len(self.entries),
            "time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "kind": kind,
            "payload": payload,
            "prev": self.head,
        }
        entry["entry_hash"] = _h(json.dumps(entry, sort_keys=True).encode())
        self.entries.append(entry)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a") as f:
            f.write(json.dumps(entry, sort_keys=True) + "\n")
        self._sign_if_campaign()
        return entry

    def _sign_if_campaign(self):
        """Sign new entries of the campaign ledger when the campaign key is
        present on this machine (see ledger_sign.py); other ledgers and
        machines without the key are left unsigned."""
        import os
        try:
            import ledger_sign
        except ImportError:
            return
        key = os.path.expanduser("~/.config/mat_theory/ledger_ed25519")
        if self.path.resolve() == ledger_sign.LEDGER.resolve() and os.path.exists(key):
            ledger_sign.cmd_sign(key)

    def verify_chain(self) -> bool:
        prev = _h(b"genesis")
        for e in self.entries:
            if e["prev"] != prev:
                return False
            body = {k: v for k, v in e.items() if k != "entry_hash"}
            if _h(json.dumps(body, sort_keys=True).encode()) != e["entry_hash"]:
                return False
            prev = e["entry_hash"]
        return True

    # Convenience: commit a batch design before execution -------------------
    def commit_batch(self, batch_id: str, rows: list[dict]) -> dict:
        leaves = [leaf_hash(r) for r in rows]
        root = merkle_root(leaves)
        return self.append("batch_commitment", {
            "batch_id": batch_id, "n_rows": len(rows), "merkle_root": root,
        })
