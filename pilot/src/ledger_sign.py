#!/usr/bin/env python3
"""Per-entry signatures and external time anchors for the campaign ledger.

Signatures: each ledger entry hash is signed with an Ed25519 key in the
OpenSSH signature format (`ssh-keygen -Y sign`, namespace
"mat-theory-ledger"). Signatures live in a sidecar file
(phase_c/ledger_signatures.jsonl) so the hash chain is untouched; the
public key is committed as an OpenSSH allowed-signers file
(phase_c/ledger_allowed_signers). The private key stays outside the
repository (default ~/.config/mat_theory/ledger_ed25519).

Anchors: the ledger head hash is sent to an RFC 3161 time-stamping
authority (default freetsa.org); only the 32-byte digest leaves the
machine. The signed token (.tsr) and the TSA certificates are committed
under phase_c/anchors/, so the time a ledger state existed is attested by
a clock the campaign does not control.

Verification needs only OpenSSH >= 8.1 (ssh-keygen -Y verify) and
OpenSSL (openssl ts -verify); both ship with Linux, macOS and Windows 10+.

Usage:
  ledger_sign.py sign   [--key PATH]      # sign every unsigned entry
  ledger_sign.py anchor                   # time-stamp the current head
  ledger_sign.py verify                   # check all signatures and anchors
"""
import argparse, json, os, shutil, subprocess, sys, tempfile, time, urllib.request
from pathlib import Path

MT = Path(__file__).resolve().parents[2]
LEDGER = MT / "phase_c" / "ledger_campaign.jsonl"
SIGS = MT / "phase_c" / "ledger_signatures.jsonl"
SIGNERS = MT / "phase_c" / "ledger_allowed_signers"
ANCH = MT / "phase_c" / "anchors"
IDENTITY = "mat-theory-campaign-agent"
NAMESPACE = "mat-theory-ledger"
TSA_URL = "https://freetsa.org/tsr"
TSA_CERTS = {"cacert.pem": "https://freetsa.org/files/cacert.pem",
             "tsa.crt": "https://freetsa.org/files/tsa.crt"}


def entries():
    return [json.loads(l) for l in LEDGER.read_text().splitlines() if l.strip()]


def signatures():
    if not SIGS.exists():
        return {}
    return {r["seq"]: r for r in map(json.loads, SIGS.read_text().splitlines()) if r}


def ssh_sign(key, message: str) -> str:
    with tempfile.TemporaryDirectory() as d:
        m = Path(d) / "m"
        m.write_text(message)
        subprocess.run(["ssh-keygen", "-q", "-Y", "sign", "-f", str(key), "-n", NAMESPACE, str(m)],
                       check=True, capture_output=True)
        return (Path(d) / "m.sig").read_text()


def ssh_verify(message: str, sig: str) -> bool:
    with tempfile.TemporaryDirectory() as d:
        s = Path(d) / "s.sig"
        s.write_text(sig)
        m = Path(d) / "m"
        m.write_bytes(message.encode())
        # message passed as a file handle, not a pipe: Windows OpenSSH waits
        # on a piped stdin from Python
        with m.open("rb") as fin:
            r = subprocess.run(["ssh-keygen", "-Y", "verify", "-f", str(SIGNERS), "-I", IDENTITY,
                                "-n", NAMESPACE, "-s", str(s)], stdin=fin,
                               capture_output=True, timeout=60)
        return r.returncode == 0


def cmd_sign(key):
    key = Path(os.path.expanduser(key))
    pub = Path(str(key) + ".pub").read_text().split()
    line = f"{IDENTITY} namespaces=\"{NAMESPACE}\" {pub[0]} {pub[1]}\n"
    if not SIGNERS.exists():
        SIGNERS.write_text(line)
    elif SIGNERS.read_text() != line:
        sys.exit("allowed-signers file names a different key; refusing to sign")
    have, n = signatures(), 0
    with SIGS.open("a") as f:
        for e in entries():
            if e["seq"] in have:
                continue
            rec = {"seq": e["seq"], "entry_hash": e["entry_hash"], "signer": IDENTITY,
                   "signed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   "entry_time": e["time"], "sig": ssh_sign(key, e["entry_hash"])}
            f.write(json.dumps(rec, sort_keys=True) + "\n")
            n += 1
    print(f"signed {n} new entries; {len(have) + n} signed in total")


def cmd_anchor():
    E = entries()
    head = E[-1]
    ANCH.mkdir(parents=True, exist_ok=True)
    for name, url in TSA_CERTS.items():
        if not (ANCH / name).exists():
            urllib.request.urlretrieve(url, ANCH / name)
    stem = ANCH / f"head_seq{head['seq']:04d}"
    subprocess.run(["openssl", "ts", "-query", "-digest", head["entry_hash"], "-sha256", "-cert",
                    "-no_nonce", "-out", f"{stem}.tsq"], check=True, capture_output=True)
    req = urllib.request.Request(TSA_URL, data=Path(f"{stem}.tsq").read_bytes(),
                                 headers={"Content-Type": "application/timestamp-query"})
    Path(f"{stem}.tsr").write_bytes(urllib.request.urlopen(req, timeout=60).read())
    gen = subprocess.run(["openssl", "ts", "-reply", "-in", f"{stem}.tsr", "-text"],
                         capture_output=True, text=True).stdout
    t = next((l.split(":", 1)[1].strip() for l in gen.splitlines() if l.startswith("Time stamp:")), "")
    meta = {"seq": head["seq"], "entry_hash": head["entry_hash"], "ledger_time": head["time"],
            "tsa": TSA_URL, "tsa_time": t}
    Path(f"{stem}.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta))


def check_anchor(meta_file: Path):
    meta = json.loads(meta_file.read_text())
    stem = str(meta_file)[:-5]
    r = subprocess.run(["openssl", "ts", "-verify", "-digest", meta["entry_hash"], "-in", f"{stem}.tsr",
                        "-CAfile", str(ANCH / "cacert.pem"), "-untrusted", str(ANCH / "tsa.crt")],
                       capture_output=True, text=True)
    return meta, r.returncode == 0 and "Verification: OK" in r.stdout


def verify(report=print):
    """Returns a list of (name, ok) checks; usable from verify_record.py."""
    out = []
    E = {e["seq"]: e for e in entries()}
    S = signatures()
    if shutil.which("ssh-keygen") is None:
        return [("ssh-keygen available for signature checks", False)]
    bad = [s for s, e in E.items() if s not in S or S[s]["entry_hash"] != e["entry_hash"]
           or not ssh_verify(e["entry_hash"], S[s]["sig"])]
    out.append((f"per-entry signatures: {len(E) - len(bad)} of {len(E)} entries signed and valid "
                f"(key in ledger_allowed_signers)", not bad))
    for mf in sorted(ANCH.glob("head_seq*.json")):
        meta, ok = check_anchor(mf)
        ok = ok and meta["seq"] in E and E[meta["seq"]]["entry_hash"] == meta["entry_hash"]
        out.append((f"RFC 3161 anchor of seq {meta['seq']} head ({meta['tsa_time']})", ok))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["sign", "anchor", "verify"])
    ap.add_argument("--key", default="~/.config/mat_theory/ledger_ed25519")
    a = ap.parse_args()
    if a.cmd == "sign":
        cmd_sign(a.key)
    elif a.cmd == "anchor":
        cmd_anchor()
    else:
        res = verify()
        for n, ok in res:
            print(("PASS " if ok else "FAIL ") + n)
        sys.exit(0 if all(ok for _, ok in res) else 1)
