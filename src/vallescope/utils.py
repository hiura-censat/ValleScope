# workflows/src/util.py
import hashlib
import gzip

EXTS = ("fasta.gz", "fa.gz", "fasta", "fa")

class NameParseError(ValueError):
    pass

def split_name(name: str):
    if not isinstance(name, str):
        raise TypeError(f"Expected str, got {type(name)}")

    n_slash = name.count("/")
    if n_slash != 1:
        raise NameParseError(
            f"Invalid name '{name}': expected exactly one '/', got {n_slash} "
            "(format must be 'sample/base.fa[.gz]')"
        )

    sample, fname = name.split("/")
    for ext in EXTS:
        suf = "." + ext
        if fname.endswith(suf):
            base = fname[: -len(suf)]
            if not base:
                raise NameParseError(f"Empty basename in '{name}'")
            return sample, base, ext

    raise NameParseError(
        f"Invalid file extension in '{name}': expected one of {', '.join('.'+e for e in EXTS)}"
    )

def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def public_pairs(public_list):
    """config['data_path'] を受け取り、[(sample, base), ...] を返す"""
    pairs = []
    for item in public_list:
        smp, base, _ext = split_name(item["name"])
        pairs.append((smp, base))
    return pairs

def find_meta(public, sample, base):
    candidates = [f"{sample}/{base}.{e}" for e in EXTS] + [f"{base}.{e}" for e in EXTS]
    for cand in candidates:
        for d in public:
            if d.get("name") == cand:
                return d
    raise ValueError(f"data_path.name not found for {sample}/{base} (tried: {candidates})")


def nonN_len(path):
    op = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with op(path, "rt") as f:
        for ln in f:
            if ln.startswith(">"):
                continue
            s = ln.strip()
            total += len(s) - s.count("N") - s.count("n")
    return total

def pick_blast_params(path):
    L = nonN_len(path)
    if  L < 10_000:       	# 0-10 kb
        return 31, "1e-2"
    elif L < 100_000:       	# 10–100 kb
        return 31, "1e-4"
    elif L < 1_000_000:	# 100kb–1 Mb
        return 31, "1e-6"
    elif L < 5_000_000:     	# 1 Mb–5 Mb
        return 31, "1e-8"
    else:                   	# >=5 Mb
        return 31, "1e-10"

