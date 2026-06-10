import torch
import esm
from typing import List
import random

# Standard 20 amino acids (order may need to match model vocab; adjust if your model uses a different order)
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"


def _get_aa_token_info(alphabet):
    """
    Get token indices and characters for amino-acid tokens from an ESM-2 alphabet.
    Returns (aa_token_ids, aa_chars) where both lists have the same length.
    """
    aa_token_ids: List[int] = []
    aa_chars: List[str] = []

    # ESM-2 alphabets expose all tokens via `all_toks`
    tokens = getattr(alphabet, "all_toks", [])
    for idx, tok in enumerate(tokens):
        if len(tok) == 1 and tok in STANDARD_AA:
            aa_token_ids.append(idx)
            aa_chars.append(tok)

    # Fallback: any single-letter alphabetic token
    if not aa_chars:
        for idx, tok in enumerate(tokens):
            if len(tok) == 1 and tok.isalpha():
                aa_token_ids.append(idx)
                aa_chars.append(tok)

    return aa_token_ids, aa_chars


def sample_peptides_with_esm2(
    model,
    alphabet,
    batch_converter,
    num_peptides: int,
    max_len: int = 50,
    temperature: float = 1.0,
    device: torch.device = None,
) -> List[str]:
    """
    Sample peptide sequences using ESM-2 logits.

    This uses an autoregressive-style procedure: for a growing prefix, we run the model,
    take logits at the last residue position, restrict to amino-acid tokens, and sample the next AA.
    Note: ESM-2 is trained as a masked LM, so this is an approximate generative use.
    """
    model.eval()
    device = device or next(model.parameters()).device

    aa_token_ids, aa_chars = _get_aa_token_info(alphabet)
    if not aa_chars:
        raise RuntimeError("Could not infer amino-acid tokens from ESM-2 alphabet.")

    # Try to identify stop tokens, e.g., EOS or <eos> among the model's tokens
    tokens = getattr(alphabet, "all_toks", [])
    stop_token_ids = []
    stop_token_strs = ["<eos>", "EOS", "[EOS]"]
    for stop_str in stop_token_strs:
        if stop_str in tokens:
            stop_token_ids.append(tokens.index(stop_str))
    # Also, sometimes eos_token or eos_idx attributes exist:
    if hasattr(alphabet, "eos_idx"):
        stop_token_ids.append(alphabet.eos_idx)
    # Remove duplicates & sanity check
    stop_token_ids = list(set([idx for idx in stop_token_ids if idx is not None and idx < len(tokens)]))
    # Map vocab token id -> amino-acid letter
    id_to_aa = {tid: ch for tid, ch in zip(aa_token_ids, aa_chars)}
    print(id_to_aa)
    peptides: List[str] = []

    with torch.no_grad():
        while len(peptides) < num_peptides:
            # Start with a single random amino acid
            seq = random.choice(aa_chars)
            stopped = False

            # Grow sequence until length < max_len (we stop at max_len - 1), or EOS token is generated
            while len(seq) < max_len - 1 and not stopped:
                data = [("seq", seq)]
                _, _, tokens_ = batch_converter(data)
                tokens_ = tokens_.to(device)

                out = model(tokens_, repr_layers=[], return_contacts=False)
                logits = out["logits"]  # [1, L_total, vocab_size]

                # Token layout is typically: [CLS] + residues + [EOS]
                # Residues occupy positions 1..len(seq), so last residue index = len(seq)
                last_res_pos = len(seq)
                # For next token: compute logits over AA and any known stop tokens
                extended_token_ids = aa_token_ids.copy()
                for stid in stop_token_ids:
                    if stid not in extended_token_ids:
                        extended_token_ids.append(stid)
                logits_last = logits[0, last_res_pos, extended_token_ids]

                if temperature <= 0:
                    next_idx_subset = logits_last.argmax().item()
                else:
                    probs = torch.softmax(logits_last.float() / temperature, dim=-1)
                    next_idx_subset = torch.multinomial(probs, 1).item()

                vocab_idx = extended_token_ids[next_idx_subset]
                if vocab_idx in stop_token_ids:
                    stopped = True
                    break
                next_aa = id_to_aa[vocab_idx]
                seq += next_aa

            peptides.append(seq)

    return peptides[:num_peptides]


def compute_perplexities(
    model,
    alphabet,
    batch_converter,
    peptides: List[str],
    device: torch.device,
    batch_size: int = 32,
) -> List[float]:
    """
    Compute (pseudo-)perplexity for each peptide under ESM-2.

    We approximate perplexity by using the model's per-position logits on the full sequence
    and computing cross-entropy over residue positions (excluding CLS/EOS/pad).
    Lower perplexity => sequence is more probable under the language model.
    """
    model.eval()
    pad_idx = alphabet.padding_idx
    perplexities: List[float] = []

    with torch.no_grad():
        for i in range(0, len(peptides), batch_size):
            batch = peptides[i : i + batch_size]
            data = [(str(j), s) for j, s in enumerate(batch)]
            _, _, tokens = batch_converter(data)
            tokens = tokens.to(device)

            out = model(tokens, repr_layers=[], return_contacts=False)
            logits = out["logits"]  # [B, L_total, vocab_size]

            # Residue positions are 1..L-2 (skip CLS at 0 and EOS at L-1)
            logits_res = logits[:, 1:-1, :]          # [B, L_res, V]
            target = tokens[:, 1:-1]                 # [B, L_res]

            log_probs = torch.log_softmax(logits_res, dim=-1)
            nll = -log_probs.gather(-1, target.unsqueeze(-1)).squeeze(-1)  # [B, L_res]

            mask = target != pad_idx
            token_count = mask.sum(dim=1)
            # Avoid division by zero
            token_count = torch.clamp(token_count, min=1)

            nll_sum = (nll * mask).sum(dim=1)
            avg_nll = nll_sum / token_count
            batch_ppl = torch.exp(avg_nll).tolist()
            perplexities.extend(batch_ppl)

    return perplexities


def main():
    num_peptides = 100
    max_len = 50

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Loading ESM-2 model on {device}...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device)
    batch_converter = alphabet.get_batch_converter()

    print(f"Sampling {num_peptides} peptides with ESM-2 (length < {max_len})...")
    peptides = sample_peptides_with_esm2(
        model,
        alphabet,
        batch_converter,
        num_peptides=num_peptides,
        max_len=max_len,
        temperature=5.0,
        device=device,
    )

    print("Computing perplexity for each peptide...")
    perplexities = compute_perplexities(
        model,
        alphabet,
        batch_converter,
        peptides,
        device=device,
        batch_size=32,
    )

    # Keep top K lowest-perplexity (highest-likelihood) peptides
    keep_k = num_peptides
    ranked = sorted(zip(peptides, perplexities), key=lambda x: x[1])  # lower perplexity is better
    top = ranked[:keep_k]

    import csv
    out_path = "sampled_esm2_peptides_with_perplexity.csv"
    with open(out_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sequence", "perplexity"])
        for seq, ppl in top:
            writer.writerow([seq, f"{ppl:.4f}"])

    print(f"Saved top {keep_k} lowest-perplexity peptides to {out_path}")


if __name__ == "__main__":
    main()

