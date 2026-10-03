import argparse
from common import ROOT, download


def main():
    p = argparse.ArgumentParser(description="Inspect the licensed upstream sampler API without vendoring it")
    p.parse_args()
    for f in ["data.py", "models.py"]:
        path = download(f"https://raw.githubusercontent.com/hmdlab/raptgen/c4986ca9fa439b9389916c05829da4ff9c30d6f3/raptgen/{f}", ROOT / f"data/external/raptgen/raptgen/{f}")
        lines = path.read_text(encoding="utf-8").splitlines()
        wanted = [i for i,s in enumerate(lines) if any(t in s for t in ["class ProfileHMMSampler", "def sample", "def calc_seq_proba", "class CNN_PHMM_VAE", "def get_dataloader", "def encoder", "def forward", "random_region_length =", "best_loss"]) ]
        for i in wanted:
            print(f, i+1, "\n".join(lines[i:i+22]))


if __name__ == "__main__":
    main()
