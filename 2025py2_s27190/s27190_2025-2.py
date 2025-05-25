from Bio import Entrez, SeqIO
import pandas as pd, matplotlib.pyplot as plt

def setup_ncbi(email, key):
    Entrez.email, Entrez.api_key = email, key

def search(taxid, min_len, max_len):
    term = f"txid{taxid}[Organism] AND {min_len}:{max_len}[SLEN]"
    r = Entrez.read(Entrez.esearch(db="nucleotide", term=term, usehistory="y"))
    return r["WebEnv"], r["QueryKey"], int(r["Count"])

def fetch(webenv, query_key, count):
    h = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                      retmax=min(count, 100), webenv=webenv, query_key=query_key)
    return list(SeqIO.parse(h, "gb"))

def to_csv(records, path):
    df = pd.DataFrame([{"Accession": r.id, "Length": len(r.seq), "Description": r.description} for r in records])
    df.to_csv(path, index=False)
    return df

def plot(df, path):
    df.sort_values("Length", ascending=False, inplace=True)
    plt.figure(figsize=(10,4))
    plt.plot(df["Accession"], df["Length"], 'o-')
    plt.xticks(rotation=90, fontsize=6)
    plt.ylabel("Length")
    plt.tight_layout()
    plt.savefig(path)

if __name__ == "__main__":
    email = input("NCBI email:")
    key = input("NCBI API key:")
    taxid = input("TaxID:")
    min_len = input("Min length:")
    max_len = input("Max length:")

    setup_ncbi(email, key)
    try:
        w, k, n = search(taxid, min_len, max_len)
        print(f"Found {n} records.")
        if n == 0: exit()

        recs = fetch(w, k, n)
        csv = f"taxid_{taxid}_filtered.csv"
        png = f"taxid_{taxid}_plot.png"
        df = to_csv(recs, csv)
        plot(df, png)
        print(f"Saved {csv} and {png}")
    except Exception as e:
        print("Error:", e)
