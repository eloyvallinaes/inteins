import re
import pandas as pd


def load_mapper():
    df = pd.read_csv("taxid2acc.tsv", sep = "\t")
    return {acc: f"{taxid}"
        for taxid, acc in df.values
    }
    return mapper


def main():
    mapper = load_mapper()
    regexp = r"(WP_[0-9]{9})\.[0-9]{1}\|intein"
    with open("rnrinteins.tree", "r") as tree:
        with open("rnrinteins_labelled.tree", "w") as new:
            for line in tree.readlines():
                if line.startswith("WP"):
                    acc = re.match(regexp, line).group(1)
                    try:
                        label = mapper[acc]
                    except KeyError:
                        label = "Unknown"
                    new.write(re.sub(regexp, label, line))
                else:
                    new.write(line)


if __name__ == "__main__":
    mapper = load_mapper()
    main()
