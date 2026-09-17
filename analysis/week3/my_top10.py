from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[2]
path = root / "data" / "airway_scaledcounts.subset.tsv"

# 이 파일은 소수점을 쉼표로 표기한다 (723,0 = 723.0)
counts = pd.read_csv(path, sep="\t", decimal=",", index_col="ensgene")

# 네 샘플의 카운트를 합친 총 발현량으로 순위를 매긴다
total = counts.sum(axis=1)
top10 = total.sort_values(ascending=False).head(10)

print(counts.loc[top10.index].assign(total=top10))
