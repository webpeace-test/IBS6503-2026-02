from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[2]
path = root / "data" / "airway_scaledcounts.subset.tsv"

# 이 파일은 소수점을 쉼표로 표기한다 (723,0 = 723.0)
counts = pd.read_csv(path, sep="\t", decimal=",", index_col="ensgene")

counts["control_mean"] = (counts["1_control"] + counts["2_control"]) / 2
counts["treated_mean"] = (counts["1_treated"] + counts["2_treated"]) / 2

# 발현량이 적은 유전자는 작은 차이로도 배수가 크게 튀므로 뺀다.
# 두 평균이 모두 0보다 커야 log2 배수가 inf/NaN이 되지 않는다.
min_mean = 10
kept = counts[(counts["control_mean"] >= min_mean) & (counts["treated_mean"] >= min_mean)].copy()

kept["log2fc"] = np.log2(kept["treated_mean"] / kept["control_mean"])

# 쌍(세포주)마다 따로 계산한 배수. 한 쌍에서만 늘어난 것인지 확인하는 용도.
kept["log2fc_pair1"] = np.log2(kept["1_treated"] / kept["1_control"])
kept["log2fc_pair2"] = np.log2(kept["2_treated"] / kept["2_control"])

up10 = kept.sort_values("log2fc", ascending=False).head(10)

print(len(counts), "->", len(kept), f"genes with both means >= {min_mean}")
print(up10[["control_mean", "treated_mean", "log2fc", "log2fc_pair1", "log2fc_pair2"]].round(2))
