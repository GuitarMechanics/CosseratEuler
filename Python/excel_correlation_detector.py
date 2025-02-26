import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_excel('Python/curvature_reginfo_xlvers.xlsx',sheet_name = 'autodetector', engine='openpyxl')

r2_threshold = 0.70
related_columns = []

# 각 열 쌍에 대해 선형 회귀 수행 및 결정계수 계산
num_cols = df.shape[1]
columns = df.columns

for i in range(num_cols):
    for j in range(i + 1, num_cols):  # 자기 자신 제외
        print(f"Processing columns {i} and {j}...")
        x = df[columns[i]].values.reshape(-1, 1)  # 독립 변수 (열 i)
        y = df[columns[j]].values  # 종속 변수 (열 j)

        model = LinearRegression().fit(x, y)  # 선형 회귀 모델 학습
        r2 = r2_score(y, model.predict(x))  # 결정계수 계산
        print(f"Columns {i} and {j} >> r2 score: {r2}")
        if r2 >= r2_threshold:
            related_columns.append((columns[i], columns[j], r2))
            plt.figure()
            sns.scatterplot(data = df, x = columns[i], y = columns[j],
                            hue = 'Fratio', palette = 'Set2')
            plt.title(f"{columns[i], columns[j],} r2 score = {r2}")
            plt.savefig(f"Python/cor_detectorfig/linreg_{i}_{j}.png")

# 결과 데이터프레임 생성
result_df = pd.DataFrame(related_columns, columns=["Column 1", "Column 2", "R²"])
# 결과 저장
output_file = "Python/linreg_columns.xlsx"
result_df.to_excel(output_file, index=False)
print(f"결정계수가 {r2_threshold} 이상인 열들을 '{output_file}'에 저장했습니다.")