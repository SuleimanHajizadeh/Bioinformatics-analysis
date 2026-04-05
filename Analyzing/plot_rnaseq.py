import pandas as pd
import matplotlib.pyplot as plt

# 1. Faylı oxuyuruq (İlk sətir kommentdir deyə onu atırıq)
df = pd.read_csv("UHR_Rep3_counts.txt", sep="\t", comment="#")

# 2. Lazımsız sütunları atırıq, yalnız Gen Adı və Say saxlayırıq
# Sonuncu sütun (bam faylının adı) say sütunudur
count_col = df.columns[-1]
data = df[['Geneid', count_col]].copy()
data.columns = ['Gene', 'Counts']

# 3. 0 olanları təmizləyirik və ən böyükdən kiçiyə düzürük
data = data[data['Counts'] > 0]
top_genes = data.sort_values(by='Counts', ascending=False).head(10)

print("Ən aktiv 10 gen:")
print(top_genes)

# 4. Qrafik çəkirik
plt.figure(figsize=(10, 6))
plt.bar(top_genes['Gene'], top_genes['Counts'], color='skyblue')
plt.xlabel('Genlər')
plt.ylabel('Oxuma Sayı (Counts)')
plt.title('RNA-Seq: Ən Yüksək İfadə Olunan 10 Gen')
plt.xticks(rotation=45)
plt.tight_layout()

# 5. Göstər və Yadda saxla
plt.savefig("rnaseq_result.png")
plt.show()
