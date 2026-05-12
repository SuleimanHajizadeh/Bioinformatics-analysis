import pandas as pd
import urllib.request
import os

print("Starting ClinVar annotation pipeline...")

# 1. Download ClinVar variant_summary
if not os.path.exists("variant_summary.txt.gz"):
    print("Downloading ClinVar data from NIH (this may take a few minutes)...")
    url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    urllib.request.urlretrieve(url, "variant_summary.txt.gz")
else:
    print("ClinVar data already downloaded.")

# 2. Load ClinVar
print("Loading ClinVar data into memory...")
# low_memory=False to avoid mixed type warnings on large datasets
clinvar = pd.read_csv("variant_summary.txt.gz", sep='\t', compression='gzip', low_memory=False)

# Keep only those with an rsID (dbSNP)
print("Filtering ClinVar for valid RSIDs...")
clinvar_rs = clinvar[clinvar['RS# (dbSNP)'] != -1]

# We are only interested in certain clinical significances for this report
interesting_clinical = clinvar_rs[clinvar_rs['ClinicalSignificance'].str.contains('Pathogenic|drug response|risk factor|protective', na=False, case=False)]

# 3. Load MyHeritage Data
print("Loading MyHeritage Data...")
dna_df = pd.read_csv('MyHeritage_raw_dna_data.csv', comment='#', quotechar='"', names=['RSID', 'CHROMOSOME', 'POSITION', 'RESULT'])
dna_df = dna_df[dna_df['RSID'] != 'RSID']

# Clean up RSIDs for merging (remove 'rs' prefix and convert to float)
dna_df['RS# (dbSNP)'] = dna_df['RSID'].str.replace('rs', '').astype(float, errors='ignore')
dna_df['RS# (dbSNP)'] = pd.to_numeric(dna_df['RS# (dbSNP)'], errors='coerce')
dna_df = dna_df.dropna(subset=['RS# (dbSNP)'])

print("Merging user DNA data with NIH ClinVar...")
# We use inner merge to find the overlap
merged = pd.merge(dna_df, interesting_clinical, on='RS# (dbSNP)', how='inner')

# Save raw overlap to CSV for deep inspection
csv_path = "MyHeritage_ClinVar_Matches.csv"
merged.to_csv(csv_path, index=False)
print(f"Raw clinical matches saved to {csv_path}")

# 4. Generate a markdown report
md_content = "# 🏥 Clinical Genomic Annotation (ClinVar)\n\n"
md_content += "> **Mənbə:** MyHeritage Raw DNA + NIH ClinVar Database\n\n"
md_content += "Bu hesabat sənin DNT faylındakı bütün 600,000+ SNP-nin ABŞ Milli Səhiyyə İnstitutunun ClinVar tibbi məlumat bazası ilə avtomatik çarpazlaşdırılması nəticəsində tapılan məlumatları əks etdirir. Bu, tamhüquqlu bir bioinformatika 'pipeline' nümunəsidir.\n\n"
md_content += "⚠️ **DİQQƏT:** Bu avtomatlaşdırılmış bioinformatika analizidir. Bu hesabatda göstərilən məlumatlar ('Pathogenic' və s.) yalnız genin potensial əhəmiyyətini göstərir və sənin konkret xəstə olduğunu bildirmir. Hər genin dominant/resessiv tərəfləri var. **Bu tibbi diaqnoz deyil!**\n\n"

md_content += "## 💊 Dərman Reaksiyaları, Risklər və Patogenik Variantlar\n\n"
md_content += "| Gen / Xəstəlik | RSID | Sənin Genotipin | Klinik Əhəmiyyəti (ClinVar) | Phenotype (Təsvir) |\n"
md_content += "|---|---|---|---|---|\n"

# Sort by clinical significance to put Pathogenic/Drug response at top
sorted_merged = merged.sort_values('ClinicalSignificance')
seen_rsids = set()

count = 0
for index, row in sorted_merged.iterrows():
    rsid = row['RSID']
    if rsid not in seen_rsids:
        gene = str(row['GeneSymbol'])
        phenotype = str(row['PhenotypeList']).split(';')[0] # Get first phenotype to keep it clean
        significance = str(row['ClinicalSignificance'])
        genotype = row['RESULT']
        
        md_content += f"| **{gene}** | `{rsid}` | **{genotype}** | {significance} | {phenotype} |\n"
        seen_rsids.add(rsid)
        count += 1
        # Limit to 50 for the markdown report so it's readable
        if count >= 50: 
            break

md_content += f"\n\n*Ümumilikdə **{len(merged)}** ədəd genetik marker ClinVar bazası ilə eşləşdi. Tam cədvəl `{csv_path}` faylında saxlanıldı.*\n"

md_file = "MyHeritage_ClinVar_Report.md"
with open(md_file, 'w', encoding='utf-8') as f:
    f.write(md_content)

print(f"Analysis complete. Found {len(merged)} total matches. Markdown report saved to {md_file}")
