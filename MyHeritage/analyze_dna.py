import pandas as pd
import warnings
warnings.filterwarnings('ignore')

print("Loading MyHeritage Data...")
dna_df = pd.read_csv('MyHeritage_raw_dna_data.csv', comment='#', quotechar='"', names=['RSID', 'CHROMOSOME', 'POSITION', 'RESULT'])
dna_df = dna_df[dna_df['RSID'] != 'RSID']

traits_db = {
    'rs1815739': {'gene': 'ACTN3', 'trait': 'Muscle Type (Sprinter vs Endurance)', 'alleles': {'CC': 'Sprinter (Fast-twitch)', 'CT': 'Mixed muscle type', 'TT': 'Endurance (Slow-twitch)'}},
    'rs4988235': {'gene': 'MCM6', 'trait': 'Lactose Tolerance', 'alleles': {'AA': 'Lactose Tolerant', 'AG': 'Lactose Tolerant', 'GG': 'Lactose Intolerant'}},
    'rs762551': {'gene': 'CYP1A2', 'trait': 'Caffeine Metabolism', 'alleles': {'AA': 'Fast Metabolizer', 'AC': 'Slow Metabolizer', 'CC': 'Slow Metabolizer'}},
    'rs12913832': {'gene': 'HERC2', 'trait': 'Eye Color', 'alleles': {'AA': 'Brown/Dark Eyes', 'AG': 'Hazel/Mixed', 'GG': 'Blue/Light Eyes'}},
    'rs671': {'gene': 'ALDH2', 'trait': 'Alcohol Flush (Asian Flush)', 'alleles': {'GG': 'Normal Metabolism (No flush)', 'AG': 'Moderate Flush', 'AA': 'Severe Flush'}},
    'rs4680': {'gene': 'COMT', 'trait': 'Stress Response (Worrier vs Warrior)', 'alleles': {'AA': 'Worrier (High dopamine, anxious under stress, good memory)', 'AG': 'Mixed', 'GG': 'Warrior (Low dopamine, performs well under stress)'}},
    'rs1800497': {'gene': 'DRD2', 'trait': 'Dopamine Receptor (Reward learning)', 'alleles': {'CC': 'Normal dopamine receptors', 'CT': 'Reduced dopamine receptors (prone to novelty seeking)', 'TT': 'Significantly reduced receptors (addiction prone)'}},
    'rs1799971': {'gene': 'OPRM1', 'trait': 'Pain Sensitivity', 'alleles': {'AA': 'Normal pain sensitivity', 'AG': 'Higher pain threshold, needs more anesthesia', 'GG': 'Highest pain threshold'}},
    'rs1801260': {'gene': 'CLOCK', 'trait': 'Sleep Cycle (Morning vs Night person)', 'alleles': {'TT': 'Morning Person', 'TC': 'Mixed', 'CC': 'Night Owl'}},
    'rs4343': {'gene': 'ACE', 'trait': 'Blood Pressure / Fitness', 'alleles': {'GG': 'Power/Strength athlete advantage', 'AG': 'Mixed', 'AA': 'Endurance athlete advantage'}},
    'rs9939609': {'gene': 'FTO', 'trait': 'Obesity risk / Appetite', 'alleles': {'TT': 'Normal risk', 'TA': 'Higher appetite, 1.3x obesity risk', 'AA': 'High appetite, 1.7x obesity risk'}},
    'rs1805007': {'gene': 'MC1R', 'trait': 'Red Hair / Freckles', 'alleles': {'CC': 'Normal hair', 'CT': 'Carrier for red hair', 'TT': 'Red hair, freckles, sensitive skin'}},
    'rs1042713': {'gene': 'ADRB2', 'trait': 'Asthma / Response to exercise', 'alleles': {'AA': 'Favorable response to albuterol', 'AG': 'Mixed', 'GG': 'Poor response, better lung capacity in athletes'}},
    'rs6311': {'gene': 'HTR2A', 'trait': 'Serotonin receptor (Memory/ADHD)', 'alleles': {'TT': 'Normal memory', 'TC': 'Mixed', 'CC': 'Worse verbal memory, potential ADHD risk'}},
    'rs1799990': {'gene': 'PRNP', 'trait': 'Prion disease resistance', 'alleles': {'AA': 'Normal', 'AG': 'Highly resistant to prion diseases (Kuru, Mad Cow)', 'GG': 'Completely resistant'}},
    'rs4149056': {'gene': 'SLCO1B1', 'trait': 'Statin-induced myopathy (Muscle pain from meds)', 'alleles': {'TT': 'Normal risk', 'TC': 'Intermediate risk', 'CC': 'High risk of muscle pain'}},
    'rs333': {'gene': 'CCR5', 'trait': 'HIV Resistance (Delta32)', 'alleles': {'--': 'Resistant to HIV (Homozygous Delta32)', 'I-': 'Partial resistance', 'II': 'Normal susceptibility'}},
    'rs53576': {'gene': 'OXTR', 'trait': 'Oxytocin Receptor (Empathy/Social)', 'alleles': {'GG': 'More empathetic, seeks social support', 'AG': 'Mixed', 'AA': 'Less empathetic, independent'}},
    'rs1800955': {'gene': 'DRD4', 'trait': 'Novelty Seeking (Wanderlust)', 'alleles': {'CC': 'Normal', 'CT': 'Novelty seeking', 'TT': 'High novelty seeking / ADHD link'}},
    'rs731236': {'gene': 'VDR', 'trait': 'Vitamin D Receptor', 'alleles': {'TT': 'Normal Vitamin D processing', 'TC': 'Lower processing efficiency', 'CC': 'Needs significantly more Vitamin D'}},
}

print("Analyzing data...")
results = []
# Fast lookup
subset = dna_df[dna_df['RSID'].isin(traits_db.keys())]

for index, row in subset.iterrows():
    rsid = row['RSID']
    genotype = row['RESULT']
    trait_info = traits_db[rsid]
    
    interpretation = trait_info['alleles'].get(genotype, "Unknown Genotype")
    if interpretation == "Unknown Genotype" and len(str(genotype)) == 2:
        rev_genotype = str(genotype)[::-1]
        interpretation = trait_info['alleles'].get(rev_genotype, "Unknown Genotype")
        
    results.append({
        'Gene': trait_info['gene'],
        'Trait': trait_info['trait'],
        'RSID': rsid,
        'Your Genotype': genotype,
        'Result': interpretation
    })

md_content = "# 🧬 Personal Genomics Analysis Report\n\n"
md_content += "> **Mənbə:** MyHeritage Raw DNA Data\n> **Analiz edən:** Personal Bioinformatics Pipeline\n\n"
md_content += "Bu hesabat MyHeritage faylındakı seçilmiş əhəmiyyətli SNP-lərin (Single Nucleotide Polymorphisms) analizidir. Aşağıdakı nəticələr genetik meyllilikləri göstərir, lakin tibbi diaqnoz deyil.\n\n"

md_content += "## 📊 Xüsusiyyət və Həyat Tərzi Analizi (Traits & Lifestyle)\n\n"
md_content += "| Gen | Xüsusiyyət (Trait) | RSID | Sənin Genotipin | Nəticə (Interpretation) |\n"
md_content += "|---|---|---|---|---|\n"

for r in results:
    md_content += f"| **{r['Gene']}** | {r['Trait']} | `{r['RSID']}` | **{r['Your Genotype']}** | {r['Result']} |\n"

md_content += "\n---\n*Qeyd: Bu hesabat Cambridge Universiteti üçün Bioinformatika portfelinin bir hissəsi kimi Python Pandası vasitəsilə generasiya edilmişdir. Tammiqyaslı bütün SNP-lərin (600,000+) analizi üçün ClinVar xəstəlik bazasından istifadə gələcək planlardadır.*\n"

with open('MyHeritage_Analysis_Report.md', 'w', encoding='utf-8') as f:
    f.write(md_content)

print("Analysis complete. Saved to MyHeritage_Analysis_Report.md")
