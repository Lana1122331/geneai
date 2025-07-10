import streamlit as st

# 🌈 CSS للتصميم الخارجي
st.markdown("""
<style>
body, .stApp {
    background: linear-gradient(to bottom right, #044B7F, #00A5A0);
    color: white;
    font-family: 'Segoe UI', sans-serif;
}
input {
    color: black !important;
}
.stTextInput > div > div > input {
    background-color: white;
    border: none;
    border-radius: 5px;
    height: 40px;
    padding-left: 10px;
    font-size: 16px;
}
.stTextInput label {
    font-weight: bold;
    color: #ffffff;
    font-size: 16px;
}
.result-box {
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 12px;
    padding: 20px;
    margin-top: 20px;
}
.result-box table {
    width: 100%;
    color: white;
    font-size: 15px;
}
</style>
""", unsafe_allow_html=True)

# 🧬 جدول الكودونات
codon_table = {
    'AUG': 'Methionine', 'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
    'UUA': 'Leucine', 'UUG': 'Leucine', 'CUU': 'Leucine', 'CUC': 'Leucine',
    'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
    'UAU': 'Tyrosine', 'UAC': 'Tyrosine', 'UGU': 'Cysteine', 'UGC': 'Cysteine',
    'UGG': 'Tryptophan', 'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
    'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine'
}

known_mutations = {
    ('GAG', 'GTG'): 'Sickle Cell Anemia'
}

def transcribe_dna_to_mrna(dna):
    return dna.upper().replace('T', 'U')

def translate_mrna_to_protein(mrna):
    protein = []
    for i in range(0, len(mrna), 3):
        codon = mrna[i:i+3]
        if len(codon) < 3:
            break
        amino_acid = codon_table.get(codon, 'Unknown')
        if amino_acid == 'STOP':
            break
        protein.append(amino_acid)
    return protein

def detect_mutation_type(normal, mutated):
    if len(normal) == len(mutated):
        diff = sum(1 for a, b in zip(normal, mutated) if a != b)
        return "Substitution" if diff > 0 else "None"
    elif len(normal) > len(mutated):
        return "Deletion"
    elif len(normal) < len(mutated):
        return "Insertion"
    return "Unknown"

def get_known_diagnosis(normal_codon, mutated_codon):
    return known_mutations.get((normal_codon, mutated_codon), "None")

def estimate_risk(protein1, protein2):
    if not protein1 or not protein2:
        return 100
    diff = sum(1 for a, b in zip(protein1, protein2) if a != b)
    length = max(len(protein1), len(protein2))
    return int((diff / length) * 100)

# 🧬 عنوان التطبيق
st.title("DNA Analyze")

# 🧬 حقول الإدخال
normal_dna = st.text_input("Normal DNA Sequence", "")
mutated_dna = st.text_input("Mutated DNA Sequence", "")

# ✅ عند الإدخال
if normal_dna and mutated_dna:
    mrna = transcribe_dna_to_mrna(mutated_dna)
    protein_normal = translate_mrna_to_protein(transcribe_dna_to_mrna(normal_dna))
    protein_mutated = translate_mrna_to_protein(mrna)

    mutation_type = detect_mutation_type(normal_dna, mutated_dna)
    diagnosis = get_known_diagnosis(normal_dna[3:6], mutated_dna[3:6])
    risk = estimate_risk(protein_normal, protein_mutated)

    # 💡 عرض النتائج
    st.markdown('<div class="result-box">', unsafe_allow_html=True)
    st.markdown(f"""
    <table>
    <tr><td><strong>Mutation Type</strong></td><td>{mutation_type}</td></tr>
    <tr><td><strong>mRNA</strong></td><td>{mrna}</td></tr>
    <tr><td><strong>Normal Protein</strong></td><td>{' – '.join(protein_normal)}</td></tr>
    <tr><td><strong>Mutated Protein</strong></td><td>{' – '.join(protein_mutated)}</td></tr>
    <tr><td><strong>Possible Diagnosis</strong></td><td>{diagnosis}</td></tr>
    <tr><td><strong>Mutation Risk</strong></td><td>{risk}%</td></tr>
    </table>
    """, unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)


