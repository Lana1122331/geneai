import streamlit as st
import matplotlib.pyplot as plt

st.markdown("""
<style>
body, .stApp {
    background: linear-gradient(to bottom right, #044B7F, #00A5A0);
    color: white;
    font-family: 'Segoe UI', sans-serif;
    font-size: 30px;  /* هذا للخط العام */
}

/* تكبير حجم نص التسميات (Labels) فوق الحقول - اللون أبيض */
.stTextInput label, 
.stTextArea label {
    font-weight: bold;
    color: white !important;
    font-size: 25px !important;
}

/* textarea لـ DNA الطبيعي */
textarea.stTextArea > div > textarea {
    background-color: white !important;
    color: black !important;
    font-size: 18px !important;  /* حجم الخط داخل الصندوق */
    border-radius: 5px !important;
    height: 45px !important;
    padding: 12px !important;
    border: none !important;
    resize: vertical;
}

/* input عادي لـ DNA المطفر */
.stTextInput > div > div > input {
    background-color: white;
    border: none;
    border-radius: 5px;
    height: 45px;
    padding-left: 12px;
    font-size: 18px;  /* حجم الخط داخل الصندوق */
    color: black !important;
}

.result-box {
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 20px;
    padding: 25px;
    margin-top: 25px;
}

/* جدول النتائج: حجم خط أصغر ونص أسود */
.result-box table {
    width: 100%;
    font-size: 24px;  /* حجم الخط داخل الجدول أصغر */
    border-collapse: collapse;
    color: black !important;
    background-color: white !important;
}

.result-box table td, .result-box table th {
    border: 1px solid #ddd;
    padding: 12px;
    text-align: center;
    font-weight: normal;
    font-family: 'Segoe UI', sans-serif;
}

/* رؤوس الجدول بالعربي (أبيض) */
.result-box table th {
    background-color: #003366;
    color: white !important;
    font-weight: bold;
    font-size: 20px;  /* حجم الخط للعناوين داخل الجدول */
}
</style>
""", unsafe_allow_html=True)

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
        return "استبدال" if diff > 0 else "لا يوجد"
    elif len(normal) > len(mutated):
        return "حذف"
    elif len(normal) < len(mutated):
        return "إدخال"
    return "غير معروف"

def get_known_diagnosis(normal_codon, mutated_codon):
    return known_mutations.get((normal_codon, mutated_codon), "No known disease associated with this mutation.")

def estimate_risk(protein1, protein2):
    if not protein1 or not protein2:
        return 100
    diff = sum(1 for a, b in zip(protein1, protein2) if a != b)
    length = max(len(protein1), len(protein2))
    return int((diff / length) * 100)

def plot_mutation_positions(dna_normal, dna_mutated):
    length = min(len(dna_normal), len(dna_mutated))
    mutation_positions = []
    for i in range(0, length, 3):
        codon_normal = dna_normal[i:i+3].upper()
        codon_mutated = dna_mutated[i:i+3].upper()
        if codon_normal != codon_mutated:
            mutation_positions.append(i // 3 + 1)
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.hlines(1, 1, len(dna_normal) // 3, colors='white', linewidth=3)
    ax.scatter(mutation_positions, [1] * len(mutation_positions), color='red', s=200, label='Mutation Site')
    ax.set_ylim(0.8, 1.2)
    ax.set_yticks([])
    ax.set_xlabel('Codon Positions', fontsize=20, color='white')
    ax.set_title('Mutation Locations on DNA Sequence', fontsize=25, color='white')
    ax.legend(loc='upper right')
    plt.tight_layout()
    st.pyplot(fig)

st.title("🧬 GeneAI - تحليل الطفرات الجينية")

normal_dna = st.text_area("🔬 أدخل تسلسل DNA الطبيعي:", height=100, max_chars=1000)
mutated_dna = st.text_input("🧬 أدخل تسلسل DNA بعد الطفرة:")

if normal_dna and mutated_dna:
    valid = all(base in "ATCGatcg" for base in normal_dna + mutated_dna)
    if not valid:
        st.error("❌ تأكد أن الحروف هي فقط A, T, C, G.")
    else:
        mrna = transcribe_dna_to_mrna(mutated_dna)
        protein_normal = translate_mrna_to_protein(transcribe_dna_to_mrna(normal_dna))
        protein_mutated = translate_mrna_to_protein(mrna)
        mutation_type = detect_mutation_type(normal_dna, mutated_dna)
        diagnosis = get_known_diagnosis(normal_dna[3:6], mutated_dna[3:6])
        risk = estimate_risk(protein_normal, protein_mutated)

        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        st.markdown(f"""
        <table>
        <thead>
            <tr>
                <th>نوع الطفرة</th>
                <th>mRNA</th>
                <th>البروتين الطبيعي</th>
                <th>البروتين المطفر</th>
                <th>التشخيص المحتمل</th>
                <th>نسبة الخطورة</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>{mutation_type}</td>
                <td>{mrna}</td>
                <td>{' – '.join(protein_normal)}</td>
                <td>{' – '.join(protein_mutated)}</td>
                <td>{diagnosis}</td>
                <td>{risk}%</td>
            </tr>
        </tbody>
        </table>
        """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

        plot_mutation_positions(normal_dna, mutated_dna)

st.markdown("---")
st.caption("🔖 هذا المشروع من تنفيذ الطالبات: **رغد الضويان** و **لانا الشهراني**")
