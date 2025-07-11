import streamlit as st
import matplotlib.pyplot as plt

# CSS تصميم جميل
st.markdown("""
<style>
body, .stApp {
    background: linear-gradient(to bottom right, #044B7F, #00A5A0);
    color: white;
    font-family: 'Segoe UI', sans-serif;
    font-size: 30px;
}
.stTextInput > div > div > input {
    background-color: white;
    border: none;
    border-radius: 5px;
    height: 45px;
    padding-left: 12px;
    font-size: 20px;
    color: black !important;
}
.stTextInput label {
    font-weight: bold;
    color: #ffffff;
    font-size: 20px;
}
.result-box {
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 25px;
    padding: 30px;
    margin-top: 30px;
}
.result-box table {
    width: 100%;
    color: white;
    font-size: 25px;
}
</style>
""", unsafe_allow_html=True)

# جدول الكودونات
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
    return known_mutations.get((normal_codon, mutated_codon), "لا يوجد مرض معروف مرتبط بهذه الطفرة.")

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
    ax.scatter(mutation_positions, [1] * len(mutation_positions), color='red', s=200, label='موقع الطفرة')
    ax.set_ylim(0.8, 1.2)
    ax.set_yticks([])
    ax.set_xlabel('مواقع الكودونات', fontsize=20)
    ax.set_title('مواقع الطفرات في تسلسل DNA', fontsize=25)
    ax.legend(loc='upper right')
    st.pyplot(fig)

# عنوان التطبيق
st.title("🧬 GeneAI - تحليل الطفرات الجينية")

# حقول الإدخال
normal_dna = st.text_input("🔬 أدخل تسلسل DNA الطبيعي:")
mutated_dna = st.text_input("🧬 أدخل تسلسل DNA بعد الطفرة:")

# عند الإدخال
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

        # عرض النتائج
        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        st.markdown(f"""
        <table>
        <tr><td><strong>نوع الطفرة</strong></td><td>{mutation_type}</td></tr>
        <tr><td><strong>mRNA</strong></td><td>{mrna}</td></tr>
        <tr><td><strong>البروتين الطبيعي</strong></td><td>{' – '.join(protein_normal)}</td></tr>
        <tr><td><strong>البروتين المطفر</strong></td><td>{' – '.join(protein_mutated)}</td></tr>
        <tr><td><strong>التشخيص المحتمل</strong></td><td>{diagnosis}</td></tr>
        <tr><td><strong>نسبة الخطورة</strong></td><td>{risk}%</td></tr>
        </table>
        """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

        # عرض الرسم البياني
        plot_mutation_positions(normal_dna, mutated_dna)

# توقيع المشروع
st.markdown("---")
st.caption("🔖 هذا المشروع من تنفيذ الطالبات: **رغد الضويان** و **لانا الشهراني**")
