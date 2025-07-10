import streamlit as st
import matplotlib.pyplot as plt

# 🌈 CSS بتعديل الخطوط حسب طلبك
st.markdown("""
<style>
/* تكبير حجم نص التسميات فوق الحقول */
.stTextInput label {
    font-size: 28px !important;
    font-weight: bold;
    color: white !important;
}

/* حجم الخط داخل الصناديق البيضاء (الكتابة بداخلها) */
.stTextInput > div > div > input {
    font-size: 18px !important;
    color: black !important;
    background-color: white !important;
    border-radius: 5px !important;
    padding-left: 10px !important;
    height: 40px !important;
}

/* لو تستخدم TextArea بدل input */
.stTextArea > div > textarea {
    font-size: 18px !important;
    color: black !important;
    background-color: white !important;
    border-radius: 5px !important;
    padding: 10px !important;
}

/* صندوق النتائج */
.result-box {
    background-color: rgba(255, 255, 255, 0.1);
    border-radius: 15px;
    padding: 20px;
    margin-top: 20px;
}

/* حجم الخط داخل جدول النتائج */
.result-box table {
    width: 100%;
    font-size: 18px;
    color: white;
}

/* الخلفية العامة */
body, .stApp {
    background: linear-gradient(to bottom right, #044B7F, #00A5A0);
    color: white;
    font-family: 'Segoe UI', sans-serif;
}
</style>
""", unsafe_allow_html=True)

# جدول الكودونات mRNA شامل
codon_table = {
    'AUG': 'Methionine',
    'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
    'UUA': 'Leucine', 'UUG': 'Leucine', 'CUU': 'Leucine', 'CUC': 'Leucine',
    'CUA': 'Leucine', 'CUG': 'Leucine',
    'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
    'AGU': 'Serine', 'AGC': 'Serine',
    'UAU': 'Tyrosine', 'UAC': 'Tyrosine',
    'UGU': 'Cysteine', 'UGC': 'Cysteine',
    'UGG': 'Tryptophan',
    'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
    'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
    'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'AAU': 'Asparagine', 'AAC': 'Asparagine',
    'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine',
    'CAU': 'Histidine', 'CAC': 'Histidine',
    'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine',
    'AAA': 'Lysine', 'AAG': 'Lysine',
    'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
    'AGA': 'Arginine', 'AGG': 'Arginine',
    'GAC': 'Aspartic Acid', 'GAU': 'Aspartic Acid',
    'ACC': 'Threonine', 'ACU': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
    'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
}

known_mutations = {
    ('GAG', 'GTG'): 'أنيميا الخلايا المنجلية (Sickle Cell Anemia)',
}

def transcribe_dna_to_mrna(dna: str) -> str:
    return dna.upper().replace('T', 'U')

def translate_mrna_to_protein(mrna: str) -> list:
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

def detect_mutation_type(normal: str, mutated: str) -> str:
    if len(normal) == len(mutated):
        diff = sum(1 for a, b in zip(normal, mutated) if a != b)
        return "استبدال" if diff > 0 else "لا يوجد"
    elif len(normal) > len(mutated):
        return "حذف"
    elif len(normal) < len(mutated):
        return "إدخال"
    return "غير معروف"

def get_known_diagnosis(normal_codon: str, mutated_codon: str) -> str:
    return known_mutations.get((normal_codon, mutated_codon), "لا يوجد مرض معروف مرتبط بهذه الطفرة.")

def estimate_risk(protein1: list, protein2: list) -> int:
    if not protein1 or not protein2:
        return 100
    diff = sum(1 for a, b in zip(protein1, protein2) if a != b)
    length = max(len(protein1), len(protein2))
    return int((diff / length) * 100)

def contains_stop_codon(mrna: str) -> bool:
    stop_codons = ['UAA', 'UAG', 'UGA']
    for i in range(0, len(mrna), 3):
        codon = mrna[i:i+3]
        if codon in stop_codons:
            return True
    return False

def plot_mutation_positions(dna_normal: str, dna_mutated: str):
    length = min(len(dna_normal), len(dna_mutated))
    codon_positions = []
    mutation_positions = []
    for i in range(0, length, 3):
        codon_normal = dna_normal[i:i+3].upper()
        codon_mutated = dna_mutated[i:i+3].upper()
        codon_positions.append(i // 3 + 1)
        if codon_normal != codon_mutated:
            mutation_positions.append(i // 3 + 1)
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.hlines(1, 1, len(dna_normal) // 3, colors='white', linewidth=3)
    ax.scatter(mutation_positions, [1] * len(mutation_positions), color='red', s=200, label='موقع الطفرة')
    ax.set_ylim(0.8, 1.2)
    ax.set_yticks([])
    ax.set_xlabel('مواقع الكودونات', fontsize=20, color='white')
    ax.set_title('مواقع الطفرات في تسلسل DNA', fontsize=25, color='white')
    ax.legend(loc='upper right', fontsize=15)
    st.pyplot(fig)

st.set_page_config(page_title="GeneAI - تحليل الطفرات مع رسم", layout="centered")
st.title("🧬 GeneAI - تحليل متقدّم لتسلسل DNA مع تمثيل رسومي")
st.markdown("أدخل تسلسل DNA الطبيعي والمطفر، وسنقوم بتحليل الطفرة مع عرض تمثيل رسومي لمواقع الطفرات.")

col1, col2 = st.columns(2)
with col1:
    dna_normal = st.text_input("🔬 أدخل تسلسل DNA الطبيعي:")
with col2:
    dna_mutated = st.text_input("🧬 أدخل تسلسل DNA بعد الطفرة:")

if dna_normal and dna_mutated:
    valid = all(base in "ATCGatcg" for base in dna_normal + dna_mutated)
    if not valid:
        st.error("❌ تأكد أن الحروف هي فقط A, T, C, G.")
    else:
        mrna_normal = transcribe_dna_to_mrna(dna_normal)
        mrna_mutated = transcribe_dna_to_mrna(dna_mutated)

        protein_normal = translate_mrna_to_protein(mrna_normal)
        protein_mutated = translate_mrna_to_protein(mrna_mutated)

        mutation_type = detect_mutation_type(dna_normal, dna_mutated)

        diagnosis = "غير معروف"
        for i in range(0, min(len(dna_normal), len(dna_mutated)) - 2, 3):
            codon_normal = dna_normal[i:i+3]
            codon_mutated = dna_mutated[i:i+3]
            if codon_normal != codon_mutated:
                diagnosis = get_known_diagnosis(codon_normal, codon_mutated)
                break

        stop_alert = contains_stop_codon(mrna_mutated)
        if stop_alert:
            diagnosis += " - تحذير: وجود كودون STOP مبكر في التسلسل المطفر!"

        risk = estimate_risk(protein_normal, protein_mutated)

        st.subheader("🔍 النتائج:")
        st.write("**نوع الطفرة:**", mutation_type)
        st.write("**mRNA الطبيعي:**", mrna_normal)
        st.write("**mRNA المطفر:**", mrna_mutated)
        st.write("**البروتين الطبيعي:**", ", ".join(protein_normal))
        st.write("**البروتين المطفر:**", ", ".join(protein_mutated))

        if protein_normal != protein_mutated:
            st.warning("⚠️ تم اكتشاف تغير في تسلسل البروتين! قد تكون الطفرة مؤثرة.")
        else:
            st.success("✅ لا يوجد تغيير في البروتين الناتج.")

        st.markdown(f"**🧠 التشخيص المحتمل:** {diagnosis}")
        st.markdown(f"**🔬 نسبة الخطورة المقدّرة:** {risk}%")

        # عرض الرسم البياني لمواقع الطفرات
        plot_mutation_positions(dna_normal, dna_mutated)

st.markdown("---")
st.caption("🔖 هذا المشروع من تنفيذ الطالبات: **رغد الضويان** و **لانا الشهراني**")
