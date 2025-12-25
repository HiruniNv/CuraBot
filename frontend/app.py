import streamlit as st
import requests
from fpdf import FPDF
from datetime import datetime

# --- PAGE CONFIG ---
st.set_page_config(page_title="CuraBot | AI Medical Research", page_icon="⚕️", layout="centered")

# --- CUSTOM STYLING ---
st.markdown("""
    <style>
    .main { background-color: #f0f2f6; }
    .stTextArea label { font-weight: bold; color: #1e3d59; }
    .report-container {
        background-color: white;
        padding: 30px;
        border-radius: 15px;
        border-top: 8px solid #007bff;
        box-shadow: 0 10px 25px rgba(0,0,0,0.05);
        margin-top: 20px;
    }
    .medical-header {
        text-align: center;
        color: #1e3d59;
        font-family: 'Helvetica Neue', sans-serif;
    }
    </style>
    """, unsafe_allow_html=True)

def sanitize_text(text):
    return text.encode('latin-1', 'ignore').decode('latin-1')

def generate_pdf(report_text):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", 'B', 16)
    pdf.cell(190, 10, "CURABOT MEDICAL RESEARCH REPORT", align='C', new_x="LMARGIN", new_y="NEXT")
    pdf.ln(10)
    pdf.set_font("helvetica", size=11)
    pdf.multi_cell(0, 7, text=sanitize_text(report_text))
    return bytes(pdf.output())

# --- UI CONTENT ---
st.markdown("<h1 class='medical-header'>⚕️ CuraBot Assistant</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>Evidence-based research for clinical symptoms.</p>", unsafe_allow_html=True)

with st.container():
    # Fixed: Added a label for accessibility
    user_input = st.text_area("Describe your symptoms:", 
                              placeholder="e.g. Sharp headache and abdominal cramps for 2 hours...", 
                              height=120, 
                              label_visibility="visible")
    
    col1, col2, col3 = st.columns([1,2,1])
    with col2:
        search_btn = st.button("🚀 Generate Analysis", type="primary")

if search_btn:
    if not user_input:
        st.warning("Please provide symptom details first.")
    else:
        with st.status("🧬 **Agents are analyzing medical data...**", expanded=True) as status:
            try:
                response = requests.get(f"http://127.0.0.1:8000/analyze", params={"symptoms": user_input}, timeout=150)
                if response.status_code == 200:
                    st.session_state.report_content = response.json().get("analysis")
                    status.update(label="Analysis Complete!", state="complete", expanded=False)
                else:
                    st.error(f"Backend Error: {response.json().get('detail')}")
            except Exception as e:
                st.error("Could not connect to the CuraBot Backend. Please ensure 'main.py' is running.")

if 'report_content' in st.session_state:
    st.markdown("---")
    st.markdown("<div class='report-container'>", unsafe_allow_html=True)
    st.markdown("### 📝 Clinical Research Summary")
    st.markdown(st.session_state.report_content)
    
    pdf_bytes = generate_pdf(st.session_state.report_content)
    st.download_button(
        label="📥 Download PDF Report",
        data=pdf_bytes,
        file_name=f"CuraBot_Report_{datetime.now().strftime('%Y%m%d')}.pdf",
        mime="application/pdf",
        use_container_width=True
    )
    st.markdown("</div>", unsafe_allow_html=True)

st.sidebar.markdown("### ⚠️ Emergency Info")
st.sidebar.error("If you are experiencing severe chest pain or difficulty breathing, call emergency services immediately.")