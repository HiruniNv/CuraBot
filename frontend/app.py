import streamlit as st
import requests
from fpdf import FPDF
from datetime import datetime
import re

# Page Config
st.set_page_config(page_title="CuraBot | AI Medical Research", page_icon="⚕️", layout="wide")

# Custom CSS for Professional Dark Dashboard
st.markdown("""
    <style>
    .main { background-color: #0e1117; }
    .medical-title { font-size: 3rem; font-weight: 800; color: #ffffff; margin-bottom: 0px; }
    .medical-subtitle { color: #ffffff; font-size: 1.4rem; font-weight: 600; font-style: italic; margin-bottom: 25px; }
    .stTextArea label { color: #ffffff; font-weight: 500; font-size: 1.1rem; }
    .stTextArea textarea { background-color: #1e2129; color: white; border: 1px solid #3e424b; border-radius: 10px; }
    .result-card { background-color: #161b22; padding: 30px; border-radius: 15px; border: 1px solid #30363d; margin-top: 25px; }
    .stButton>button { background-color: #ff4b4b; color: white; border-radius: 8px; font-weight: bold; width: 100%; height: 3.5em; }
    section[data-testid="stSidebar"] { background-color: #161b22; border-right: 1px solid #30363d; }
    </style>
    """, unsafe_allow_html=True)

def sanitize_text(text):
    """Cleans Markdown and AI artifacts for PDF generation."""

    text = re.sub(r"I now can give a great answer.*", "", text)
    clean = re.sub(r'[*_#]', '', text)
    clean = clean.replace('bullet --', '').replace('bullet', '').strip()
    return clean.encode('latin-1', 'ignore').decode('latin-1')

def generate_pdf(report_text):
    """Generates a professional clinical PDF report."""
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=25)
    pdf.add_page()
    
    # Header Banner
    pdf.set_fill_color(41, 128, 185) 
    pdf.rect(0, 0, 210, 45, 'F')
    pdf.set_y(12); pdf.set_font("helvetica", 'B', 24); pdf.set_text_color(255, 255, 255)
    pdf.cell(0, 10, "CURABOT", align='C', ln=True)
    pdf.set_font("helvetica", 'B', 14); pdf.cell(0, 10, "CLINICAL RESEARCH SUMMARY", align='C', ln=True)
    pdf.set_font("helvetica", 'I', 9); pdf.cell(0, 8, f"ID: CB-{datetime.now().strftime('%Y%m%d%H%M')} | {datetime.now().strftime('%B %d, %Y')}", align='C', ln=True)
    
    # Body Styling
    pdf.set_y(55); pdf.set_text_color(44, 62, 80)
    
    report_text = re.sub(r"(?i)I now can give a great answer.*", "", report_text).strip()
    
    lines = report_text.split('\n')
    for line in lines:
        line = line.strip()
        if not line or line.startswith('---'): continue
        
        upper_line = line.upper()
        
        if "REPORTED SYMPTOMS" in upper_line:
            pdf.ln(5); pdf.set_font("helvetica", 'B', 12); pdf.set_text_color(22, 160, 133)
            pdf.cell(0, 10, "REPORTED SYMPTOMS", ln=True)
            pdf.set_draw_color(22, 160, 133); pdf.line(10, pdf.get_y(), 200, pdf.get_y()); pdf.ln(3)
        elif "RESEARCH FINDINGS" in upper_line:
            pdf.ln(5); pdf.set_font("helvetica", 'B', 12); pdf.set_text_color(41, 128, 185)
            pdf.cell(0, 10, "KEY CLINICAL RESEARCH FINDINGS", ln=True)
            pdf.set_draw_color(41, 128, 185); pdf.line(10, pdf.get_y(), 200, pdf.get_y()); pdf.ln(3)
        elif "CLINICAL CONSIDERATIONS" in upper_line:
            pdf.ln(5); pdf.set_font("helvetica", 'B', 12); pdf.set_text_color(211, 84, 0)
            pdf.cell(0, 10, "CLINICAL CONSIDERATIONS", ln=True)
            pdf.set_draw_color(211, 84, 0); pdf.line(10, pdf.get_y(), 200, pdf.get_y()); pdf.ln(3)
        elif "NEXT STEPS" in upper_line:
            pdf.ln(5); pdf.set_font("helvetica", 'B', 12); pdf.set_text_color(44, 62, 80)
            pdf.cell(0, 10, "RECOMMENDED NEXT STEPS", ln=True); pdf.ln(1)
        elif "DISCLAIMER" in upper_line and "IMPORTANT" not in upper_line:
           
            continue
        else:
            pdf.set_text_color(60, 60, 60); pdf.set_font("helvetica", '', 10)
            if line.startswith(('*', '-', '•')):
                pdf.set_x(15); pdf.write(6, chr(149) + " "); pdf.multi_cell(0, 6, sanitize_text(line[1:].strip())); pdf.ln(1)
            else:
                pdf.multi_cell(0, 6, sanitize_text(line)); pdf.ln(1)

    # Disclaimer Section
    if pdf.get_y() > 250: pdf.add_page()
    pdf.set_y(-40); pdf.set_fill_color(253, 237, 236); pdf.set_draw_color(231, 76, 60)
    pdf.set_text_color(192, 57, 43); pdf.set_font("helvetica", 'B', 8)
    pdf.cell(0, 6, "IMPORTANT MEDICAL DISCLAIMER", align='C', ln=True, border='LTR', fill=True)
    pdf.set_font("helvetica", 'I', 8)
    disclaimer = ("CuraBot is an AI research tool for informational purposes. It is NOT a medical diagnosis. "
                  "Consult a licensed physician before making health decisions.")
    pdf.multi_cell(0, 5, disclaimer, border='LBR', align='C', fill=True)
    
    return bytes(pdf.output())

# Sidebar
with st.sidebar:
    st.markdown("### 🔍 About CuraBot")
    st.info("Uses AI agents to analyze symptoms and cross-reference PubMed research.")
    st.markdown("---")
    st.markdown("### ⚠️ Emergency Info")
    st.error("Call emergency services for severe chest pain or breathing issues.")

# Main UI
st.markdown("<h1 class='medical-title'>⚕️ CuraBot</h1>", unsafe_allow_html=True)
st.markdown("<p class='medical-subtitle'>Your Intelligent Health Research Partner</p>", unsafe_allow_html=True)

user_input = st.text_area("What symptoms are you experiencing?", placeholder="e.g., severe right side stomach pain...", height=150)

col1, col2, col3 = st.columns([1,1,1])
with col2:
    search_btn = st.button("Generate Medical Analysis")

if search_btn and user_input:
    with st.status("🩺 Running Clinical Agents...", expanded=True) as status:
        try:
            response = requests.get("http://127.0.0.1:8000/analyze", params={"symptoms": user_input}, timeout=180)
            if response.status_code == 200:
                st.session_state.report_content = response.json().get("analysis")
                status.update(label="Analysis Complete!", state="complete", expanded=False)
            else:
                st.error("Backend connection failed.")
        except Exception as e:
            st.error(f"Error: {e}")

if 'report_content' in st.session_state:
    st.markdown("---")
    st.markdown("### 📋 Clinical Research Summary")
    st.markdown('<div class="result-card">', unsafe_allow_html=True)
    st.markdown(st.session_state.report_content)
    st.markdown('</div>', unsafe_allow_html=True)
    
    pdf_bytes = generate_pdf(st.session_state.report_content)
   
    st.download_button("📥 Download Clinical Report", data=pdf_bytes, 
                       file_name=f"CuraBot_Research_{datetime.now().strftime('%Y%m%d')}.pdf", mime="application/pdf", use_container_width=True)