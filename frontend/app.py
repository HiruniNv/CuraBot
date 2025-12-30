import streamlit as st
import streamlit.components.v1 as components
import requests
import pandas as pd
from fpdf import FPDF
from datetime import datetime
import time
import re

# --- 1. PAGE CONFIGURATION & SESSION STATE ---
st.set_page_config(
    page_title="CuraBot | AI Medical Research", 
    page_icon="⚕️", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# PAYHERE PAYMENT SUCCESS HANDLER
# ============================================================================
query_params = st.query_params
if query_params.get("payment") == "success":
    st.session_state.is_premium = True
    st.success("🎉 Payment Successful! Welcome to CuraBot Premium.")
    st.query_params.clear()
elif query_params.get("payment") == "cancel":
    st.warning("⚠️ Payment was cancelled. You can try again anytime.")
    st.query_params.clear()

# Initialize Session States
states = {
    'is_logged_in': False,
    'is_premium': False,
    'history': [],
    'loading': False,
    'report_content': None,
    'user_email': None 
}
for key, val in states.items():
    if key not in st.session_state:
        st.session_state[key] = val

# --- 2. HELPER FUNCTIONS ---
def is_valid_email(email):
    return re.match(r"[^@]+@[^@]+\.[^@]+", email)

def generate_pdf(report_text, symptoms, demographics):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(200, 10, "CuraBot Clinical Research Report", ln=True, align='C')
    pdf.set_font("Arial", 'I', 10)
    pdf.cell(200, 10, f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}", ln=True, align='C')
    pdf.ln(10)
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Case Summary:", ln=True)
    pdf.set_font("Arial", '', 11)
    pdf.multi_cell(0, 10, f"Symptoms: {symptoms}\nContext: {demographics}")
    pdf.ln(5)
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Clinical Synthesis:", ln=True)
    pdf.set_font("Arial", '', 11)
    clean_text = report_text.replace("**", "")
    pdf.multi_cell(0, 8, clean_text)
    pdf.ln(10)
    pdf.set_font("Arial", 'I', 8)
    pdf.multi_cell(0, 5, "DISCLAIMER: This report is AI-generated for educational purposes and does not constitute medical advice.")
    return pdf.output(dest='S').encode('latin-1', 'replace')

# ============================================================================
# PREMIUM WINDOW WITH PAYHERE CHECKOUT INTEGRATION (POPUP UPDATE)
# ============================================================================
def render_premium_content():
    st.markdown("### 🔐 Premium Features")
    st.caption("Unlock full clinical reports and PDF downloads.")
    
    if not st.session_state.is_logged_in:
        st.warning("⚠️ Please Sign In first to upgrade.")
        return
    
    st.markdown('<div class="premium-feature">📄 <b>Download Reports:</b> Save PDF health reports.</div>', unsafe_allow_html=True)
    st.markdown('<div class="premium-feature">📊 <b>Comparison View:</b> Advanced condition analysis.</div>', unsafe_allow_html=True)
    st.markdown('<div class="premium-feature">🗂️ <b>History Tracking:</b> Secure long-term tracking.</div>', unsafe_allow_html=True)
    st.markdown("---")
    
    plan = st.radio(
        "Choose Plan", 
        ["Monthly - Rs. 999", "Yearly - Rs. 8,999 (Save 25%)"], 
        key="p_radio"
    )
    
    if st.button("🚀 Proceed to Checkout", use_container_width=True):
        try:
            # INCREASED TIMEOUT to 30s to prevent 'Read timed out' errors
            res = requests.post(
                "http://127.0.0.1:8000/create-payment",
                json={
                    "plan": plan,
                    "user_email": st.session_state.user_email,
                    "user_name": "CuraBot User"
                },
                timeout=30
            )
            
            if res.status_code == 200:
                p = res.json()
                
                # --- POPUP WINDOW LOGIC ---
                # target="payhere_popup" matches the name in window.open
                form_html = f"""
                <html>
                <body>
                <form method="post" action="https://sandbox.payhere.lk/pay/checkout" id="payform" target="payhere_popup">
                    <input type="hidden" name="merchant_id" value="{p['merchant_id']}">
                    <input type="hidden" name="return_url" value="{p['return_url']}">
                    <input type="hidden" name="cancel_url" value="{p['cancel_url']}">
                    <input type="hidden" name="notify_url" value="{p['notify_url']}">
                    <input type="hidden" name="order_id" value="{p['order_id']}">
                    <input type="hidden" name="items" value="{p['items']}">
                    <input type="hidden" name="currency" value="LKR">
                    <input type="hidden" name="amount" value="{p['amount']}">
                    <input type="hidden" name="first_name" value="{p['first_name']}">
                    <input type="hidden" name="last_name" value="User">
                    <input type="hidden" name="email" value="{p['email']}">
                    <input type="hidden" name="hash" value="{p['hash']}">
                </form>
                <script>
                    // Calculate center position
                    var w = 500;
                    var h = 600;
                    var left = (screen.width/2)-(w/2);
                    var top = (screen.height/2)-(h/2);
                    
                    // Open centered popup
                    window.open("", "payhere_popup", "width="+w+",height="+h+",top="+top+",left="+left+",location=no,toolbar=no,menubar=no");
                    
                    // Submit form into the popup
                    document.getElementById('payform').submit();
                </script>
                </body>
                </html>
                """
                components.html(form_html, height=0)
                st.success("✅ Checkout popup opened. Complete payment to activate Premium!")
            else:
                st.error(f"❌ Failed to initialize payment: {res.text}")
                
        except Exception as e:
            st.error(f"❌ Connection Error: {str(e)}")

# --- 3. CUSTOM CSS STYLING ---
st.markdown(f"""
    <style>
    .main {{ background-color: #0e1117; }}
    .medical-title {{ font-size: 2.8rem; font-weight: 800; color: #ffffff; margin-bottom: 0px; }}
    .medical-subtitle {{ color: #4b9fff; font-size: 1.1rem; font-weight: 500; margin-bottom: 25px; }}
    
    div.stButton > button:first-child[kind="primary"] {{
        background-color: #28a745 !important;
        color: white !important;
        border: none !important;
        height: 3.5em !important;
        width: 100% !important;
        font-weight: bold !important;
    }}
    
    .premium-feature {{
        padding: 10px;
        border-radius: 8px;
        background: #1e2129;
        margin-bottom: 10px;
        border-left: 3px solid #d4af37;
    }}

    .result-card {{
        background-color: #1e2129;
        padding: 25px;
        border-radius: 15px;
        border-left: 5px solid #4b9fff;
        margin-top: 20px;
        line-height: 1.6;
    }}

    .disclaimer-box {{
        background-color: #2c1a1a;
        padding: 15px;
        border-radius: 10px;
        border: 1px solid #dc3545;
        color: #ffbaba;
        font-size: 0.85rem;
        margin-top: 30px;
    }}
    
    .premium-prompt-text {{
        color: #ffffff;
        font-size: 0.95rem;
        margin-bottom: 8px;
    }}
    </style>
    """, unsafe_allow_html=True)

# --- 4. TOP NAVIGATION ---
nav_left, nav_right = st.columns([0.5, 0.5])
with nav_left:
    st.markdown("<h1 class='medical-title'>⚕️ CuraBot</h1>", unsafe_allow_html=True)
    st.markdown("<p class='medical-subtitle'>Intelligent Evidence-Based Health Insights</p>", unsafe_allow_html=True)

with nav_right:
    c1, c2 = st.columns(2)
    with c1:
        if not st.session_state.is_logged_in:
            with st.popover("👤 Sign In", use_container_width=True, disabled=st.session_state.loading):
                tab1, tab2 = st.tabs(["Login", "Create Account"])
                with tab1:
                    lemail = st.text_input("Email", key="l_email")
                    lpass = st.text_input("Password", type="password", key="l_pass")
                    if st.button("Log In", use_container_width=True):
                        if is_valid_email(lemail) and lpass:
                            st.session_state.is_logged_in = True
                            st.session_state.user_email = lemail
                            st.rerun()
                with tab2:
                    remail = st.text_input("Email", key="r_email")
                    rpass = st.text_input("Password", type="password", key="r_pass")
                    rconf = st.text_input("Confirm Password", type="password")
                    if st.button("Register", use_container_width=True):
                        if rpass == rconf and is_valid_email(remail) and len(rpass) >= 8:
                            st.session_state.is_logged_in = True
                            st.session_state.user_email = remail
                            st.rerun()
        else:
            st.markdown("<div style='text-align:center; padding:10px; color:#28a745;'>✅ <b>Logged In</b></div>", unsafe_allow_html=True)

    with c2:
        if not st.session_state.is_premium:
            with st.popover("👑 Go Premium", use_container_width=True):
                render_premium_content()
        else:
            st.markdown("<div style='text-align:center; padding:10px; color:#d4af37;'>⭐ <b>Premium Active</b></div>", unsafe_allow_html=True)

# --- 5. MAIN CONTENT: ASSESSMENT FORM ---
left_panel, right_panel = st.columns([2, 1])

with left_panel:
    with st.expander("🚨 EMERGENCY TRIAGE", expanded=False):
        st.error("Seek immediate care for severe symptoms such as chest pain or breathing difficulties.")

    st.markdown("### 🔍 Symptom Assessment")
    user_input = st.text_area("Describe your symptoms:", height=180, disabled=st.session_state.loading, placeholder="E.g., Severe headache with light sensitivity...")

    p1, p2, p3 = st.columns(3)
    with p1: age = st.selectbox("Age Group", ["Child (0-12)", "Teen (13-17)", "Adult (18-64)", "Senior (65+)"], index=2, disabled=st.session_state.loading)
    with p2: gender = st.selectbox("Gender", ["Male", "Female", "Other"], disabled=st.session_state.loading)
    with p3: severity = st.select_slider("Severity", options=["Mild", "Moderate", "Severe"], disabled=st.session_state.loading)

    b1, b2 = st.columns([0.7, 0.3])
    with b1: analyze_btn = st.button("🚀 Analyze Symptoms", type="primary", disabled=st.session_state.loading, use_container_width=True)
    with b2:
        if st.button("✖ Cancel", use_container_width=True):
            st.session_state.loading = False
            st.rerun()

# --- 6. ANALYSIS LOGIC ---
if analyze_btn and user_input:
    st.session_state.loading = True
    st.rerun()

if st.session_state.loading:
    st.markdown("---")
    st.progress(0)
    st.info("🧬 **Searching Clinical Databases...** (Deep research can take 1-3 minutes)")
    try:
        payload = {"symptoms": f"{user_input} ({age}, {gender}, {severity})"}
        response = requests.get("http://127.0.0.1:8000/analyze", params=payload, timeout=600)
        
        if response.status_code == 200:
            st.session_state.report_content = response.json().get("analysis")
            if st.session_state.is_logged_in:
                st.session_state.history.append(f"{datetime.now().strftime('%H:%M')} - {user_input[:15]}...")
            st.session_state.loading = False
            st.rerun()
        else:
            st.session_state.loading = False
            st.error(f"⚠️ **Server Error:** (Code {response.status_code})")
    except Exception as e:
        st.session_state.loading = False
        st.error(f"❌ Connection Failed: {e}")

# --- 7. SIDEBAR & RESULTS DISPLAY ---
with right_panel:
    st.markdown("### 🔒 Privacy First")
    st.info("Symptoms remain anonymous. No personal data stored.")
    if st.session_state.is_logged_in:
        st.markdown("### 🕒 Recent Activity")
        for h in st.session_state.history[-5:]: st.caption(f"📅 {h}")

if st.session_state.report_content:
    st.markdown("---")
    res_col, pdf_col = st.columns([0.65, 0.35])
    with res_col: st.markdown("### 📋 Clinical Findings")
    with pdf_col:
        if not st.session_state.is_premium:
            st.markdown('<p class="premium-prompt-text">👑 Want to save or share this report?</p>', unsafe_allow_html=True)
            if st.button("✨ Upgrade to Premium", use_container_width=True):
                st.info("💡 Click '👑 Go Premium' button in the top navigation to upgrade.")
        else:
            demo_string = f"{age}, {gender}, {severity} severity"
            pdf_data = generate_pdf(st.session_state.report_content, user_input, demo_string)
            st.download_button(
                label="📄 Download PDF Report",
                data=pdf_data,
                file_name=f"CuraBot_Report_{datetime.now().strftime('%Y%m%d')}.pdf",
                mime="application/pdf",
                use_container_width=True
            )

    st.markdown(f'<div class="result-card">{st.session_state.report_content}</div>', unsafe_allow_html=True)
    st.markdown("""
        <div class="disclaimer-box">
            <b>⚠️ MEDICAL DISCLAIMER:</b> Educational purposes only. This AI-synthesized report is not a substitute 
            for professional medical advice, diagnosis, or treatment. Always consult a physician.
        </div>
    """, unsafe_allow_html=True)